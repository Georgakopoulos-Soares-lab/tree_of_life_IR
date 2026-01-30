# Null Hypothesis Definition
#
import itertools
import os
import pickle
import shutil
from abc import abstractmethod
from typing import ClassVar, Optional

from Bio.Seq import Seq
from process_VCF_step_2 import *

@attr.s(slots=True, frozen=True)
class Prediction:
    base_rate: float = field()
    event_rate: dict[str, float] = field()

class NotFittedError(Exception):
    def __init__(
        self,
        message="Model must be fitted before making predictions. Call fit() method first.",
    ):
        super().__init__(message)

class PolymorphismModel:
    def __init__(self):
        self._data: Optional[pl.DataFrame] = None

    @abstractmethod
    def fit(self, VCF_in: str, fasta: str) -> "PolymorphismModel":
        raise NotImplementedError()

    @abstractmethod
    def predict(self, sequence: str) -> float:
        raise NotImplementedError()

    def saveas(self, outdir: Optional[Path | str]) -> None:
        if outdir is None:
            outdir = Path("models")
            outdir.mkdir(exist_ok=True)
        else:
            outdir = Path(outdir)
        if not outdir.is_dir():
            raise FileNotFoundError(f"Invalid output directory `{outdir}`.")


class NaiveModel(PolymorphismModel):
    def __init__(self):
        super().__init__()
        self.base_rate: Optional[float] = None

    def fit(self, VCF_in: str, fasta: str) -> "NaiveModel":
        intervals = []
        vcf = pysam.VariantFile(VCF_in)
        gs = Fasta(fasta).extract_gs()
        for record in vcf:
            for idx, _ in enumerate(record.alts, start=0):
                if record.info["TYPE"][idx] != "snp":
                    continue
                intervals.append(
                    Interval(
                        chrom=record.chrom,
                        start=record.start,
                        end=record.stop,
                        score=str(round(record.info["AF"][idx], 4)),
                    )
                )
        vcf.close()
        data = BedTool(intervals)
        base_rate = data.sort().merge().count() / gs.genome_size
        if base_rate == 0:
            raise ValueError(f"Empty vcf file?")
        self._data = data
        self.base_rate = base_rate
        return self

    def predict(self, sequence: int) -> float:
        return len(sequence) * self.base_rate

    def saveas(self, outdir: Optional[str | Path]) -> None:
        super().saveas(outdir)
        with open(outdir / f"model_naive_{self.offset}.pkl", "wb") as fout:
            pickle.dump(self, fout)


class TrinucleotideModel(PolymorphismModel):
    FIELDS: ClassVar[list[str]] = ["seqID", "start", "end", "transition"]

    def __init__(self, offset: int = 1, threads: int = 8):
        self.nucleotides = "agct".upper()
        self.offset = offset
        self.context_length = self.offset * 2 + 1
        self.threads = threads
        self.event_rate: Optional[dict[str, float]] = None
        self.mutation_counts: Optional[dict[str, float]] = None
        self.gw_occurrences: Optional[dict[str, int]] = None
        self.canonical = self.collect(self.context_length)
        self.spectrum = {
            "A>T": "AT>TA",
            "T>A": "AT>TA",
            "A>G": "AT>GC",
            "T>C": "AT>GC",
            "A>C": "AT>CG",
            "T>G": "AT>CG",
            "G>A": "GC>AT",
            "C>T": "GC>AT",
            "G>C": "GC>CG",
            "C>G": "GC>CG",
            "G>T": "GC>TA",
            "C>A": "GC>TA",
        }
        self.spectrum_ids = {
            "C>T": 1,
            "G>A": 1,
            "T>C": 2,
            "A>G": 2,
            "C>A": 3,
            "G>T": 3,
            "C>G": 4,
            "G>C": 4,
            "T>A": 5,
            "A>T": 5,
            "T>G": 6,
            "A>C": 6,
        }

    def collect(self, kmer_length: int) -> dict[str, str]:
        canonical = dict()
        trinucleotides = map(
            lambda x: "".join(x),
            itertools.product(self.nucleotides, repeat=kmer_length),
        )
        for tri in trinucleotides:
            tri_rev = Seq(tri).reverse_complement()._data.decode("utf-8")
            mini = min(tri, tri_rev)
            canonical[tri] = mini
            canonical[tri_rev] = mini
        return canonical

    def saveas(
        self, outdir: Optional[Path | str] = None, species_taxid: Optional[int] = None
    ) -> None:
        super().saveas(outdir)
        # if species_taxid is None:
        #    outfile = outdir / f"model_trinucleotide_context_{self.offset}.pkl"
        # else:
        #    outfile = f"{outdir}/model_trinucleotide_context_{self.offset}_{species_taxid}.pkl"
        # with open(outfile, "wb") as fout:
        #    pickle.dump(self, fout)
        if species_taxid is None:
            outfile = outdir / f"trinucleotide_model_event_rate_{self.offset}.pkl"
        else:
            outfile = f"{outdir}/trinucleotide_model_event_rate_{self.offset}_{species_taxid}.pkl"
        # SAVE events
        with open(outfile, "wb") as fout:
            # pickle.dump(self, fout)
            pickle.dump(self.event_rate | {"baseline": self.base_rate}, fout)

    def count_gw_occurrences(
        self, fasta: str, kmer_length: int = 3, mem: str = "12M", threads: int = 8
    ) -> dict[str, int]:
        if shutil.which("jellyfish") is None:
            raise ValueError(
                "jellyfish is required. Please install using micromamba or conda."
            )
        occurrences = defaultdict(int)
        tmp_file1 = tempfile.NamedTemporaryFile("w", suffix=".jf", delete=False)
        tmp_file2 = tempfile.NamedTemporaryFile("w", suffix=".counts.jf", delete=False)
        tmp_file1.close()
        tmp_file2.close()
        try:
            VCFProcessor.run_cmd(
                f"jellyfish count -m {kmer_length} -t {threads} -s {mem} -o {tmp_file1.name} {fasta}"
            )
            VCFProcessor.run_cmd(
                f"jellyfish dump -c {tmp_file1.name} > {tmp_file2.name}"
            )
            canonical = self.collect(kmer_length=kmer_length)
            with open(tmp_file2.name, mode="r", encoding="UTF-8") as fin:
                for line in fin:
                    line = line.strip().split(" ")
                    occurrences[canonical[line[0]]] += int(line[1])
        finally:
            os.unlink(tmp_file1.name)
            os.unlink(tmp_file2.name)
        return occurrences

    def predict(self, sequence: str) -> float:
        if self.event_rate is None:
            raise NotFittedError()
        content = defaultdict(int)
        for i in range(self.offset, len(sequence) - self.offset):
            content[self.canonical[sequence[i - self.offset : i + self.offset]]] += 1
        event_burden = sum(content[kmer] * self.event_rate[kmer] for kmer in content)
        return event_burden

    def fit(self, VCF_in: str, fasta: str) -> "TrinucleotideModel":
        VCF_with_context: pl.DataFrame = self.extend_context(fasta=fasta, VCF_in=VCF_in)
        self._data = VCF_with_context
        self.gw_occurrences = self.count_gw_occurrences(
            fasta=fasta, kmer_length=self.context_length, threads=self.threads
        )
        # spectrum raw
        self.spectrum_rate = VCF_with_context.group_by(["transition"]).agg(
            total=pl.len()
        )
        # event rate with context and spectrum
        mutation_counts = (
            VCF_with_context.group_by(["context", "transition"])
            .agg(total=pl.len())
            .with_columns(transition=pl.col("context") + ";" + pl.col("transition"))
        )
        mutation_counts = dict(
            zip(mutation_counts["transition"], mutation_counts["total"])
        )
        self.mutation_counts = {
            transition: occurrences / self.gw_occurrences[transition.split(";")[0]]
            for transition, occurrences in mutation_counts.items()
        }
        # event rate with context
        event_rate = (
            VCF_with_context.unique(["seqID", "start", "end", "context"])
            .group_by("context")
            .agg(occurrences=pl.len())
        )
        event_rate = dict(zip(event_rate["context"], event_rate["occurrences"]))
        event_rate = {
            context: occurrences / self.gw_occurrences[context]
            for context, occurrences in event_rate.items()
        }
        self.event_rate = event_rate
        return self

    def extend_context(self, fasta: str, VCF_in: str) -> pl.DataFrame:
        # Regenerate FASTA index if needed to avoid warnings
        fasta_path = Path(fasta)
        index_path = Path(f"{fasta}.fai")

        if not index_path.is_file() or (
            index_path.is_file()
            and fasta_path.stat().st_mtime > index_path.stat().st_mtime
        ):
            VCFProcessor.run_cmd(f"samtools faidx {fasta}")

        intervals = []
        vcf = pysam.VariantFile(VCF_in)
        gs = Fasta(fasta).extract_gs()
        for record in vcf:
            ref = record.ref
            for idx, alt in enumerate(record.alts, start=0):
                if record.info["TYPE"][idx] != "snp":
                    continue
                transition = self.spectrum[f"{ref}>{alt}"]
                intervals.append(
                    Interval(
                        chrom=record.chrom,
                        start=max(0, record.start - self.offset),
                        end=min(gs.size[record.chrom], record.stop + self.offset),
                        name=transition,
                        score=str(round(record.info["AF"][idx], 4)),
                    )
                )
        vcf.close()
        extended_offset_context = BedTool(intervals)
        # extended_offset_context = extended_loci_bed.sort().merge()
        base_rate = extended_offset_context.sort().merge().count() / gs.genome_size
        if base_rate == 0:
            raise ValueError(f"Empty vcf file?")
        self.base_rate = base_rate
        extended_offset_context = extended_offset_context.sequence(
            fi=fasta, name=True, tab=True
        )
        canonical = self.collect(kmer_length=self.context_length)
        VCF_with_context = []
        with open(extended_offset_context.seqfn, mode="r", encoding="UTF-8") as f:
            for line in f:
                line = line.strip().split()
                data = line[0]
                spectrum, data = data.split("::")
                seqID, coords = data.split(":")
                start, end = list(map(int, coords.split("-")))
                if end - start < self.context_length:
                    continue
                context = canonical[line[1]]
                VCF_with_context.append(
                    {
                        "seqID": seqID,
                        "start": start,
                        "end": end,
                        "transition": spectrum,
                        "context": context,
                    }
                )
        VCF_with_context = pl.DataFrame(VCF_with_context)
        VCF_with_context = VCF_with_context.join(
            pl.read_csv(
                extended_offset_context.fn,
                has_header=False,
                separator="\t",
                new_columns=TrinucleotideModel.FIELDS + ["AF"],
                columns=list(range(5)),
            ),
            on=TrinucleotideModel.FIELDS,
            how="left",
        )
        return VCF_with_context


def main():
    import argparse

    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("--vcf", type=str)
    parser.add_argument("--fasta", type=str)
    parser.add_argument("-o", type=str, default="models")
    parser.add_argument("--offset", type=int, default=1)
    parser.add_argument("-t", "--threads", type=int, default=8)
    args = parser.parse_args()
    vcf = args.vcf
    fasta = args.fasta
    outdir = args.o
    if outdir is None:
        outdir = Path().cwd()
    else:
        outdir = Path(outdir).resolve()
        outdir.mkdir(exist_ok=True)
    model = TrinucleotideModel(offset=args.offset, threads=args.threads)
    model.fit(VCF_in=args.vcf, fasta=args.fasta)
    model.saveas(outdir=outdir)


if __name__ == "__main__":
    main()
