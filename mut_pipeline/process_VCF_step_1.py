import gzip
import lzma
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from subprocess import run
from typing import Generator, Iterable, Iterator, Optional
import attr
import pandas as pd
import polars as pl
import pysam
from attr import field
from Bio.SeqIO.FastaIO import SimpleFastaParser
from termcolor import colored

class MultiallelicRecordException(Exception):
    pass

@attr.s
class Fasta:
    genome: str = field()
    genome_size: int = field(init=False)
    size: dict[str, int] = field(init=False)
    plasmids: dict[str, str] = field(init=False)

    def __attrs_post_init__(self):
        self.genome = Path(self.genome).resolve()
        if not self.genome.is_file():
            raise FileNotFoundError(f"Fasta file `{self.genome}` does not exist.")
        self.size = dict()

    @staticmethod
    def parse_fasta(genome: str) -> Iterable[tuple[str, str]]:
        if Path(genome).name.endswith(".xz"):
            fin = lzma.open(genome, "r")
        elif Path(genome).name.endswith(".gz"):
            fin = gzip.open(genome, "rt")
        else:
            fin = open(genome, "r", encoding="UTF-8")
        for record in SimpleFastaParser(fin):
            seqID, seq = record
            yield seqID, seq
        fin.close()

    def _parse_fasta(self) -> Iterable[tuple[str, str]]:
        for seqID, seq in Fasta.parse_fasta(self.genome):
            yield seqID, seq

    def extract_gs(self):
        self.genome_size = 0
        self.plasmids = dict()
        for seqID, seq in self._parse_fasta():
            sequence_length = len(seq)
            self.genome_size += sequence_length
            seqID_name = seqID.split(" ")[0]
            self.size[seqID_name] = sequence_length
            if "plasmid" in seqID.lower():
                self.plasmids[seqID_name] = "plasmid"
            elif "chromosome" in seqID.lower():
                self.plasmids[seqID_name] = "chromosome"
            else:
                self.plasmids[seqID_name] = "unknown"
        return self


def determine_TYPE(ref_len: int, alt_len: int, MAX_SMALL_INDEL_SIZE: int = 50) -> str:
    if ref_len == alt_len == 1:
        return "snp"
    elif ref_len > alt_len and ref_len < MAX_SMALL_INDEL_SIZE:
        return "smalldel"
    elif ref_len > alt_len:
        return "del"
    elif alt_len > ref_len and alt_len < MAX_SMALL_INDEL_SIZE:
        return "smallins"
    elif alt_len > ref_len:
        return "ins"
    return "complex"


@attr.s(slots=True)
class VCFProcessor:
    MQ: Optional[int] = field(converter=attr.converters.optional(int), default=None)
    QUAL: Optional[int] = field(
        converter=attr.converters.optional(int),
        default=None,
    )
    DP: Optional[int] = field(converter=attr.converters.optional(int), default=None)

    @staticmethod
    def run_cmd(*args, **kwargs):
        try:
            return run(*args, **kwargs, shell=True, check=True)
        except Exception as err:
            return

    def __attrs_post_init__(self):
        pass

    def perform_quality_control(
        self, VCF_in: str, reference: str, annotate: bool = True
    ) -> str:
        if not shutil.which("bcftools"):
            raise RuntimeError(
                "bcftools is not available in PATH. Please install bcftools."
            )
        with tempfile.NamedTemporaryFile(
            delete=False, mode="wb", suffix=".normalized.vcf.gz"
        ) as tmp_file:
            conditions = []
            if self.QUAL:
                conditions.append(f"QUAL>={self.QUAL}")
            if self.MQ:
                conditions.append(f"INFO/MQ>={self.MQ}")
            if self.DP:
                conditions.append(f"INFO/DP>={self.DP}")
            condition = " && ".join(conditions)
            if condition:
                command = (
                    f"bcftools norm -m-any -f {reference} {VCF_in}"
                    f"| bcftools filter -i '{condition}'"
                    f"| bcftools +fill-tags -- -t AF,AC,AN,TYPE"
                    f"| bcftools view -Oz -o {tmp_file.name}"
                )
            else:
                command = (
                    f"bcftools norm -m-any -f {reference} {VCF_in}"
                    f"| bcftools +fill-tags -- -t AF,AC,AN,TYPE"
                    f"| bcftools view -Oz -o {tmp_file.name}"
                )
                # if annotate:
                #    command += "bcftools +fill-tags - -- -t AF,AC,AN"
            result = VCFProcessor.run_cmd(command)
            if result is None:
                return
            if not Path(tmp_file.name).is_file():
                raise FileNotFoundError(
                    f"Failed to produce a normalized VCF file from reference `{reference}`."
                )
            self.ensure_vcf_indexed(tmp_file.name)
        return tmp_file.name

    @staticmethod
    def ensure_vcf_indexed(vcf_path: str):
        """Ensure VCF file has tabix index"""
        vcf_path = Path(vcf_path)
        index_path = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
        if index_path.exists():
            os.remove(index_path)
        if not index_path.exists():
            print(f"Creating index for {vcf_path}")
            try:
                pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
            except OSError as e:
                print(colored(f"Failed to index {vcf_path}: {e}.", "red"))
                return
        return str(vcf_path)

    def process(
        self,
        VCF_in: str,
        output_vcf: str,
        reference: str,
        quality_control: bool = True,
        max_indel_size: int = 10,
        MAX_SMALL_INDEL_SIZE: int = 50,
    ):
        output_vcf = Path(output_vcf).resolve()
        if quality_control:
            VCF_in = self.perform_quality_control(VCF_in, reference)
            if VCF_in is None:
                return
        vcf_in = pysam.VariantFile(VCF_in)
        vcf_out = pysam.VariantFile(output_vcf, "wz", header=vcf_in.header)

        for record in vcf_in:
            ref_len = len(record.ref)
            total_alts = len(record.alts)
            if total_alts > 1:
                raise MultiallelicRecordException()
            alt_len = len(str(record.alts[0]))
            # size_diff = abs(alt_len - ref_len)
            variant_type = determine_TYPE(ref_len, alt_len, MAX_SMALL_INDEL_SIZE)
            record.info["TYPE"] = variant_type
            vcf_out.write(record)
        # VCFProcessor.run_cmd(f"tabix {vcf_out}")
        if not output_vcf.is_file():
            raise FileNotFoundError(f"Failed to create new VCF `{vcf_out}`.")
        # self.ensure_vcf_indexed(output_vcf)
        return output_vcf

def main():
    import argparse
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("--VCF_in", type=str, default=".")
    parser.add_argument("--reference", type=str)
    parser.add_argument("--MQ", type=int, default=50)
    parser.add_argument("--QUAL", type=int, default=60)
    parser.add_argument("--DP", type=int)
    parser.add_argument("--outdir", type=str, default="processed_VCF")
    args = parser.parse_args()

    print(colored(f"Maximum MQ: {args.MQ}.", "magenta"))
    print(colored(f"Maximum QUAL: {args.QUAL}.", "magenta"))
    print(colored(f"Maximum DP: {args.DP}.", "magenta"))
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)
    output_vcf = outdir.joinpath(
        Path(args.VCF_in).stem
        + f"_{args.MQ}_{args.QUAL}_{args.DP}_single_allele_filtered.vcf.gz"
    )
    VCFProcessor(MQ=args.MQ, DP=args.DQ, QUAL=args.QUAL).process(
        VCF_in=args.VCF_in, reference=args.reference, output_vcf=output_vcf
    )


if __name__ == "__main__":
    main()