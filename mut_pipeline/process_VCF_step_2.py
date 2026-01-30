import csv
import os
import shutil
import tempfile
from collections import defaultdict

import polars as pl
from pybedtools import BedTool, Interval
from tqdm import tqdm

from process_VCF_step_1 import *


class BatchProcessor:
    def __init__(self, location: str, reference: str, outdir: str) -> None:
        self.mutation_types = ["snp", "smalldel", "smallins"]
        self.outdir = Path(outdir).resolve()
        self.outdir.mkdir(exist_ok=True)
        self.location = Path(location).resolve()
        if not self.location.is_dir():
            raise FileNotFoundError(f"Location `{self.location}` doesn't exist.")
        self.infiles = [
            file for file in self.location.glob("*.vcf.gz") if file.is_file()
        ]
        self.total_VCF = len(self.infiles)
        if self.total_VCF == 0:
            raise ValueError("No VCF files were detected.")
        print(colored(f"Total VCF detected: {self.total_VCF}.", "magenta"))
        self.reference = Path(reference).resolve()
        if not self.reference.is_file():
            raise FileNotFoundError(
                f"Reference files `{self.reference}` doesn't exist."
            )
        self.ref: dict = self._read_reference()
        self.extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
        self.genomes_loc = Path("genomes")
        self.genomes_loc.mkdir(exist_ok=True)
        self.download_reference()
        self.fasta_to_species = {value: key for key, value in self.ref.items()}

        self.fasta_files = {
            self.fasta_to_species[self.extract_id(file)]: file
            for file in self.genomes_loc.glob("*.fna")
            if self.extract_id(file) in self.fasta_to_species
        }
        self.total_genomes = len(self.fasta_files)
        print(colored(f"Total fasta files: `{self.total_genomes}`.", "green"))
        self.density_out = Path("VCF_stats").joinpath("density_VCF_chrom_stats.tsv")
        self.density_out.parent.mkdir(exist_ok=True)

    def _read_reference(self) -> dict[int, str]:
        ref = dict()
        with open(self.reference, mode="r", encoding="UTF-8") as f:
            for line in f:
                species_taxid, name = line.strip().split()
                species_taxid = int(species_taxid)
                name = "_".join(name.split("/")[-1].split("_")[:2])
                ref[species_taxid] = name
        return ref

    def process(
        self, MQ: Optional[int] = 50, QUAL: Optional[int] = 60, DP: Optional[int] = None
    ):
        processed_VCF = dict()
        for VCF_in in tqdm(self.infiles):
            species_taxid = int(VCF_in.name.split("_")[0])
            if species_taxid not in self.fasta_files:
                continue
            output_vcf = self.outdir.joinpath(
                Path(VCF_in).stem + f"_{MQ}_{QUAL}_{DP}_monoallelic_filtered.vcf.gz"
            )
            result = VCFProcessor(MQ=MQ, DP=DP, QUAL=QUAL).process(
                VCF_in=VCF_in,
                reference=self.fasta_files[species_taxid],
                output_vcf=output_vcf,
            )
            if result is not None:
                VCFProcessor.run_cmd(f"tabix -p vcf '{output_vcf}'")
            # processed_VCF[result] = None
        # for result in processed_VCF:

    def download_reference(self):
        files = {self.extract_id(file): file for file in self.genomes_loc.glob("*.fna")}
        total = 0
        with tempfile.NamedTemporaryFile("w", delete=True) as tmp_file:
            name = tmp_file.name
            for value in self.ref.values():
                if value not in files:
                    tmp_file.write(value + "\n")
                    total += 1
            tmp_file.flush()
            os.system(f"cat {name}")
            if total > 0:
                os.system(
                    f"datasets download genome accession --inputfile {name} --include genome,gff3"
                )
                VCFProcessor.run_cmd("unzip -n ncbi_dataset")
                VCFProcessor.run_cmd("rm ncbi_dataset.zip")
        ref = Path("ncbi_dataset")
        if ref.is_dir():
            files = [file for file in ref.glob("**/*.fna")]
            for file in files:
                shutil.move(file, self.genomes_loc)
            shutil.rmtree(ref)
        files = {
            self.extract_id(file): file for file in self.genomes_loc.glob("**/*.fna")
        }
        for value in self.ref.values():
            if value not in files:
                print(files)
                raise ValueError(f"Missing file afted downloads {value}.")
        if len(files) == 0:
            raise ValueError(f"No downloads detected.")
        print(colored("CHECKS PASSED!", "green"))
        print(colored("DONE!", "green"))

    def fetch_files(self) -> dict[str, Path]:
        files = {
            int(file.name.split("_")[0]): file
            for file in self.location.glob("*.vcf.gz")
        }
        return files

    def extract_density(self, overwrite: bool = False) -> pl.DataFrame:
        files = self.fetch_files()
        mutation_types = ["snp", "smalldel", "smallins"]
        FIELDS = ["species_taxid", "mutation", "chrom", "size", "group", "events_n"]
        erroneous_records = Path("VCF_stats").joinpath("failed_VCF.txt")
        total_failed = 0
        if self.density_out.is_file() and not overwrite:
            return self
        with (
            erroneous_records.open("w", encoding="UTF-8") as fout1,
            self.density_out.open("w", encoding="UTF-8") as fout2,
        ):
            writer = csv.DictWriter(fout2, delimiter="\t", fieldnames=FIELDS)
            writer.writeheader()
            for species_taxid, VCFin in tqdm(files.items()):
                intervals_by_type = {mut: defaultdict(list) for mut in mutation_types}
                result = VCFProcessor.ensure_vcf_indexed(VCFin)
                if result is None:
                    total_failed += 1
                    fout1.write(f"{species_taxid}\t{VCFin}\n")
                    continue

                genome = Fasta(genome=self.fasta_files[species_taxid])
                genome.extract_gs()
                with pysam.VariantFile(VCFin) as vcf:
                    for record in vcf:
                        variant_type = record.info["TYPE"]
                        if len(variant_type) > 1:
                            raise ValueError(
                                f"Monoallelic assumption not satisfied for `{VCFin}`."
                            )
                        variant_type = variant_type[0]
                        if variant_type not in intervals_by_type:
                            continue
                        intervals_by_type[variant_type][record.chrom].append(
                            # tricky part < we need events not per base
                            # otherwise, place; end=record.end
                            Interval(
                                chrom=record.chrom,
                                start=record.start,
                                end=record.start + 1,
                            )
                        )

                for mutation in mutation_types:
                    if intervals_by_type[mutation]:
                        for chrom in intervals_by_type[mutation]:
                            if intervals_by_type[mutation][chrom]:
                                events_n = (
                                    BedTool(intervals_by_type[mutation][chrom])
                                    .sort()
                                    .merge()
                                    .total_coverage()
                                )
                            else:
                                events_n = 0

                            writer.writerow(
                                {
                                    "species_taxid": int(species_taxid),
                                    "mutation": mutation,
                                    "chrom": chrom,
                                    "size": genome.size[chrom],
                                    "group": genome.plasmids[chrom],
                                    "events_n": events_n,
                                }
                            )
        if total_failed > 0:
            print(colored(f"Failed to process {total_failed} VCF files.", "red"))
        return self

    def annotate_density(self, summary: str):
        df_summary = (
            pl.read_csv(summary, separator="\t", null_values=["na"])
            .select(
                [
                    "#assembly_accession",
                    "species_taxid",
                    "genome_size",
                    "genome_size_ungapped",
                    "gc_percent",
                    "organism_name",
                    "total_gene_count",
                    "non_coding_gene_count",
                    "protein_coding_gene_count",
                ]
            )
            .filter(
                pl.col("#assembly_accession").is_in(set(self.fasta_to_species.keys()))
            )
            .with_columns(
                pl.col("species_taxid").cast(pl.UInt32).alias("species_taxid"),
                pl.col("genome_size").cast(pl.UInt32).alias("genome_size"),
            )
        )

        (
            pl.read_csv(self.density_out, separator="\t")
            .with_columns(
                pl.col("species_taxid").cast(pl.UInt32).alias("species_taxid")
            )
            .group_by(["species_taxid", "mutation"], maintain_order=True)
            .agg(pl.col("events_n").sum(), pl.col("size").sum())
            .join(df_summary, on="species_taxid", how="left")
            .with_columns(density=pl.col("events_n") * 1e3 / pl.col("genome_size"))
            .write_csv(
                Path("VCF_stats").joinpath(
                    "density_VCF_chrom_stats_additional.species_level.tsv"
                ),
                separator="\t",
                include_header=True,
                float_precision=4,
            )
        )

        (
            pl.read_csv(self.density_out, separator="\t")
            .with_columns(
                pl.col("species_taxid").cast(pl.UInt32).alias("species_taxid")
            )
            .group_by(["species_taxid", "mutation", "group"], maintain_order=True)
            .agg(pl.col("events_n").sum(), pl.col("size").sum())
            .join(df_summary, on="species_taxid", how="left")
            .with_columns(density=pl.col("events_n") * 1e3 / pl.col("size"))
            .write_csv(
                Path("VCF_stats").joinpath(
                    "density_VCF_chrom_stats_additional.group_level.tsv"
                ),
                separator="\t",
                include_header=True,
                float_precision=4,
            )
        )

        # return df

def main():
    import argparse
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("--location", type=str, default="../pangenomes/results_snp")
    parser.add_argument("--reference", type=str, default="reference_files.txt")
    parser.add_argument("--MQ", type=int, default=50)
    parser.add_argument("--QUAL", type=int, default=60)
    parser.add_argument("--DP", type=int)
    parser.add_argument("--outdir", type=str, default="processed_VCF")
    parser.add_argument("--overwrite", type=int, default=0, choices=[0, 1])
    args = parser.parse_args()

    print(colored(f"Maximum MQ: {args.MQ}.", "magenta"))
    print(colored(f"Maximum QUAL: {args.QUAL}.", "magenta"))
    print(colored(f"Maximum DP: {args.DP}.", "magenta"))
    summary = Path("assembly_summary_with_tree.csv.gz")
    processor = BatchProcessor(
        location=args.location, reference=args.reference, outdir=args.outdir
    )
    # processor.process(MQ=args.MQ, QUAL=args.QUAL, DP=args.DP)
    processor.extract_density(overwrite=args.overwrite)
    processor.annotate_density(summary=summary)


if __name__ == "__main__":
    main()
