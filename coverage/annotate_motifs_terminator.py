import csv
import gzip
import lzma
import os
from pathlib import Path
from typing import Iterable

import attr
import pandas as pd
import polars as pl
from attr import field
from Bio.SeqIO.FastaIO import SimpleFastaParser
from gff_utils import GFFExtractor
from pybedtools import BedTool
from tqdm import tqdm

from utils import load_bucket

extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])


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


def annotate_motifs(
    genome_files: dict[str, str],
    motif_mapping: dict[str, str],
    gff_files: dict[str, str],
    outdir: str,
):
    outdir = Path(outdir).resolve()
    outdir.mkdir(exist_ok=True)
    # GFF_FIELDS = ["chrom", "ann_start", "ann_end", "strand", "gene_length", "biotype"]
    e = GFFExtractor(compartments=["gene"])
    TERMINATOR_KB = 50
    for gff_file in tqdm(gff_files, leave=True):
        accession_id = extract_id(gff_file)
        if accession_id not in motif_mapping:
            print(f"Skipping {accession_id}.")
            continue

        # fetch genome
        infile = Path(genome_files[accession_id])
        accession_name = infile.name.split(".fna")[0]
        fasta = Fasta(infile).extract_gs()

        # load gff
        gff_table = e.read_gff(gff_file, overload_biotype=False, parse_ID=True)
        gff_table = gff_table.with_columns(
            pl.col("seqID")
            .map_elements(lambda y: fasta.size[y], return_dtype=pl.Int32)
            .alias("region_end"),
            pl.lit(0).alias("region_start"),
        )
        # # #
        terminators_df = e.parse_terminators(
            gff_table,
            use_names=True,
            terminator_kb=TERMINATOR_KB,
        ).select(
            [
                "seqID",
                "start",
                "end",
                "strand",
                "gene_length",
                "biotype",
            ]
        )
        # rRNA = terminators_df.filter(pl.col("biotype") == "rRNA")
        # rRNA_bed = BedTool.from_dataframe(rRNA.to_pandas()).sort()
        gff_bed = BedTool.from_dataframe(terminators_df.to_pandas()).sort()

        # fetch & Load Motifs
        motif_infile = motif_mapping[accession_id]
        motif_df = pl.read_csv(motif_infile, separator="\t")
        motif_bed = BedTool.from_dataframe(motif_df.to_pandas()).sort()
        terminator_motifs_bed = motif_bed.intersect(gff_bed, u=True, f=0.7)
        # rRNA_terminators_bed = motif_bed.intersect(rRNA_bed, u=True, f=0.7)
        # terminator_motifs = pl.read_csv(
        #     rRNA_terminators_bed.fn,
        #     has_header=False,
        #     separator="\t",
        #     new_columns=list(motif_df.columns),
        # ).with_columns(biotype=pl.lit("rRNA"))

        terminator_motifs = pl.read_csv(
            terminator_motifs_bed.fn,
            has_header=False,
            separator="\t",
            new_columns=list(motif_df.columns),
        ).with_columns(compartment=pl.lit("Terminator"))
        # non_motifs = pl.read_csv(
        #    non_motifs_bed.fn,
        #    has_header=False,
        #    separator="\t",
        #    new_columns=list(motif_df.columns),
        # ).with_columns(compartment=pl.lit("Other"))
        size = motif_df.shape[0]
        motif_df = motif_df.join(
            terminator_motifs,
            on=list(motif_df.columns),
            how="left",
        ).with_columns(
            compartment=pl.when(pl.col("compartment").is_null())
            .then(pl.lit("Other"))
            .otherwise(pl.col("compartment"))
        )
        assert size == motif_df.shape[0], "Invalid size!"
        motif_df.write_csv(
            f"{outdir}/{accession_name}.processed.annotated.csv",
            separator="\t",
            include_header=True,
        )

    return


def download_gff(data: str):
    data = Path(data)
    if not data.is_dir():
        raise ValueError(f"{data} is not a directory!")
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    fasta_files = {
        extract_id(infile): infile for infile in data.glob("*.fna") if infile.is_file()
    }
    inputfile = "inputfile.txt"
    with open(inputfile, mode="w", encoding="UTF-8") as f:
        for accession_id in fasta_files:
            f.write(f"{accession_id}\n")
    command = (
        f"datasets download genome accession --inputfile {inputfile} --include gff3"
    )
    os.system(command)
    os.system("unzip -n ncbi_dataset")
    return


def main():
    import argparse

    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("schedule", type=str)
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--design", type=str)
    parser.add_argument("--outdir", type=str, default="annotated_motifs")
    args = parser.parse_args()
    # download_gff(data=args.data)
    #
    motif_mapping = dict()
    gff_files = dict()
    with open(args.design, mode="r", encoding="UTF-8") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            gff_files[extract_id(row["#assembly_accession"])] = row["gff_file"]
            motif_mapping[extract_id(row["#assembly_accession"])] = row["extraction"]

    files = {
        extract_id(infile): infile
        for infile in load_bucket(args.schedule, args.bucket_id)
    }
    annotate_motifs(
        genome_files=files,
        motif_mapping=motif_mapping,
        gff_files=gff_files,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    main()
