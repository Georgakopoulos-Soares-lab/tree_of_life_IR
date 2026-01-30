import os
import shutil
from pathlib import Path
import pandas as pd
import polars as pl
from pybedtools import BedTool
from tqdm import tqdm
from gff_utils import GFFExtractor
from process_VCF_step_1 import Fasta

extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])

def annotate_motifs(data: str, motifs: str, outdir: str):
    data = Path(data)
    if not data.is_dir():
        raise ValueError(f"Directory `{data}` does not exist!")
    motifs = Path(motifs)
    if not motifs.is_dir():
        raise ValueError(f"Directory `{motifs}` does not exist!")
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    gff_files = {
        infile.parent.name: infile
        for infile in Path("ncbi_dataset/data").glob("**/*.gff")
        if infile.is_file()
    }
    fasta_files = {
        extract_id(infile): infile
        for infile in data.glob("**/*.fna")
        if infile.is_file()
    }
    motif_files = {
        extract_id(infile): infile for infile in motifs.glob("**/*.processed.tsv")
    }
    GFF_FIELDS = ["chrom", "ann_start", "ann_end", "strand", "gene_length", "biotype"]
    e = GFFExtractor(compartments=["gene"])
    TERMINATOR_KB = 50
    for accession_id, infile in tqdm(gff_files.items()):
        if accession_id not in motif_files:
            print(f"Skipping {accession_id}.")
            continue
        paired_fasta = fasta_files[accession_id].name.split(".fna")[0]
        dest_gff = f"{outdir}/{paired_fasta}.gff"
        shutil.copy(infile, dest_gff)

        fasta = Fasta(fasta_files[accession_id]).extract_gs()
        gff_table = e.read_gff(dest_gff, overload_biotype=False, parse_ID=True)
        gff_table = gff_table.with_columns(
            pl.col("seqID")
            .map_elements(lambda y: fasta.size[y], return_dtype=pl.Int32)
            .alias("region_end"),
            pl.lit(0).alias("region_start"),
        )
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
        gff_bed = BedTool.from_dataframe(terminators_df.to_pandas()).sort()
        motif_infile = motif_files[accession_id]
        motif_df = pl.read_csv(motif_infile, separator="\t")
        motif_bed = BedTool.from_dataframe(motif_df.to_pandas()).sort()
        terminator_motifs_bed = motif_bed.intersect(gff_bed, u=True, f=0.7)
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
            f"{outdir}/{paired_fasta}.processed.annotated.csv",
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
    parser.add_argument("--data", type=str, default="genomes")
    parser.add_argument("--motifs", type=str, default="motifs_IR")
    parser.add_argument("--outdir", type=str, default="annotated")
    args = parser.parse_args()
    # download_gff(data=args.data)
    annotate_motifs(data=args.data, outdir=args.outdir, motifs=args.motifs)

if __name__ == "__main__":
    main()
