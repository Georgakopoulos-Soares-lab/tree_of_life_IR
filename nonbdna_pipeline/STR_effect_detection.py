from pybedtools import BedTool
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import logging
import csv
import threading
from termcolor import colored
from typing import ClassVar
from attr import field
import attr
from nonbdna_pipeline.logger import Logger
from nonbdna_pipeline.utils import load_bucket
from nonbdna_pipeline.stream_and_merge_bucket import StreamAndMerge

@attr.s 
class STREffectDetector(StreamAndMerge):

    STR_indir: Path = field(factory=lambda: Path("extractions_STR").resolve())
    outdir: Path = field(init=False)
    COVERAGE_FIELDS: ClassVar[list[str]] = ["total_hits", "overlapping_bp", "compartment_length", "coverage"]

    def __attrs_post_init__(self) -> None:
        self.STR_indir = Path(self.STR_indir).resolve()
        if not self.STR_indir.is_dir():
            raise FileNotFoundError(f"Missing STR input directory `{self.STR_indir}`.")
        self.outdir = self.indir.joinpath("STR_effects")
        self.outdir.mkdir(exist_ok=True)
        return super().__attrs_post_init__()

    def process_bucket(self, bucket_id: int, 
                       pattern: str,
                       partition_col: str,
                       min_partition: int = 0,
                       max_partition: int = 9
                       ) -> None:
        fieldnames = [
            "#assembly_accession",
            "partition_col",
            "partition",
            "motif_explained_by_STR",
            "motif_region_bp",
            "motif_explained_by_STR_perc",
            "STR_explained_by_motif",
            "STR_region_bp",
            "STR_explained_by_motif_perc",
            "total_IR_bp",
            "total_STR_bp",
            "overlap_ratio",
        ]
        # Fetch validated files
        infiles = self.load_validated_files_from_log(bucket_id=bucket_id, pattern=pattern)
        total_files_processed = 0
        outfile = self.outdir.joinpath(f"bucket_{bucket_id}_STR_effects.tsv")
        partitions = ["."] + list(range(min_partition, max_partition))

        # Setup Logger
        logger = Logger(total_files=len(infiles))
        logger._setup_logging(bucket_id=bucket_id)
        thread = threading.Thread(target=logger._log_progress, daemon=True)
        thread.start()

        # Start
        fout = open(outfile, mode="w", encoding="utf-8", newline="")
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for file_idx, infile in enumerate(infiles, start=1):
            infile = Path(infile)
            accession_name = StreamAndMerge.extract_name(infile, pattern=pattern)
            accession_id = StreamAndMerge.extract_id(infile)
            # motif_file = self.indir / f"{name}_{pattern}.tsv.gz"
            STR_file = self.STR_indir / f"{accession_name}_STR.tsv.gz"
            if not STR_file.is_file():
                logging.warning(f"Absent STR file for {accession_id} within {self.STR_indir}. (bucket {bucket_id}).")
                continue
            if not infile.is_file():
                logging.warning(f"Absent motif file for {accession_id} within {self.indir}. (bucket {bucket_id}).")
                continue
            motif_df = pd.read_table(infile, 
                                     usecols=["seqID", "start", "end", partition_col], 
                                     dtype={partition_col: int}
                                     )
            STR_df = pd.read_table(STR_file, usecols=["seqID", "start", "end"])
            STR_bed = BedTool.from_dataframe(STR_df).sort().merge()

            # To what extent do STRs "explain" motif regions?
            for partition in partitions:
                if partition != ".":
                    motif_df_temp = motif_df[motif_df[partition_col] == partition].copy()
                else:
                    motif_df_temp = motif_df.copy()

                motif_bed = BedTool.from_dataframe(motif_df_temp).sort().merge()
                motif_cov = pd.read_table(
                    motif_bed.coverage(STR_bed).fn,
                    header=None,
                    names=["seqID", "start", "end"] + STREffectDetector.COVERAGE_FIELDS,
                )
                motif_explained_by_STR = motif_cov["overlapping_bp"].sum()
                motif_region_bp = motif_cov["compartment_length"].sum()
                motif_explained_by_STR_perc = round(1e2 * motif_explained_by_STR / motif_region_bp, 2) if motif_region_bp > 0 else 0.0

                # To what extent do motif regions overlap STRs?
                str_cov = pd.read_table(
                    STR_bed.coverage(motif_bed).fn,
                    header=None,
                    names=["seqID", "start", "end"] + STREffectDetector.COVERAGE_FIELDS,
                )
                STR_explained_by_motif = str_cov["overlapping_bp"].sum()
                STR_region_bp = str_cov["compartment_length"].sum()
                STR_explained_by_motif_perc = round(1e2 * STR_explained_by_motif / STR_region_bp, 2) if STR_region_bp > 0 else 0.0

                writer.writerow({
                    "#assembly_accession": accession_id,
                    "partition_col": partition_col,
                    "partition": partition,
                    "motif_explained_by_STR": motif_explained_by_STR,
                    "motif_region_bp": motif_region_bp,
                    "motif_explained_by_STR_perc": motif_explained_by_STR_perc,
                    "STR_explained_by_motif": STR_explained_by_motif,
                    "STR_region_bp": STR_region_bp,
                    "STR_explained_by_motif_perc": STR_explained_by_motif_perc,
                })
                logger.files_processed += 1
        fout.close()
        print(colored(f"Motif coverage processing for pattern {pattern} in bucket {bucket_id}.", "green"))
        logging.info(f"Process has been completed succesfully (bucket {bucket_id}).")
        return

def main():
    import argparse
    parser = argparse.ArgumentParser(description="STR effect detection per bucket")
    parser.add_argument("schedule", type=str, default="schedule.json")
    parser.add_argument("--indir", "-i", type=str, default="extractions_IR")
    parser.add_argument("--pattern", type=str, default="IR", choices=["IR"])
    parser.add_argument("--STR_indir", type=str, default="extractions_STR")
    parser.add_argument("--min_partition", type=int, default=0)
    parser.add_argument("--max_partition", type=int, default=8)
    parser.add_argument("--outdir", type=str)
    parser.add_argument("--bucket_id", type=int)
    parser.add_argument("--partition", type=str, default="spacer_length")
    
    # # # # # # #
    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)
    indir = Path(args.indir).resolve()
    if not indir.is_dir():
        raise FileNotFoundError(f"Missing motif input directory `{indir}`.")
    STR_indir = Path(args.STR_indir).resolve()
    if not STR_indir.is_dir():
        raise FileNotFoundError(f"Missing STR input directory `{STR_indir}`.")

    STREffectDetector(
        schedule=args.schedule,
        indir=indir,
        STR_indir=STR_indir,
        outdir=outdir,
    ).process_bucket(
        bucket_id=args.bucket_id,
        pattern=args.pattern,
        partition_col=args.partition_col,
        min_partition=args.min_partition,
        max_partition=args.max_partition,
    )

if __name__ == "__main__": main()