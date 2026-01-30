import time
import polars as pl 
from pathlib import Path 
from typing import ClassVar
from termcolor import colored
from typing import Optional 
import pandas as pd
import numpy as np
import pyranges as pr
import threading 
import gzip
import logging
from typing import Optional, ClassVar 
import attr
from attr import field
from pybedtools import BedTool
from nonbdna_pipeline.stream_and_merge_bucket import StreamAndMerge 
from nonbdna_pipeline.pwm_density import PWMExtractor 
from nonbdna_pipeline.tss_tes_processing import * 

@attr.s 
class GFFMotifCoverageProcessor(StreamAndMerge):
    files_processed: int = field(init=False, default=0)
    log_dir: Path = field(init=False)
    total_files: int = field(init=False, default=0)
    outdir: Path = field(init=False)
    polarities: list[str] = field(init=False, factory=lambda: ["Template", "Non-Template"])
    biotypes: list[str] = field(init=False, factory=lambda: ["protein_coding", "non_coding", "."])
    valid_compartments: set[str] = field(factory=lambda: {"gene", "exon", "CDS", "five_prime_UTR", "three_prime_UTR"})
    COVERAGE_FIELDS: ClassVar[list[str]] = ["total_hits", "total_bases", "compartment_length", "coverage"]
    LOG_INTERVAL: int = 240
    FIELDS: list[str] = ["#assembly_accession",
                         "pattern",
                         "compartment",
                         "biotype",
                         "total_bases",
                         "compartment_length",
                         "total_compartments",
                         "pct_at_least_one",
                         ]
    def __attrs_post_init__(self) -> None:
        super().__attrs_post_init__()
        self.outdir = self.indir.joinpath("gff_motif_coverage")
        self.outdir.mkdir(exist_ok=True, parents=True)
        self.log_dir = self.log_indir.joinpath("gff_motif_coverage_logs")
        self.log_dir.mkdir(exist_ok=True, parents=False)
        return

    @staticmethod 
    def _get_min_max_partition(pattern: str) -> tuple[Optional[int], Optional[int]]:
        if pattern == "IR" or pattern == "DR":
            min_partition, max_partition = 0, 9
        elif pattern == "MR" or pattern == "HDNA" or pattern == "H-DNA" or pattern == "GT":
            min_partition, max_partition = 0, 8
        elif pattern == "STR":
            min_partition, max_partition = 1, 10
        else:
            raise ValueError(f"Invalid pattern `{pattern}`.")
        return min_partition, max_partition

    class _TrackProgress:
        def __init__(self, bucket_id: int, 
                           total_records: int,
                           log_indir: Path) -> None:
                        
            self.track = 0
            self.total_records = total_records
            self.bucket_id = bucket_id
            self.log_indir = log_indir
            self._setup_logging()
        def _setup_logging(self) -> None:
            DATE = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s [%(levelname)s] %(message)s",
                filename=self.log_indir.joinpath(f"gff_motif_coverage_{DATE}_{self.bucket_id}.log"),
            )
        def _log_progress(self) -> None:
            while True:
                progress = self.track * 1e2 / self.total_records if self.total_records > 0 else 0.0
                logging.info(f"Current progress for bucket `{self.bucket_id}`: {progress:.2f}.")
                time.sleep(GFFMotifCoverageProcessor.LOG_INTERVAL)

    def process_bucket(self, bucket_id: int, 
                            pattern: str, 
                            gff_indir: Path,
                            use_biotype: bool = False,
                            pseudogenes_to_genes: bool = True,
                            partition_col: Optional[str] = None,
                            assembly_summary: Optional[str] = None,
                            polarity: bool = False
                            ) -> None:
        gff_indir = Path(gff_indir).resolve()
        # LOGGING
        tracker = GFFMotifCoverageProcessor._TrackProgress(
            bucket_id=bucket_id,
            total_records=0,
            log_indir=self.log_dir,
        )
        logger = threading.Thread(target=tracker._log_progress, 
                                  daemon=True, 
                                  name="LoggerDaemon")
        logger.start()

        # LOAD FILES
        FIELDS = ["seqID", "start", "end", "compartment"] + GFFMotifCoverageProcessor.COVERAGE_FIELDS
        EXPECTED_FIELDS = ["#assembly_accession", 
                    "biotype", 
                    "pattern", 
                    "compartment", 
                    "total_bases", 
                    "compartment_length", 
                    "total_coverage",
                    "min_coverage",
                    "max_coverage",
                    "avg_coverage",
                    "median_coverage",
                    "std_coverage",
                    "total_compartments", 
                    "pct_at_least_one"]

        infiles = self.load_validated_files_from_log(bucket_id=bucket_id, pattern=pattern)
        tracker.total_records = len(infiles)
        outfile = self.outdir.joinpath(f"gff_motif_coverage_{pattern}_bucket_{bucket_id}.tsv.gz")
        logging.info(f"Initializing extraction process for pattern {pattern} (bucket {bucket_id}).")
        fin = gzip.open(outfile, mode="wt")
        wrote_header = False

        for file_idx, infile in enumerate(infiles, start=1):
            tracker.track += 1
            accession_name = StreamAndMerge.extract_name(infile, pattern)
            accession_id = StreamAndMerge.extract_id(infile)
            extraction_file = self.indir.joinpath(accession_name + f"_{pattern}.processed.tsv")
            gff_file = gff_indir.joinpath(accession_name + ".gff")
            if not gff_file.is_file():
                gff_file = gff_file.with_suffix(".gff.gz")
                if not gff_file.is_file():
                    continue 
            gff_df = read_gff(gff_file, 
                              end_to_one_bp=False, 
                              pseudogenes_to_genes=pseudogenes_to_genes,
                              parse_biotype=use_biotype,
                              filter_on=self.valid_compartments
            ).rename(columns={"start": "Start",
                              "end": "End"})
            if gff_df.shape[0] == 0:
                logging.warning(f"No features found in GFF file `{gff_file}`. Skipping accession `{accession_id}`.")
                continue 

            df = pd.read_table(extraction_file)
            if partition_col:
                # partitions += list(range(min_partition, max_partition))
                raise NotImplementedError("Partitioning not yet implemented.")
            if df.shape[0] == 0:
                logging.warning(f"No motifs found in file `{extraction_file}`.")
                # continue 
            df = df.rename(columns={"seqID": "Chromosome", 
                                    "start": "Start", 
                                    "end": "End", 
                                    "strand": "Strand"})
            # if polarity:
            #   df = calculate_strand_polarity(df, pattern=pattern)
            df_gr = pr.PyRanges(df).merge(strand=False).as_df()
            df_bed = BedTool.from_dataframe(df_gr).sort()
            biotypes = self.biotypes if use_biotype else ["."]
            # Why we need to loop
            # We are attempting to calculate the proportion of:
            # - Total Genic Region covered by motif base pairs 
            # Thus, we want this number to be merged, otherwise we will be counting overlapping regions multiple times
            # Instead of (X/Y) we will computing (X+ε/Y+μ)
            # Also, when we estimate protein coding regions, they may overlap with non coding regions
            # they usually do (lncRNAs on the opposite strand), thus we have to perform the computation separately

            # # #
            # # #
            for biotype in biotypes:
                if biotype != ".":
                    gff_df_temp = gff_df[gff_df["biotype"] == biotype].copy()
                else:
                    gff_df_temp = gff_df.copy()
                if gff_df_temp.shape[0] == 0:
                    logging.warning(f"No features found for biotype `{biotype}` in GFF file `{gff_file}`. Skipping accession `{accession_id}`.")
                    continue
                gff_df_temp.loc[:, "Chromosome"] = gff_df_temp["seqID"] + ";" + gff_df_temp["compartment"]
                gff_gr = (
                        pr.PyRanges(gff_df_temp)
                        .merge(strand=False)
                        .as_df()
                        .assign(
                            seqID=lambda ds: ds["Chromosome"].str.split(";", expand=True)[0],
                            compartment=lambda ds: ds["Chromosome"].str.split(";", expand=True)[1]
                        )[["seqID", "Start", "End", "compartment"]]
                )
                gff_bed = BedTool.from_dataframe(gff_gr).sort() 
                coverage_df = ( 
                        pl.read_csv( 
                            gff_bed.coverage(df_bed).fn,
                            has_header=False,
                            separator="\t",
                            new_columns=FIELDS
                        )
                        .with_columns(
                            at_least_one=(pl.col("total_bases") > 0).cast(pl.UInt8)
                        )
                        .group_by(["compartment"], maintain_order=True)
                        .agg(
                            pl.col("total_bases").sum().alias("total_bases"),
                            pl.col("compartment_length").sum().alias("compartment_length"),
                            pl.col("total_bases").count().alias("total_compartments"),
                            (1e2 * pl.col("at_least_one").mean()).alias("pct_at_least_one"),
                            pl.col("coverage").min().alias("min_coverage"),
                            pl.col("coverage").max().alias("max_coverage"),
                            pl.col("coverage").mean().alias("avg_coverage"),
                            pl.col("coverage").median().alias("median_coverage"),
                            pl.col("coverage").std().alias("std_coverage")
                        )
                        .with_columns(
                            (1e3 * pl.col("total_bases") / pl.col("compartment_length")).alias("total_coverage")
                        )
                        .with_columns(
                            pl.lit(accession_id).alias("#assembly_accession"),
                            biotype=pl.lit(biotype),
                            pattern=pl.lit(pattern),
                        )
                )
                (
                    coverage_df.select(EXPECTED_FIELDS)
                    .to_pandas()
                    .to_csv(fin, 
                            sep="\t", 
                            index=False, 
                            header=not wrote_header)
                )
                wrote_header = True
                # coverage_df.write_csv(fin, separator="\t", include_header=file_idx==1)
            self.files_processed += 1

        # Proper thread join
        logger.join(timeout=1.0)
        fin.close()
        if isinstance(assembly_summary, str) and Path(assembly_summary).is_file():    
            logging.info(f"Found assembly summary file `{assembly_summary}`. Merging density data with assembly information.")
            assembly_summary = Path(assembly_summary).resolve()
            # headers = list(pd.read_table("headers.txt").columns)
            try:
                assembly_df = pd.read_table(assembly_summary, 
                                            dtype={"species_taxid": int,
                                                   "gc_percent": float,
                                                   "genome_size": int},
                                                   low_memory=False
                                                    )
            except Exception as e:
                logging.error(f"Error reading assembly summary file `{assembly_summary}`: {e}. Skipping merge with assembly data.")
                print(colored(f"Error reading assembly summary file `{assembly_summary}`: {e}. Skipping merge with assembly data.", "red"))
                return
            coverage_df = pd.read_table(outfile)
            coverage_df = coverage_df\
                        .merge(
                                    assembly_df,
                                    left_on="#assembly_accession",
                                    right_on="#assembly_accession",
                                    how="left"
                        )
            merged_outfile = self.outdir.joinpath(f"tss_tes_density_{pattern}_bucket_{bucket_id}_with_assembly_data.tsv.gz")
            coverage_df.to_csv(merged_outfile, mode="w", sep="\t", index=False, compression="gzip")
            logging.info(f"Merged density data with assembly summary and saved to `{merged_outfile}`.")
        print(colored(f"Motif coverage processing for pattern {pattern} in bucket {bucket_id}.", "green"))
        logging.info(f"Process has been completed succesfully (bucket {bucket_id}).")
        return
    
def main():
    import argparse 
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("schedule", type=str)
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--pattern", "-p", type=str, default="IR", choices=["IR", "MR", "DR", "STR"])
    parser.add_argument("--indir", "-i", type=str)
    parser.add_argument("--gff_indir", "-g", type=str)
    parser.add_argument("--use_biotype", action="store_true", default=False)
    parser.add_argument("--pseudogenes_to_genes", action="store_true", default=True)
    parser.add_argument("--assembly_summary", type=str, default="data/assembly_summary_with_tree.csv.gz")
    args = parser.parse_args()
    GFFMotifCoverageProcessor(
        indir=args.indir,
        schedule=args.schedule,
    ).process_bucket(
                                bucket_id=args.bucket_id,
                                pattern=args.pattern,
                                gff_indir=args.gff_indir,
                                use_biotype=args.use_biotype,
                                pseudogenes_to_genes=args.pseudogenes_to_genes
                            )