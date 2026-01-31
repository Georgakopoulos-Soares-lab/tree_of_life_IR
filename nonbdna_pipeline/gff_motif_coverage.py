import time
import polars as pl 
from pathlib import Path 
from termcolor import colored
import pandas as pd
import pyranges as pr
import threading 
import gzip
import logging
from typing import Optional, ClassVar 
import attr
from attr import field
from pybedtools import BedTool
import pybedtools
from nonbdna_pipeline.stream_and_merge_bucket import StreamAndMerge 
from nonbdna_pipeline.pwm_density import PWMExtractor 
from nonbdna_pipeline.tss_tes_processing import * 

def merge_with_summary(assembly_summary: str, outfile: str):
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
        return coverage_df 

@attr.s 
class GFFMotifCoverageProcessor(StreamAndMerge):
    files_processed: int = field(init=False, default=0)
    log_dir: Path = field(init=False)
    total_files: int = field(init=False, default=0)
    outdir: Path = field(init=False)
    gff_suffix: str = field(init=True, default=".agat.gff")
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
                            polarity: bool = False,
                            min_partition: Optional[int] = None,
                            max_partition: Optional[int] = None,
                            ) -> None:
        gff_indir = Path(gff_indir).resolve()
        partition_list = ["."] 
        motif_columns = ["Chromosome", "Start", "End"]
        if polarity:
            motif_columns.append("Strand")
            raise NotImplementedError("Strand-specific processing not yet implemented.")
        if isinstance(partition_col, str) and isinstance(min_partition, int) and isinstance(max_partition, int):
            partition_list += list(range(min_partition, max_partition))
            motif_columns.append("partition")
        else:
            partition_col = None
            min_partition = None 
            max_partition = None
        # LOGGING
        tracker = GFFMotifCoverageProcessor._TrackProgress(
            bucket_id=bucket_id,
            total_records=0,
            log_indir=self.log_dir,
        )
        logger = threading.Thread(target=tracker._log_progress, 
                                  daemon=True,
                                  name="LoggerDaemon")
        # LOAD FILES
        FIELDS = ["seqID", "start", "end", "compartment", "merged_count"] + GFFMotifCoverageProcessor.COVERAGE_FIELDS
        EXPECTED_FIELDS = ["#assembly_accession", 
                    "biotype", 
                    "partition",
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
                    "pct_at_least_one",
                    "total_merged_members"]
        infiles = self.load_validated_files_from_log(bucket_id=bucket_id, pattern=pattern)
        tracker.total_records = len(infiles)
        outfile = self.outdir.joinpath(f"gff_motif_coverage_{pattern}_bucket_{bucket_id}.tsv.gz")
        logging.info(f"Initializing extraction process for pattern {pattern} (bucket {bucket_id}).")
        fin = gzip.open(outfile, mode="wt")
        wrote_header = False 
        logger.start()
        for file_idx, infile in enumerate(infiles, start=1):
            tracker.track += 1
            accession_name = StreamAndMerge.extract_name(infile, pattern)
            accession_id = StreamAndMerge.extract_id(infile)
            extraction_file = self.indir.joinpath(accession_name + f"_{pattern}.processed.tsv")
            gff_file = gff_indir.joinpath(accession_name + self.gff_suffix)
            if not gff_file.is_file():
                gff_file = gff_file.with_suffix(f"{self.gff_suffix}.gz")
                if not gff_file.is_file():
                    continue 
            gff_df = read_gff(gff_file, 
                              end_to_one_bp=False, 
                              pseudogenes_to_genes=pseudogenes_to_genes,
                              parse_biotype=use_biotype,
                              filter_on=self.valid_compartments
            ).rename(columns={
                              "start": "Start",
                              "end": "End"
                              })
            if gff_df.shape[0] == 0:
                logging.warning(f"No features found in GFF file `{gff_file}`. Skipping accession `{accession_id}`.")
                continue 
            biotypes = self.biotypes if use_biotype else ["."]
            df = self.read_motifs(extraction_file)
            if df.shape[0] == 0:
                gff_df = gff_df.rename(columns={"seqID": "Chromosome"})
                logging.warning(f"No motifs found in file `{extraction_file}`.")
                for partition in partition_list:
                    for biotype in biotypes:
                        if biotype != ".":
                            gff_biotype_df = gff_df[gff_df["biotype"] == biotype].copy()
                        else:
                            gff_biotype_df = gff_df.copy()
                        if gff_biotype_df.shape[0] == 0:
                            continue
                        compartments = gff_biotype_df["compartment"].unique().tolist()
                        for compartment in compartments:
                            gff_df_temp = gff_biotype_df[gff_biotype_df["compartment"] == compartment].copy()
                            total_merged_members = gff_df_temp.shape[0]
                            if total_merged_members == 0:
                                continue
                            merged_ranges = pr.PyRanges(gff_df_temp[["Chromosome", "Start", "End"]]).merge(strand=False)
                            merged_df = merged_ranges.as_df()
                            lengths = (merged_df["End"] - merged_df["Start"]).clip(lower=0)
                            compartment_length = int(lengths.sum())
                            total_compartments = merged_df.shape[0]
                            coverage_df = pd.DataFrame({
                                        "#assembly_accession": [accession_id],
                                        "pattern": [pattern],
                                        "partition": [partition],
                                        "compartment": [compartment],
                                        "biotype": [biotype],
                                        "total_bases": [0],
                                        "compartment_length": [compartment_length],
                                        "total_compartments": [total_compartments],
                                        "pct_at_least_one": [0.0],
                                        "min_coverage": [0.0],
                                        "max_coverage": [0.0],
                                        "avg_coverage": [0.0],
                                        "median_coverage": [0.0],
                                        "std_coverage": [0.0],
                                        "total_merged_members": [total_merged_members],
                                        "total_coverage": [0.0]
                                    })
                            coverage_df[EXPECTED_FIELDS].to_csv(
                                fin,
                                sep="\t",
                                index=False,
                                header=not wrote_header,
                            )
                            wrote_header = True
                continue
            if partition_col not in df.columns:
                raise KeyError(f"Partition column `{partition_col}` not found in motif extraction file `{extraction_file}` for accession `{accession_id}`.")
            if partition_col:
                df = df.rename(columns={partition_col: "partition"})
            df = df.rename(columns={"seqID": "Chromosome", 
                                    "start": "Start", 
                                    "end": "End", 
                                    "strand": "Strand"})
            df_gr = pr.PyRanges(df[motif_columns])
            # Why we need to loop
            # We are attempting to calculate the proportion of:
            # - Total Genic Region covered by motif base pairs 
            # Thus, we want this number to be merged, otherwise we will be counting overlapping regions multiple times
            # Instead of (X/Y) we will computing (X+ε/Y+μ)
            # Also, when we estimate protein coding regions, they may overlap with non coding regions
            # they usually do (lncRNAs on the opposite strand), thus we have to perform the computation separately
            for biotype in biotypes:
                if biotype != ".":
                    gff_df_temp = gff_df[gff_df["biotype"] == biotype].copy()
                else:
                    gff_df_temp = gff_df.copy()
                if gff_df_temp.shape[0] == 0:
                    logging.warning(f"No features found for biotype `{biotype}` in GFF file `{gff_file}`. Skipping accession `{accession_id}`.")
                    continue
                compartments = gff_df_temp["compartment"].unique().tolist()
                if biotype != ".":
                    if len(compartments) > 1:
                        raise ValueError(f"Multiple compartments found for biotype `{biotype}` in GFF file `{gff_file}`. Expected a single compartment per biotype.")
                    compartment = compartments[0]
                    if compartment != "Gene":
                        raise ValueError(f"Unexpected compartment `{compartment}` for biotype `{biotype}` in GFF file `{gff_file}`. Expected `gene` compartment.")
                gff_df_temp.loc[:, "Chromosome"] = gff_df_temp["seqID"] + ";" + gff_df_temp["compartment"]
                # Merge Starts
                orig_pr = pr.PyRanges(gff_df_temp[["Chromosome", "Start", "End"]])
                merged_pr = orig_pr.merge(strand=False)
                merged_df = merged_pr.as_df()
                # Merge Ends
                compartment_counts = gff_df_temp["compartment"].value_counts().to_dict()
                merged_df["compartment"] = merged_df["Chromosome"].str.split(";", expand=True)[1]
                merged_df["seqID"] = merged_df["Chromosome"].str.split(";", expand=True)[0]
                merged_df.loc[:, "merged_count"] = merged_df["compartment"].map(compartment_counts)
                gff_gr = merged_df[["seqID", "Start", "End", "compartment", "merged_count"]]
                gff_bed = BedTool.from_dataframe(gff_gr).sort() 
                for partition in partition_list:
                    if partition != ".":
                        df_partition = df_gr[df_gr.partition == partition]
                    else:
                        df_partition = df_gr
                    df_bed = BedTool.from_dataframe(
                                                    df_partition
                                                    .merge(strand=False)
                                                    .as_df()
                                                    ).sort()
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
                                (1e2 * pl.col("at_least_one").mean()).round(3).alias("pct_at_least_one"),
                                pl.col("coverage").min().round(3).alias("min_coverage"),
                                pl.col("coverage").max().round(3).alias("max_coverage"),
                                pl.col("coverage").mean().round(3).alias("avg_coverage"),
                                pl.col("coverage").median().round(3).alias("median_coverage"),
                                pl.col("coverage").std().round(3).alias("std_coverage"),
                                pl.col("merged_count").sum().alias("total_merged_members")
                            )
                            .with_columns(
                                (1e3 * pl.col("total_bases") / pl.col("compartment_length")).round(3).alias("total_coverage")
                            )
                        .with_columns(
                            pl.lit(accession_id).alias("#assembly_accession"),
                            biotype=pl.lit(biotype),
                            pattern=pl.lit(pattern),
                            partition=pl.lit(partition),
                        )
                        .select(EXPECTED_FIELDS)
                        .to_pandas()
                    )
                    (
                        coverage_df.to_csv(
                                fin, 
                                sep="\t", 
                                index=False, 
                                header=not wrote_header)
                    )
                    wrote_header = True
            self.files_processed += 1
        logger.join(timeout=1.0)
        fin.close()
        coverage_df = merge_with_summary(assembly_summary=assembly_summary,
                                          outfile=outfile)
        if coverage_df is not None:
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
    parser.add_argument("--gff_suffix", type=str, default=".agat.gff")
    parser.add_argument("--use_biotype", action="store_true", default=False)
    parser.add_argument("--partition_col", type=str, default=None)
    parser.add_argument("--min_partition", type=int, default=None)
    parser.add_argument("--max_partition", type=int, default=None)
    parser.add_argument("--tmpdir", type=str, default="garbage")
    parser.add_argument("--pseudogenes_to_genes", action="store_true", default=True)
    parser.add_argument("--assembly_summary", type=str, default="data/assembly_summary_with_tree.csv.gz")
    args = parser.parse_args()

    tmpdir = Path(args.tmpdir)
    tmpdir.mkdir(exist_ok=True)
    pybedtools.helpers.set_tempdir(tmpdir)

    GFFMotifCoverageProcessor(
        indir=args.indir,
        schedule=args.schedule,
        gff_suffix=args.gff_suffix
    ).process_bucket(
                                bucket_id=args.bucket_id,
                                pattern=args.pattern,
                                gff_indir=args.gff_indir,
                                use_biotype=args.use_biotype,
                                pseudogenes_to_genes=args.pseudogenes_to_genes,
                                partition_col=args.partition_col,
                                min_partition=args.min_partition,
                                max_partition=args.max_partition,
        )
    pybedtools.helpers.cleanup(remove_all=True)
