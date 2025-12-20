# Enrichment Pipepline - IGS Lab
# Nikol <3

import numpy as np
import pandas as pd
import polars as pl
import json
import csv
import threading
from tqdm import tqdm
from termcolor import colored
import logging
import attr
from attr import field
from scipy.optimize import curve_fit
from functools import partial
from pybedtools import BedTool
import pybedtools
from gff_utils import CoverageExtractor, GFFExtractor, Expander, evaluate_motif_strand, PolarityError, extract_id
from pwm_density import PWMExtractor
from typing import Optional, Callable, Iterable
from pathlib import Path

def distance_from_origin(intersect_df: pd.DataFrame, 
                         polynomial: Callable,
                         window_size: int,
                         step: int = 500) -> tuple[pd.DataFrame, pd.Series, np.ndarray]:
    # intersect_df = intersect_df.query("overlap > 0")
    intersect_df["origin"] = intersect_df["start"] + window_size 
    intersect_df["distance"] = np.minimum(
                            np.abs(intersect_df["motif_start"] - intersect_df["origin"]),
                            np.abs(intersect_df["motif_end"] - 1 - intersect_df["origin"])
                    )
    intersect_df["distance_bin"] = intersect_df["distance"].apply(lambda x: x//step+1)
    counts_per_bin = intersect_df.groupby("distance_bin")\
                                 .agg(totalCounts=("seqID", "count"),)
    x_data = counts_per_bin["distance_bin"]
    y_data = counts_per_bin["totalCounts"]
    params, _ = curve_fit(polynomial, x_data, y_data)
    y_pred = polynomial(x_data, *params)
    return counts_per_bin, x_data, y_pred

def bootstrap(intersect_df: pd.DataFrame, 
              window_size: int,
              N: int = 1_000,
              lower_q: float = 0.025,
              upper_q: float = 0.975,
              ) -> tuple[pd.Series, pd.Series, pd.Series]:
    if isinstance(intersect_df, pl.DataFrame):
        intersect_df = intersect_df.to_pandas()
    bootstrapped_df = []
    extractor = PWMExtractor()
    for _ in range(N):
        sample_df = intersect_df.sample(frac=1.0, replace=True)
        density = extractor.extract_density(sample_df, 
                                            window_size=window_size,
                                            return_array=True,
                                            enrichment=True
                                            )
        bootstrapped_df.append(density)
    bootstrapped_df = pd.DataFrame(bootstrapped_df)
    average = bootstrapped_df.mean()
    ci_lower = bootstrapped_df.quantile(lower_q)
    ci_upper = bootstrapped_df.quantile(upper_q)
    return average, ci_lower, ci_upper

@attr.s
class DensityExtractor:
    
    mode: str = field(default="density", converter=str)
    format: str = field(default="txt")
    compression: str = field(default="snappy")
    use_pyarrow: bool = field(default=False, converter=bool)
    statistics: bool = field(default=True, converter=bool)
    tempdir: Optional[str] = field(default=None)
    schedule: Optional[str] = field(default=None)
    design: Optional[str] = field(default=None)
    transcription_site_loci: tuple[str, str] = field(init=False)
    float_precision: int = field(default=3, converter=int)
    biotypes: Optional[list[str]] = field(factory=lambda : ["protein_coding",
                                                            "non_coding",
                                                            "."]
                                          )
    extractions: dict[str, str] = field(init=False, repr=False)
    window_size: int = field(default=500, 
                             converter=int, 
                             validator=attr.validators.instance_of(int))
    empty_accessions: set[str] = field(init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        if self.mode != "polarity" and self.mode != "density":
            raise ValueError(f"Invalid mode `{self.mode}` detected.")
        # compression mode for parquet files
        # {‘lz4’, ‘uncompressed’, ‘snappy’, ‘gzip’, ‘lzo’, ‘brotli’, ‘zstd’}
        if self.compression != "snappy" \
                and self.compression != "gzip" \
                and self.compression != "lz4" \
                and self.compression != "uncompressed" \
                and self.compression != "zstd" \
                and self.compression != "lzo" \
                and self.compression != "brotli":
                    raise ValueError(f"Invalid parquet compression selected: `{self.compression}`.")

        if self.tempdir is not None:
            pybedtools.helpers.set_tempdir(self.tempdir)
            pybedtools.set_tempdir(self.tempdir)
            self.tempdir = Path(self.tempdir).resolve()
        if self.biotypes is None:
            self.biotypes = ["."]
        if self.schedule:
            self.schedule = Path(self.schedule).resolve()

        self.transcription_site_loci = ("start", "end")
        self.extractions = dict()
        self.empty_accessions = set()

        if not self.design:
            return

        delimiter = CoverageExtractor._sniff_delimiter(self.design)
        with open(self.design, mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                self.extractions[row["accession_id"]] = row["extraction"]
            
    def preprocess_extraction_df(self, extractions_df: pl.DataFrame, partition_col: Optional[str] = None, polarity_mode: str = "GC") -> pl.DataFrame:
        # Check if chromosomes are OK
        if "seqID" not in extractions_df.columns and "chromosome" in extractions_df.columns:
            extractions_df = extractions_df.rename({"chromosome": "seqID"})
        elif "seqID" not in extractions_df.columns and "chrom" in extractions_df.columns:
            extractions_df = extractions_df.rename({"chrom": "seqID"})
        elif "seqID" not in extractions_df.columns and "sequence_name" in extractions_df.columns:
            extractions_df = extractions_df.rename({"sequence_name": "seqID"})
        elif "seqID" not in extractions_df.columns and "Sequence_name" in extractions_df.columns:
            extractions_df = extractions_df.rename({"Sequence_name": "seqID"})
        elif "seqID" not in extractions_df:
            raise KeyError(f"Column `seqID` is not present in the extractions dataframe.")

        # end coordinate exists
        if "end" not in extractions_df.columns:
            extractions_df = extractions_df.rename({"stop": "end"})
        if "end" not in extractions_df:
            raise KeyError(f"Column `end` is not present in the extractions dataframe ({extraction}).")

        # Check if polarity can be determined correctly
        if "motif_strand" not in extractions_df.columns \
                and "strand" not in extractions_df.columns \
                and "sequence" not in extractions_df.columns \
                and self.mode == "polarity":
            raise PolarityError(f"Failure to detect polarity. Either `strand` or `sequence` must be present in the dataframe when selected mode is `template`.")
        elif "strand" not in extractions_df and "motif_strand" not in extractions_df.columns and self.mode == "polarity":
            # resolve strand if sequence is present when mode is `template`
            extractions_df = extractions_df.with_columns(
                                pl.col("sequence").str.to_uppercase()
                                .map_elements(
                                    partial(evaluate_motif_strand, polarity_mode), return_dtype=str
                                    )
                                .alias("motif_strand")
                            )
        elif "strand" in extractions_df.columns:
            extractions_df = extractions_df.rename({"strand": "motif_strand"})
        
        # select relevant columns
        if self.mode == "polarity":
            COLS = ["seqID", "start", "end", "motif_strand"]
        else:
            COLS = ["seqID", "start", "end"]
        if partition_col:
            COLS.append(partition_col)
        return extractions_df.select(COLS)

    def parse_density(self, A_df: pl.DataFrame, 
                              B_df: pl.DataFrame, 
                              site: Optional[str] = None,
                              merge_A: bool = False,
                              merge_B: bool = True,
                              enrichment: bool = False) -> pl.DataFrame:
        A_bed = BedTool.from_dataframe(A_df.to_pandas())
        if self.mode == "polarity" and "motif_strand" not in B_df:
            B_df = B_df.rename({"strand": "motif_strand"})
        B_bed = BedTool.from_dataframe(B_df.to_pandas())
        if B_df.shape[0] > 0:
            B_columns = B_df.columns[3:]
        else:
            B_columns = []
        # merge if applicable
        if merge_A:
            A_bed = A_bed.sort().merge()
        if merge_B and not self.mode == "polarity":
            B_bed = B_bed.sort().merge()
            B_columns = []
        elif merge_B and self.mode == "polarity":
            if B_df.columns[5] != "strand":
                raise PolarityError("Invalid strand location.")
            B_bed = B_bed.sort().merge(s=True, c="6", o="distinct")
            B_columns = ["motif_strand"]
        # sort bed
        A_bed = A_bed.sort()
        B_bed = B_bed.sort()
        # fetch columns
        A_columns = A_df.columns
        # Calculate intersection between target: A_df and query: B_df
        new_columns = A_columns + ["chrom", "motif_start", "motif_end"] + B_columns + ["overlap"]
        try:
            intersect_df = pl.read_csv(
                                    A_bed.intersect(B_bed, wao=True).fn,
                                    has_header=False,
                                    separator="\t",
                                    new_columns=new_columns,
                                    null_values="."
                                )
        except pl.errors.ComputeError as e:
            print(new_columns)
            raise pl.errors.ComputeError(e)
        # total queries & targets
        total_queries = A_bed.count()
        total_targets = B_bed.count()
        # Estimate:
        # - total number of compartments from A (unique)
        # - the percentage that have at least one overlap
        matched_queries = ( 
                           intersect_df.group_by(A_columns)
                                       .agg(
                                            # pl.col("overlap").count().alias("total_matches"),
                                            pl.col("overlap").sum().alias("match_bp")
                                        )
                                        .with_columns(
                                                (pl.col("match_bp") > 0).cast(pl.Int8)
                                                                        .alias("atLeastOneMatch")
                                        )
                                        .get_column("atLeastOneMatch")
                                        .value_counts()
                                        .filter(pl.col("atLeastOneMatch") == 1)
                                        .get_column("count")
                    )
        if len(matched_queries) == 0:
            total_matched_queries = 0
        else:
            total_matched_queries = matched_queries[0]
        
        # keep only overlapping hits
        intersect_df = intersect_df.filter(pl.col("overlap") > 0)
        total_motifs = B_df.shape[0]
        # keep compartments from A succesfully mapped to an element of B
        # drop duplicates to ensure that a query from B has not been mapped to more 
        # than one compartment from A; since here we are only interested to estimate 
        # the percentage of motifs from B that map to AT LEAST one compartment.
        # In the previous question, we examined the percentage of compartments from A 
        # that map to AT LEAST ONE from B.
        total_matched_targets = intersect_df.unique(subset=['chrom', 'motif_start', 'motif_end']).shape[0]
        matched_targets_template: Optional[float]
        matched_targets_non_template: Optional[float]
        if self.mode == "polarity":
            if intersect_df.shape[0] > 0:
                intersect_df = (
                            intersect_df.filter(
                                                (pl.col("strand") == "+") | (pl.col("strand") == "-"),
                                                (pl.col("motif_strand") == "+") | (pl.col("motif_strand") == "-")
                                                )
                                        .with_columns(
                                                    pl.when(pl.col("strand") == pl.col("motif_strand"))
                                                        .then(pl.lit("Non-Template"))
                                                        .otherwise(pl.lit("Template"))
                                                        .alias("strand_polarity")
                                            )
                            )
                # non template proportion
                matched_targets_homopolarity = round(1e2 * intersect_df.with_columns(
                                                                        pl.when(pl.col("strand_polarity") == "Non-Template")
                                                                          .then(1)
                                                                          .otherwise(0)
                                                                          .alias("is_non_template")
                                                                    )["is_non_template"].mean(), 
                                                     2)
                # template proportion
                matched_targets_heteropolarity = round(1e2 - matched_targets_homopolarity, 2)
            else:
                matched_targets_homopolarity = 0.0
                matched_targets_heteropolarity = 0.0
        else:
                matched_targets_homopolarity = float("nan")
                matched_targets_heteropolarity = float("nan")

        target_perc: Optional[float] = round(1e2 * total_matched_targets / total_targets, 2) if total_targets > 0 else float("nan")
        query_perc: Optional[float] = round(1e2 * total_matched_queries / total_queries, 2) if total_queries > 0 else float("nan")

        # calculate density relative to an origin around a prespecified symmetrical window
        pwm = PWMExtractor()
        if self.mode == "polarity":
            density = pwm.extract_template_density(intersect_df, 
                                                   window_size=self.window_size, 
                                                   return_frame=False,
                                                   enrichment=enrichment,
                                                )
        else:
            density = pwm.extract_density(intersect_df, 
                                          window_size=self.window_size, 
                                          return_frame=False,
                                          enrichment=enrichment)
        stats = {
                "total_targets": total_motifs,
                "matched_targets": total_matched_targets,
                "matched_targets_heteropolarity": matched_targets_heteropolarity,
                "matched_targets_homopolarity": matched_targets_homopolarity,
                "target_perc": target_perc,
                "total_queries": total_queries,
                "matched_queries": total_matched_queries,
                "query_perc": query_perc,
                }
        if site:
            stats.update({"site": site})
        if self.mode == "polarity":
            density_df = pl.DataFrame({
                                "strand_polarity": ["Template", "Non-Template"],
                                **{str(int(i) - self.window_size): [density["Template"][i], density["Non-Template"][i]] for i in range(2*self.window_size+1)},
                                **{attr: [val, val] for attr, val in stats.items()},
                                })
        else:
            density_df = pl.DataFrame({
                                    **{str(int(i) - self.window_size): [density[i]] for i in range(2*self.window_size+1)},
                                    **{attr: [val] for attr, val in stats.items()},
                                    })
        return density_df

    def process_site(self, extractions_table: pl.DataFrame,
                            gff_table: pl.DataFrame,
                            loci: str,
                            enrichment: bool = False,
                            merge_motifs: bool = True,
                            return_df: bool = True,
                            partition_col: Optional[str] = None,
                            partition_groups: Optional[list[str]] = None,
                            accession_id: Optional[str] = None,
                            ) -> pl.DataFrame:
        if loci != "start" and loci != "end":
            raise ValueError(f"Invalid transcription site loci `{loci}`. Select either `start` or `end`.")

        densities_df = []
        window_expander = Expander(window_size=self.window_size)
        gff_expanded = window_expander.expand_windows(gff_table, loci=loci)
        site = "TSS" if loci == "start" else "TES"
        densities_df: list[pl.DataFrame] = []
        for biotype in self.biotypes:
            if biotype != ".":
                gff_expanded_temp = gff_expanded.filter(pl.col("biotype") == biotype)
            else:
                gff_expanded_temp = gff_expanded

            if gff_expanded_temp.shape[0] == 0:
                logging.warning(f"GFF file `{accession_id}` for biotype `{biotype}` is empty.")
                continue

            if partition_col:
                for partition in set(partition_groups):
                    if extractions_table.shape[0] > 0:
                        extractions_df_temp = extractions_table.filter(pl.col(partition_col) == partition)
                    else:
                        extractions_df_temp = extractions_table
                    # Transcription Start Site (TSS) from GFF
                    site_density = self.parse_density(gff_expanded_temp, 
                                                     extractions_df_temp, 
                                                     site=site,
                                                     merge_B=merge_motifs,
                                                     enrichment=enrichment)\
                                            .with_columns(
                                                        pl.lit(biotype).alias("biotype"),
                                                        pl.lit(site).alias("site"),
                                                )
                    if partition_col:
                        site_density = site_density.with_columns(
                                                       pl.lit(str(partition)).alias(partition_col)
                                                       )
                    if accession_id:
                        site_density = site_density.with_columns(
                                                            pl.lit(accession_id).alias("#assembly_accession")
                                                        )
                    densities_df.append(site_density)

            # Transcription Start or End Site (TSS/TES) from GFF
            site_density = self.parse_density(gff_expanded_temp, 
                                              extractions_table,
                                              site=site,
                                              merge_B=merge_motifs,
                                              enrichment=enrichment)\
                                        .with_columns(
                                                pl.lit(biotype).alias("biotype"),
                                                pl.lit(site).alias("site")
                                        )
            if partition_col:
                site_density = site_density.with_columns( 
                                                pl.lit("generic").alias(partition_col)
                                        )
            if accession_id:
                site_density = site_density.with_columns(
                                                pl.lit(accession_id).alias("#assembly_accession")
                                            )
            densities_df.append(site_density)
        if return_df:
            densities_df = pl.concat(densities_df)
        return densities_df

    def process(self, extraction: str,
                        gff_file: str,
                        enrichment: bool = False,
                        compartment: str = "Gene",
                        partition_col: Optional[str] = None,
                        partition_groups: Optional[list] = None,
                        accession_id: Optional[str] = None,
                        polarity_mode: str = "GC",
                        biotype: bool = True,
                        return_df: bool = False
                    ) -> pl.DataFrame | list[pl.DataFrame]:
        densities_df = []
        compartment = compartment.lower()
        gff_reader = GFFExtractor(compartments=[compartment])
        delimiter = CoverageExtractor._sniff_delimiter(extraction)
        accession_id = extract_id(gff_file)
        try:
            extractions_df = pl.read_csv(extraction, 
                                         separator=delimiter, 
                                         comment_prefix="#")
            if extractions_df.shape[0] > 0:
                extractions_df = extractions_df.filter(pl.col("sequence_length") >= 10)
                # extractions_df = extractions_df.with_columns(STR_coverage_sequence=pl.col("STR_coverage_sequence").cast(int))
            if extractions_df.shape[0] == 0:
                is_empty = True
            else:
                is_empty = False
            # change column names to match
            # extractions_df.columns = [col[0].lower() + col[1:] for col in extractions_df.columns]
        except pl.exceptions.NoDataError as e:
            # if there are no extractions, then density is equal to 0
            # return an empty list
            # since this doesn't affect the results ?
            logging.info(f"Accession ID `{accession_id}` was found empty.\n{e}")
            is_empty = True
        except ValueError as e:
            # probably delimiter issue
            logging.info(f"Invalid delimiter. Accession ID `{accession_id}` was found empty.\n{e}")
            is_empty = True

        if is_empty:
            logging.warning(f"Failed to process gff file `{gff_file}`. Reason: empty dataset.")
            self.empty_accessions.add(accession_id)
            schema = {
                        "seqID": pl.Utf8,
                        "start": pl.Int32,
                        "end": pl.Int32
                    }
            if partition_col is not None:
                schema.update({partition_col: pl.Utf8})
            if self.mode == "polarity":
                schema.update({"motif_strand": pl.Utf8})
            extractions_df = pl.DataFrame([], schema=schema)
            self.empty_accessions.add(accession_id)
            # return densities_df
        extractions_df = self.preprocess_extraction_df(extractions_df, polarity_mode=polarity_mode, partition_col=partition_col)
        # Check if there are partitions
        if partition_col is not None and not isinstance(partition_groups, list):
            raise TypeError(f"Partition groups must be a list when a partition column `{partition_col}` has been chosen.")
        if partition_col is not None and len(partition_groups) == 0:
            raise ValueError("Partition groups cannot be an empty list when a partition column `{partition_col}` has been chosen.")
        if partition_col is not None and partition_col not in extractions_df:
            raise KeyError(f"Column `{partition_col}` is not present in the extractions dataframe ({extraction}).")
        
        # Continue with density extraction
        gff_df = gff_reader.read_gff(gff_file,
                                     change_names=True,
                                     end_one_base=True,
                                     parse_biotype=biotype, 
                                     # Please do not remove this because viral annotations will be processed 
                                     # erroneously if they have multiple regions
                                     # it will generate duplicate matches!
                                     join_region=False,
                                    )
        # handle that some GFF files may not have genes or other compartments
        if gff_df is None or gff_df.shape[0] == 0:
            logging.warning(f"Failed to process gff file `{gff_file}`. Reason: no records of `{compartment}` were found in the annotation file.")
            # in this case return an empty list
            # since we cannot make a statistical assesement about the density
            # there are no compartments to make such assesement
            return densities_df
        if compartment == "exon":
            # Need to filter out TSS/TES from Exon collection
            gff_df = gff_reader.drop_first_last(gff_df, compartment=compartment)

        gff_df = gff_df.select(["seqID", "start", "end", "biotype", "phase", "strand"])
        for loci in self.transcription_site_loci:
            densities_table = self.process_site(extractions_table=extractions_df,
                                                gff_table=gff_df,
                                                loci=loci,
                                                enrichment=enrichment,
                                                return_df=False,
                                                partition_col=partition_col,
                                                partition_groups=partition_groups,
                                                accession_id=accession_id
                                                )
            densities_df.extend(densities_table)
        # parse densities 
        if return_df:
            densities_df = pl.concat(densities_df)
        return densities_df

    def load_bucket(self, bucket_id: int) -> list[str]:
        with open(self.schedule, mode="r", encoding="utf-8") as f:
            return json.load(f)[str(bucket_id)]

    def process_bucket(self, bucket_id: int,
                            out: str = "",
                            sleeping_time: float = 200.0,
                            compartment: str = "gene",
                            enrichment: bool = False,
                            polarity_mode: str = "GC",
                            biotype: bool = True,
                            partition_col: Optional[str] = None,
                            partition_groups: Optional[list] = None,
                       ) -> None:
        bucket = self.load_bucket(bucket_id=bucket_id)
        logging.info(f"Processing bucket `{bucket_id}`...")
        print(f"Processing bucket `{bucket_id}`...")
        tracker = CoverageExtractor._TrackProgress(bucket_id=bucket_id,
                                                   total_records=len(bucket),
                                                   sleeping_time=sleeping_time
                                                   )
        daemon = threading.Thread(target=tracker.start, daemon=True, name="LoggingDensityDaemon")
        daemon.start()
    
        # LOAD
        MAXIMUM_STREAMING_LOAD = 10
        densities = []

        # COLLECTION
        collection_dest = dict()
        collection_handlers = dict()
        first_headers = dict()
        sites = ["TSS", "TES"]
        if self.mode == "density":
            for site in sites:
                for _biotype in self.biotypes:
                    collection_dest[site, _biotype] = f"{out}/enrichment_bucket_{bucket_id}_{window_size}_{compartment}_{self.mode}_{polarity_mode}.{_biotype}.{site}.{self.format}"
            dest_empty = f"{out}/empty_accessions_bucket_{bucket_id}_{window_size}_{compartment}_{self.mode}.txt"
        else:
            for site in sites:
                for _biotype in self.biotypes:
                    if _biotype == ".":
                        _biotype = "generic"
                    collection_dest[site, _biotype] = f"{out}/enrichment_bucket_{bucket_id}_{window_size}_{compartment}_{self.mode}_{polarity_mode}.{_biotype.replace('.', 'generic')}.{site}.{self.format}"
            dest_empty = f"{out}/empty_accessions_bucket_{bucket_id}_{window_size}_{compartment}_{self.mode}_{polarity_mode}.txt"
        # # Initialize Collection # # 
        for site in sites:
            for _biotype in self.biotypes:
                collection_handlers[site, _biotype] = open(collection_dest[site, _biotype], mode="w", encoding="UTF-8")
                first_headers[site, _biotype] = True
                at_ = collection_dest[site, _biotype]
                logging.info(f"Saving density output for {site} & biotype={_biotype} at `{at_}`... (bucket {bucket_id}).")
                print(colored(f"Saving density output for {site} & biotype={_biotype} at `{at_}`... (bucket {bucket_id}).", "green"))

        first_write = True
        for gff in tqdm(bucket):
            logging.info(f"Processing accession `{gff}` (bucket {bucket_id}).")
            accession_id = extract_id(gff)
            extraction_filename = self.extractions.get(accession_id)
            if extraction_filename is None:
                logging.info(f"Failed to find extraction file for accession id `{accession_id}`.")
                tqdm.write(f"Failed to find extraction file for accession id `{accession_id}`.")
                continue
            tracker.track += 1
            densities_table = self.process(
                                    extraction=extraction_filename,
                                    gff_file=gff,
                                    enrichment=enrichment,
                                    compartment=compartment,
                                    polarity_mode=polarity_mode,
                                    partition_col=partition_col,
                                    partition_groups=partition_groups,
                                    biotype=biotype,
                                    accession_id=accession_id,
                                    return_df=False
                                    # genome=genome,
                            )
            densities.extend(densities_table)
            EXCEEDED_STREAMING_LOAD = len(densities) > MAXIMUM_STREAMING_LOAD
            if EXCEEDED_STREAMING_LOAD:
                densities_df = pl.concat(densities)
                for site in sites:
                    for _biotype in self.biotypes:
                        (
                            densities_df
                                .filter(pl.col("site") == site, 
                                        pl.col("biotype") == _biotype)
                                .write_csv(
                                    collection_handlers[site, _biotype],
                                    separator="\t",
                                    include_header=first_headers[site, _biotype],
                                    float_precision=self.float_precision
                                )
                        )
                        first_headers[site, _biotype] = False
                densities = []
        if densities:
            densities_df = pl.concat(densities)
            for site in sites:
                for _biotype in self.biotypes:
                    (
                        densities_df
                        .filter(
                            pl.col("site") == site,
                            pl.col("biotype") == _biotype
                        ).write_csv(
                            collection_handlers[site, _biotype],
                            separator="\t",
                            include_header=first_headers[site, _biotype],
                            float_precision=self.float_precision
                        )
                    )
            densities = []
        # elif self.format == "parquet":
        #    densities_df.filter(pl.col("site") == "TSS")\
        #                .write_parquet(dest_TSS,
        # compression=self.compression,
        #                                use_pyarrow=self.use_pyarrow,
        #                                statistics=self.statistics
        #                                )
        #    densities_df.filter(pl.col("site") == "TES")\
        #           .write_parquet(dest_TES,
        #                                compression=self.compression,
        #                                use_pyarrow=self.use_pyarrow,
        #                                statistics=self.statistics
        #                                )
        logging.info(f"Closing files... (bucket {bucket_id}).")
        for site in sites:
            for _biotype in self.biotypes:
                collection_handlers[site, _biotype].close()
            logging.info(f"Closing file for {site} & biotype={_biotype} at `{at_}`... (bucket {bucket_id}).")
            print(colored(f"Closing file for {site} & biotype={_biotype} at `{at_}`... (bucket {bucket_id}).", "green"))
        with open(dest_empty, mode="w", encoding="utf-8") as f:
            for accession_id in self.empty_accessions:
                f.write(accession_id + "\n")
        logging.info(f"Output has been saved succesfully (bucket {bucket_id}).")
        logging.info(f"Bucket `{bucket_id}` has been processed succesfully.")

        print(colored(f"Output has been saved succesfully (bucket {bucket_id}).", "green"))
        print(colored(f"Bucket `{bucket_id}` has been processed succesfully.", "green"))
        # empty the empty accessions
        self.empty_accessions = set()
        return 


if __name__ == "__main__":
    import argparse
    parser  = argparse.ArgumentParser(description="""Extracts motif density across genomic subcompartment of interest around a specified window.""")
    parser.add_argument("schedule", type=str, default="schedule.json")
    parser.add_argument("--gff", type=str)
    parser.add_argument("--out", type=str, default="enrichment_out")
    parser.add_argument("--extraction", type=str)
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--design", type=str, default=None)
    parser.add_argument("--format", type=str, default="txt", choices=["parquet", "txt"])
    parser.add_argument("--compression", type=str, default="snappy")
    parser.add_argument("--statistics", type=int, default=1, choices=[0, 1])
    parser.add_argument("--use_pyarrow", type=int, default=0, choices=[0, 1])
    parser.add_argument("--window_size", "-w", type=int, default=500)
    parser.add_argument("--mode", "-m", type=str, choices=["polarity", "density"], default="density")
    parser.add_argument("--partition_col", default=None)
    parser.add_argument("--partition_groups", default=list(range(0, 9)), type=list, nargs="+")
    parser.add_argument("--float_precision", type=int, default=2)
    parser.add_argument("--polarity_mode", default="GC", choices=["GC", "GA", "GT"], type=str)
    parser.add_argument("--compartment", "-c", type=str, choices=["gene", "exon"], default="gene")
    parser.add_argument("--biotype", choices=[0, 1], default=1, type=int)
    parser.add_argument("--sleeping_time", default=200.0, type=float)
    parser.add_argument("--tempdir", type=str, default=None)

    args = parser.parse_args()
    schedule = args.schedule
    bucket_id = args.bucket_id
    gff = args.gff
    extraction = args.extraction
    mode = args.mode
    design = args.design
    window_size = args.window_size
    biotype = bool(args.biotype)
    compartment = args.compartment
    out = Path(args.out).resolve()
    from pybedtools.helpers import set_tempdir
    set_tempdir("/scratch/10904/nikolchanchan/tmp")
    if args.partition_col:
        out = out.joinpath(f"mode_{mode}_partition_{args.partition_col}")
        biolog_file = out.joinpath(f"biologs/enrichment_bucket_{bucket_id}_mode_{mode}_partition_{args.partition_col}.log")
    else:
        out = out.joinpath(f"mode_{mode}")
        biolog_file = out.joinpath(f"biologs/enrichment_bucket_{bucket_id}_mode_{mode}.log")
    out.mkdir(exist_ok=True, parents=True)
    biolog_file.parent.mkdir(exist_ok=True)

    logging.basicConfig(
                        level=logging.INFO,
                        filemode="w",
                        format="%(asctime)s:%(levelname)s:%(message)s",
                        filename=biolog_file
                        )
    extractor = DensityExtractor(schedule=schedule, 
                                 mode=mode, 
                                 window_size=window_size, 
                                 float_precision=args.float_precision,
                                 format=args.format,
                                 compression=args.compression,
                                 use_pyarrow=args.use_pyarrow,
                                 statistics=args.statistics,
                                 design=design,
                                 tempdir=args.tempdir)
#    densities_df = extractor.process(extraction=extraction, 
#                      gff_file=gff,
#                      compartment=compartment,
#                      partition_col=partition_col,
#                      return_df=True
#                    )
    extractor.process_bucket(bucket_id=bucket_id,
                            compartment=compartment,
                            out=out,
                            partition_col=args.partition_col,
                            partition_groups=args.partition_groups,
                            sleeping_time=args.sleeping_time,
                            polarity_mode=args.polarity_mode,
                        )
