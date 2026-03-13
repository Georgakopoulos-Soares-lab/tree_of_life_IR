import time
import polars as pl 
from pathlib import Path 
import pybedtools
from termcolor import colored
from typing import Optional 
import pandas as pd
import numpy as np
import pyranges as pr
import threading 
import gzip
import csv
import logging
from typing import Optional 
import attr
from attr import field
from nonbdna_pipeline.stream_and_merge_bucket import StreamAndMerge 
from nonbdna_pipeline.pwm_density import PWMExtractor 

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

def read_gff(gff_file: str, 
             end_to_one_bp: bool = False,
             pseudogenes_to_genes: bool = True, 
             filter_on: Optional[str] = None,
             parse_biotype: bool = False,
             binary_biotype: bool = True,
             change_compartment_names: bool = True) -> pl.DataFrame:
    def _parse_attributes(attributes: str, 
                          attribute: str,
                          early_stop: bool = True,
                          ) -> str:
        attrs = attributes.split(";")
        attrs_dict = dict()
        for attr in attrs:
            if "=" in attr:
                key, val = attr.split("=", 1)
                attrs_dict[key] = val
                if early_stop and key == attribute:
                    break
        return attrs_dict.get(attribute, np.nan)
    def _map_biotype(biotype: str) -> str:
        # Treat missing or unknown biotypes as non_coding so genes that
        # don't explicitly state protein_coding are classified as non_coding
        if biotype == ".":
            return "."
        if pd.isna(biotype):
            return "non_coding"
        if biotype == "protein_coding":
            return "protein_coding"
        if biotype == "pseudogene":
            return "pseudogene"
        return "non_coding"
    GFF_FIELDS = ["seqID", "source", "compartment", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.read_table(gff_file, 
                       comment="#", 
                       header=None, 
                       names=GFF_FIELDS, 
                       dtype={"start": np.int32, 
                              "end": np.int32})
    df = df[df["start"] < df["end"]].reset_index(drop=True)
    if df.shape[0] == 0:
        return df
    df.loc[:, "start"] = df["start"] - 1
    if end_to_one_bp:
        df.loc[:, "end"] = df["end"] - 1
    if pseudogenes_to_genes:
        df["compartment"] = df["compartment"].replace({"pseudogene": "gene"})
    if filter_on:
        if isinstance(filter_on, str):
            df = df[df["compartment"] == filter_on]
        else:
            filter_on = set(filter_on)
            df = df[df["compartment"].isin(filter_on)]
        df = df.reset_index(drop=True)
    if df.shape[0] == 0:
        return df
    selected_columns = ["seqID", "start", "end", "compartment", "strand"]
    if parse_biotype:
        df.loc[df["compartment"] != "gene", "biotype"] = "."
        df.loc[df["compartment"] == "gene", "biotype"] = df.loc[df["compartment"] == "gene", "attributes"].apply(lambda y: _parse_attributes(y, attribute="gene_biotype"))
        df = df[~df["biotype"].isna()].reset_index(drop=True)
        df["biotype"] = df["biotype"].astype("string")
        selected_columns.append("biotype")
        if binary_biotype:
            df["biotype"] = df["biotype"].apply(_map_biotype)
    if change_compartment_names:
        compartment_mapping = {
            "exon": "Exon",
            "intron": "Intron",
            "pseudogene": "Pseudogene",
            "gene": "Gene",
            "five_prime_UTR": "5'UTR",
            "three_prime_UTR": "3'UTR",
        }
        for compartment, mapped_name in compartment_mapping.items():
            df["compartment"] = df["compartment"].replace(compartment, mapped_name)
    return df[selected_columns]

def expand_gff(gff_df, on: str, 
               window_size: int, 
               ignore_biotype: bool = False) -> pd.DataFrame:
    positive_charge = gff_df[gff_df["strand"] == "+"].copy()
    negative_charge = gff_df[gff_df["strand"] == "-"].copy()
    if on == "TSS":
        positive_charge.loc[:, "expanded_start"] = np.maximum(0, positive_charge["start"] - window_size)
        positive_charge.loc[:, "expanded_end"] = positive_charge["start"] + window_size + 1

        negative_charge.loc[:, "expanded_start"] = np.maximum(0, negative_charge["end"] - window_size)
        negative_charge.loc[:, "expanded_end"] = negative_charge["end"] + window_size + 1
    elif on == "TES":
        positive_charge.loc[:, "expanded_start"] = np.maximum(0, positive_charge["end"] - window_size)
        positive_charge.loc[:, "expanded_end"] = positive_charge["end"] + window_size + 1

        negative_charge.loc[:, "expanded_start"] = np.maximum(0, negative_charge["start"] - window_size)
        negative_charge.loc[:, "expanded_end"] = negative_charge["start"] + window_size + 1
    else:
        raise ValueError(f"Invalid `on` value `{on}`. Must be one of 'TSS' or 'TES'.")
    gff_df = pd.concat([positive_charge, negative_charge], ignore_index=True)
    selected_cols = ["seqID", "expanded_start", "expanded_end", "compartment", "strand"]
    if "biotype" in gff_df.columns and not ignore_biotype:
        selected_cols.append("biotype")
    gff_df = gff_df[selected_cols].rename(columns={"seqID": "Chromosome", 
                                                    "expanded_start": "Start", 
                                                    "expanded_end": "End",
                                                    "strand": "Strand"})
    return gff_df 

def calculate_strand_polarity(df, pattern: str, THRESHOLD: float = 0.8) -> pd.DataFrame:
    df["arm_length"] = df["sequence_of_arm"].str.len()
    if pattern == "HDNA":
        df["ga_proportion"] = df["sequence_of_arm"].str.count("[ga]") / df["arm_length"]
        # df["ct_proportion"] = df["sequence_of_arm"].str.count("[ct]") / df["arm_length"]
        df["Strand"] = np.where(df["ga_proportion"] >= THRESHOLD, "+", "-")
    elif pattern == "GT":
        df["gt_proportion"] = df["sequence_of_arm"].str.count("[gt]") / df["arm_length"]
        # df["ca_proportion"] = df["sequence_of_arm"].str.count("[ca]") / df["arm_length"]
        df["Strand"] = np.where(df["gt_proportion"] >= THRESHOLD, "+", "-")
    elif pattern == "G4":
        df["g_proportion"] = df["sequence_of_arm"].str.count("g") / df["arm_length"]
        df["c_proportion"] = df["sequence_of_arm"].str.count("c") / df["arm_length"]
        df["Strand"] = np.where(df["g_proportion"] >= df["c_proportion"], "+", "-")
    elif pattern == "IR":
        # Does not really make sense; for testing purposes
        df["gc_proportion"] = df["sequence_of_arm"].str.count("[gc]") / df["arm_length"]
        df["Strand"] = np.where(df["gc_proportion"] >= 0.5, "+", "-")
    else:
        raise ValueError(f"Strand polarity calculation not implemented for pattern `{pattern}`.")
    return df

@attr.s
class TSSTESProcessor(StreamAndMerge):
    window_size: int = field(converter=int, default=500)
    files_processed: int = field(init=False, default=0)
    log_dir: Path = field(init=False)
    total_files: int = field(init=False, default=0)
    sites: list[str] = field(init=False, factory=lambda: ["TSS", "TES"])
    outdir: Path = field(init=False)
    polarities: list[str] = field(init=False, factory=lambda: ["Template", "Non-Template"])
    biotypes: list[str] = field(init=False, factory=lambda: ["protein_coding", "non_coding", "."])
    LOG_INTERVAL: int = 240
    FIELDS: list[str] = ["#assembly_accession", "pattern", "site", "biotype", "polarity", "partition", "overlapping_genes", "pct_gene", "total_genes"]
    def __attrs_post_init__(self) -> None:
        super().__attrs_post_init__()
        self.outdir = self.indir.joinpath("tss_tes_density")
        self.outdir.mkdir(exist_ok=True, parents=True)
        self.log_dir = self.log_indir.joinpath("tss_tes_density_logs")
        self.log_dir.mkdir(exist_ok=True, parents=False)
        return
    def _setup_logging(self, bucket_id: int) -> None:
        DATE = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            filename=self.log_dir.joinpath(f"tss_tes_processing_{DATE}_{bucket_id}.log"),
            # handlers=[
            #     logging.FileHandler(self.log_dir.joinpath(f"tss_tes_processing_{DATE}.log")),
            #     logging.StreamHandler()
            # ]
        )
    def _log_progress(self) -> None:
        while True:
            prc = self.files_processed / self.total_files * 1e2 if self.total_files > 0 else 0.0
            logging.info(f"Progress: {prc:.2f}% files processed.")
            time.sleep(TSSTESProcessor.LOG_INTERVAL)
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

    def process_bucket(self, bucket_id: int, 
                       gff_indir: str, 
                       pattern: str, 
                       partition_col: str, 
                       polarity: bool = False,
                       min_partition: Optional[int] = None, 
                       max_partition: Optional[int] = None,
                       pseudogenes_to_genes: bool = True,
                       use_biotype: bool = False,
                       assembly_summary: Optional[str] = None,
                       ignore_errors: bool = False) -> None:

        if use_biotype:
            biotypes = tuple(self.biotypes)
        else:
            biotypes = (".", )
        infiles = self.load_validated_files_from_log(bucket_id=bucket_id, pattern=pattern, ignore_errors=ignore_errors)
        print(colored(f"Total files to process in bucket {bucket_id}: {len(infiles)}.", "green"))
        # extract_name = lambda x: Path(x).name.split(".fna")[0]
        gff_indir = Path(gff_indir).resolve()
        if min_partition is None or max_partition is None:
            min_partition, max_partition = TSSTESProcessor._get_min_max_partition(pattern)
        if not gff_indir.is_dir():
            raise ValueError(f"Invalid directory `{gff_indir}`.")
        # Dummy for now
        partitions = ["."] 
        pwm = PWMExtractor()
        self.total_files = len(infiles)
        self._setup_logging(bucket_id=bucket_id)
        if polarity:
            logging.info(f"Strand polarity calculation enabled for pattern `{pattern}`.")
            print(colored(f"Strand polarity calculation enabled for pattern `{pattern}`.", "yellow"))
        logging.info(f"Initializing processing of {self.total_files} files in bucket {bucket_id}.")
        thread = threading.Thread(target=self._log_progress, daemon=True)
        thread.start()
        
        # Process initialization
        outfile = self.outdir.joinpath(f"tss_tes_density_{pattern}_bucket_{bucket_id}.tsv.gz")
        fin = gzip.open(outfile, mode="wt")
        writer = csv.DictWriter(fin, fieldnames=TSSTESProcessor.FIELDS + list(map(str, range(-self.window_size, self.window_size+1))), delimiter="\t")
        window_range = list(map(str, range(-self.window_size, self.window_size+1)))
        writer.writeheader()
        for file_idx, infile in enumerate(infiles, start=1):
            accession_name = StreamAndMerge.extract_name(infile, pattern=pattern)
            accession_id = StreamAndMerge.extract_id(infile)
            extraction_file = self.indir.joinpath(accession_name + f"_{pattern}.processed.tsv")
            gff_file = gff_indir.joinpath(accession_name + ".gff")
            if not gff_file.is_file():
                gff_file = gff_file.with_suffix(".gff.gz")
                if not gff_file.is_file():
                    continue 
            gff_df = read_gff(gff_file, 
                              end_to_one_bp=True, 
                              pseudogenes_to_genes=pseudogenes_to_genes,
                              parse_biotype=use_biotype,
                              filter_on="gene")
            total_genes = {".": gff_df.shape[0]}
            if use_biotype:
                total_genes.update({"protein_coding": gff_df[gff_df["biotype"] == "protein_coding"].shape[0],
                                    "non_coding": gff_df[gff_df["biotype"] == "non_coding"].shape[0]})
            if total_genes["."] == 0:
                logging.warning(f"No features found in GFF file `{gff_file}`. Skipping accession `{accession_id}`.")
                continue 

            df = self.read_motifs(extraction_file)
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
            if polarity:
                df = calculate_strand_polarity(df, pattern=pattern)
                df_gr = pr.PyRanges(df)
                if not df_gr.stranded:
                    raise ValueError(f"Invalid data: Strand polarity calculation enabled but motifs are not stranded in file `{extraction_file}`.")
                df_gr = df_gr.merge(strand=True)
            else:
                df_gr = pr.PyRanges(df)
                # if df_gr.stranded:
                #    logging.warning(f"Strand polarity calculation disabled but motifs are stranded in file `{extraction_file}`.")
                df_gr = df_gr.merge(strand=False)

            for site in self.sites:
                for biotype in biotypes:
                    if biotype != ".":
                        gff_biotype_df = gff_df[gff_df["biotype"] == biotype].copy()
                    else:
                        gff_biotype_df = gff_df.copy()
                    if gff_biotype_df.shape[0] == 0:
                        logging.warning(f"No features found for biotype `{biotype}` in GFF file `{gff_file}`. Skipping.")
                        continue
                    expanded_gff_df = expand_gff(gff_biotype_df, on=site, window_size=self.window_size)
                    gff_gr = pr.PyRanges(expanded_gff_df)
                    df_joined = (
                                    gff_gr.join(df_gr)
                                          .as_df()
                                          .rename(columns={
                                                            "Chromosome": "seqID",
                                                            "Start": "start",
                                                            "End": "end",
                                                            "Strand": "strand",
                                                            "Start_b": "motif_start",
                                                            "End_b": "motif_end",
                                                            "Strand_b": "motif_strand"}
                                          )
                            )
                    # In that case, directly emit zero vectors for this site across all partitions.
                    if df_joined.shape[0] == 0:
                        zero_vec = {locus: 0 for locus in window_range}
                        for partition in partitions:
                            base_data = {
                                "#assembly_accession": accession_id,
                                "pattern": pattern,
                                "site": site,
                                "biotype": biotype,
                                "partition": str(partition),
                                "overlapping_genes": 0,
                                "pct_gene": 0.0,
                                "total_genes": total_genes[biotype]
                            }
                            if polarity:
                                for charge in self.polarities:
                                    writer.writerow(base_data | {"polarity": charge} | zero_vec)
                            else:
                                writer.writerow(base_data | {"polarity": "."} | zero_vec)
                        continue
                    assert df.shape[0] > 0, f"No motifs found after joining with GFF for file `{extraction_file}` and site `{site}`."
                    df_joined["overlap"] = (np.minimum(df_joined["end"], df_joined["motif_end"]) - np.maximum(df_joined["start"], df_joined["motif_start"])).clip(lower=0)

                    if polarity:
                        df_joined["strand_polarity"] = np.where(df_joined["strand"] == df_joined["motif_strand"], "Non-Template", "Template")

                    for partition in partitions:
                        if partition != ".":
                            df_partitioned = df_joined[df_joined[partition_col] == partition]
                        else:
                            df_partitioned = df_joined 
                        df_partitioned = pl.from_pandas(df_partitioned)
                        if polarity:
                            density_df = pwm.extract_template_density(df_partitioned, window_size=self.window_size, return_frame=False)
                        else:
                            density_df = pwm.extract_density(df_partitioned, window_size=self.window_size)

                        gene_overlap = (
                            df_partitioned
                            .select(["seqID", "start", "end"])
                            .unique()
                            .height
                        )
                        pct_gene_overlap = round(1e2 * gene_overlap / total_genes[biotype], 2) if total_genes[biotype] > 0 else np.nan
                        data = {
                            "#assembly_accession": accession_id,
                            "pattern": pattern,
                            "site": site,
                            "biotype": biotype,
                            "partition": str(partition),
                            "overlapping_genes": gene_overlap,
                            "pct_gene": pct_gene_overlap,
                            "total_genes": total_genes[biotype],
                        }
                        if polarity:
                            for charge in self.polarities:
                                data.update({"polarity": charge} | {locus: int(counts) for locus, counts in zip(window_range, density_df[charge])})
                                writer.writerow(data)
                        else:
                            data.update({"polarity": "."} | {locus: int(counts) for locus, counts in zip(window_range, density_df)})
                            writer.writerow(data)
            self.files_processed = file_idx
        thread.join(timeout=1)
        fin.close()
        density_df = merge_with_summary(assembly_summary=assembly_summary,
                                          outfile=outfile)
        if density_df is not None:
            merged_outfile = self.outdir.joinpath(f"tss_tes_density_{pattern}_bucket_{bucket_id}_with_assembly_data.tsv.gz")
            density_df.to_csv(merged_outfile, mode="w", sep="\t", index=False, compression="gzip")
            logging.info(f"Merged density data with assembly summary and saved to `{merged_outfile}`.")
        logging.info(f"Process has been completed succesfully (bucket {bucket_id}).")
        return

def main():
    import argparse
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("schedule", type=str)
    parser.add_argument("-i", "--indir", type=str, required=True)
    parser.add_argument("--gff_indir", type=str, required=True)
    parser.add_argument("--window_size", type=int, default=500)
    parser.add_argument("--strand_polarity", "-s", default=0, type=int, choices=[0,1], help="Whether to calculate strand polarity of motifs.")
    parser.add_argument("--bucket_id", "-bid", type=int, default=0)
    parser.add_argument("-p", "--pattern", type=str, default='IR', choices=['IR', 'MR', 'STR'])
    parser.add_argument("--partition_col", type=str, default=None)
    parser.add_argument("--use_biotype", "-b", action="store_true", default=False)
    parser.add_argument("--tmpdir", type=str, default="garbage")
    parser.add_argument("--assembly_summary", "-asm", type=str, 
                        default="data/assembly_summary_with_tree.csv.gz", 
                        help="Path to assembly summary file to merge with density data.")
    parser.add_argument("--ignore_errors", action="store_true", help="Whether to ignore errors when loading validated files from log.")
    args = parser.parse_args()
    tmpdir = Path(args.tmpdir)
    tmpdir.mkdir(exist_ok=True)
    pybedtools.helpers.set_tempdir(tmpdir)

    TSSTESProcessor(indir=args.indir, 
                    schedule=args.schedule,
                    window_size=args.window_size).process_bucket(bucket_id=args.bucket_id,
                                                           pattern=args.pattern,
                                                           gff_indir=args.gff_indir,
                                                           polarity=args.strand_polarity,
                                                           partition_col=args.partition_col,
                                                           assembly_summary=args.assembly_summary,
                                                           ignore_errors=args.ignore_errors,
                                                           use_biotype=args.use_biotype,
                                                           )
    pybedtools.helpers.cleanup(remove_all=True)
if __name__ == "__main__": main()
