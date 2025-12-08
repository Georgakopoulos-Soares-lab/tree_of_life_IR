import pandas as pd 
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass, field
import polars as pl
import json
from gff_utils import CoverageExtractor
from termcolor import colored
from tqdm import tqdm
from typing import Optional
from scipy.stats import ks_2samp
from pwm_density import PWMExtractor

# plotting
import matplotlib.pyplot as plt
import seaborn as sns

@dataclass(frozen=True)
class params:
    N: int = field(default=1000)
    window_size: int = field(default=500)
    alpha: float = field(default=0.05)

def bootstrap_density(intersect_df: pd.DataFrame, 
                        window_size: int, 
                        nsamples: int = 1000, 
                        alpha: float = 0.05) -> tuple:
    df_bootstrap = []
    extractor = PWMExtractor()
    for _ in tqdm(range(nsamples), leave=True, position=0):
        df_sample = intersect_df.sample(frac=1.0, replace=True)
        df_sample = extractor.extract_density(df_sample, 
                                              window_size=window_size, 
                                              enrichment=True,
                                              return_array=True)
        df_bootstrap.append(df_sample)
    df_bootstrap = pd.DataFrame(df_bootstrap)
    mean = df_bootstrap.mean()
    ci_lower = df_bootstrap.quantile(alpha/2)
    ci_upper = df_bootstrap.quantile(1-alpha/2)
    return mean, ci_lower, ci_upper

def bootstrap(df: pl.DataFrame, N: int = 1_000, alpha: float = 0.05) -> tuple:
    bootstrapped_samples = []
    df = df.to_pandas()
    df: pd.DataFrame 
    for _ in tqdm(range(N), leave=True, position=0):
        sample_df = df.sample(frac=1.0, replace=True).sum()
        # poor performance
        # sample_df = df.sample(fraction=total_samples, with_replacement=True).sum()
        # sample_df = sample_df / sample_df.mean_horizontal()
        sample_df = sample_df / sample_df.mean()
        bootstrapped_samples.append(sample_df)
    bootstrapped_samples = pd.concat(bootstrapped_samples, axis=1).T
    average = bootstrapped_samples.mean()
    # two-tailed interval (1-a)%
    lower_ci = bootstrapped_samples.quantile(alpha/2)
    upper_ci = bootstrapped_samples.quantile(1 - alpha/2)
    return average, lower_ci, upper_ci

class Bootstrapper:
    
    def __init__(self, enrichment_file: str, 
                        design: str, 
                        params=params, 
                        biotypes: Optional[list[str]] = None,
                        bootstrap_levels: Optional[list[str]] = None,
                        taxonomic_ranks: Optional[list[str]] = None,
                        float_precision: int = 3) -> None:
        self.enrichment_df = None
        self.enrichment_file = Path(enrichment_file).resolve()
        self.design = Path(design).resolve()
        self.params = params
        self.float_precision = float_precision
        if not taxonomic_ranks:
            self.taxonomic_ranks = ["phylum", "kingdom", "domain"]
        else:
            self.taxonomic_ranks = taxonomic_ranks
        if not bootstrap_levels:
            self.bootstrap_levels = ["species_taxid", "family", "order", "class"]
        else:
            self.bootstrap_levels = bootstrap_levels
        if not biotypes:
            self.biotypes = ["protein_coding", "non_coding", "."]
        if not self.design.is_file():
            raise FileNotFoundError(f"Could not detect design file `{design}`.") 
        if not self.enrichment_file.is_file():
            raise FileNotFoundError(f"Could not detect enrichment file `{enrichment_file}`.") 

    def load_table(self, taxonomic_ranks: Optional[list[str]] = None):
        if self.enrichment_df is not None:
            return self
        if taxonomic_ranks is None:
            taxonomic_ranks = self.taxonomic_ranks
        if not isinstance(taxonomic_ranks, list):
            raise TypeError(f"Invalid type for taxonomic ranks. Expected list, but received {type(taxonomic_ranks)}.")

        delimiter = CoverageExtractor._sniff_delimiter(self.design)
        if not self.bootstrap_levels:
            columns = ["accession_id"] + taxonomic_ranks
        else:
            columns = ["accession_id"] + self.bootstrap_levels + taxonomic_ranks
        design_df = pl.read_csv(self.design, 
                                columns=columns,
                                separator=delimiter
                                )
        print(f"Loaded design file with columns: {columns}.")
        if "parquet" in self.enrichment_file.name:
            enrichment_df = pl.read_parquet(self.enrichment_file)
        else:
            enrichment_df = pl.read_csv(self.enrichment_file, separator="\t")
        enrichment_df = enrichment_df\
                        .join(
                                design_df,
                                right_on="accession_id",
                                left_on="#assembly_accession",
                                how="inner"
                              )
        self.enrichment_df = enrichment_df
        return self

    def bootstrap_enrichment(self, taxonomic_rank: str, 
                                   rank: str, 
                                   output: str, 
                                   polarity: bool = True,
                                   partition_col: Optional[str] = None,
                                   merge_species_taxid: bool = True,
                                   partition_list: Optional[list] = None,
                                  combinations: Optional[list[tuple]] = None) -> pl.DataFrame:
        if isinstance(partition_col, str) and len(partition_col) == 0:
            partition_col = None
        if partition_col is None:
            partition_list = None
        print(f"Enrichment file detected: '{self.enrichment_file}'.")
        enrichment_df = self.load_table().enrichment_df
        if taxonomic_rank not in enrichment_df:
            raise KeyError(f"Invalid specified taxonomic rank `{taxonomic_rank}`.")
        if (partition_col is None and partition_list is not None) or (partition_col is not None and partition_list is None):
            raise TypeError("Partition column and partition list variables must simultaneously be None.")
        enrichment_df = enrichment_df.filter(pl.col(taxonomic_rank) == rank)

        print(f"Initializing bootstrap for taxonomic rank {taxonomic_rank} with value=`{rank}`.")
        print(f"Specified confidence: {1e2 * (1-self.params.alpha):.2f}")
        print(f"Total resampling iterations: {self.params.N}")
        if enrichment_df.shape[0] == 0:
            raise ValueError(f"Empty dataframe for taxonomic rank `{taxonomic_rank}` with value=`{rank}`.")
        selected_columns = ["biotype"]
        if partition_col and "." not in partition_list:
            print(f"Chosen partition on column `{partition_col}'. Partition values: {partition_list}.")
            partition_list.append("generic")
            selected_columns.append(partition_col)
        else:
            print(f"No partition column detected. Won't perform any partitioning.")
            partition_list = [None]
        if polarity:
            polarities = ["Template", "Non-Template"]
            selected_columns.append("strand_polarity")
        else:
            polarities = [None]
        combinations = [(polarity, biotype, partition_value) for biotype in self.biotypes for polarity in polarities for partition_value in partition_list]
        ## Domain level bootstrap
        confidence_intervals = defaultdict(list)
        for comb in combinations: 
            polarity, biotype, partition_value = comb
            temp_df = enrichment_df.filter(pl.col("biotype") == biotype)
            if polarity:
                temp_df = temp_df.filter(pl.col("strand_polarity") == polarity)
            if partition_col:
                temp_df = temp_df.filter(pl.col(partition_col) == str(partition_value))

            for level in self.bootstrap_levels:
                print(f"Bootstrapping for level `{level}`.")
                if merge_species_taxid:
                    if level not in temp_df.columns:
                        raise KeyError(f"Level {level} doesn't exist.")
                    before = temp_df.shape[0]
                    temp_df_agg = (
                                temp_df
                                .group_by(level, maintain_order=True)
                                .agg(*[pl.col(str(i)).sum() for i in range(-params.window_size, params.window_size+1)])
                                .select([str(i) for i in range(-params.window_size, params.window_size+1)])
                            )
                    print(f"Merging species taxid. Total rows prior to reduction {before}. Total rows after reduction: {temp_df_agg.shape[0]}.")
                else:
                    temp_df_agg = temp_df.select([str(i) for i in range(-params.window_size, params.window_size+1)])
                average, lower_ci, upper_ci = bootstrap(temp_df_agg, N=params.N, alpha=params.alpha)
                average = average.tolist()
                lower_ci = lower_ci.tolist()
                upper_ci = upper_ci.tolist()
                biotype = biotype.replace("_coding", "-Coding")
                for i in range(2 * params.window_size+1):
                    confidence_intervals[str(i - params.window_size)].append(average[i])
                    confidence_intervals[str(i - params.window_size)].append(lower_ci[i])
                    confidence_intervals[str(i - params.window_size)].append(upper_ci[i])
                values = ["average", "lowerCI", "upperCI"]
                for i in range(3):
                    confidence_intervals["statistic"].append(values[i])
                    confidence_intervals["biotype"].append(biotype)
                    if polarity:
                        confidence_intervals["strand_polarity"].append(polarity)
                    if partition_value is not None:
                        confidence_intervals[partition_col].append(str(partition_value))
                    confidence_intervals["bootstrap_level"].append(level)
                    confidence_intervals[taxonomic_rank].append(rank)
                if not merge_species_taxid:
                    break

        confidence_intervals = pl.DataFrame(confidence_intervals)\
                                        .select([taxonomic_rank, "statistic", "bootstrap_level"] + selected_columns + list(map(str, range(-params.window_size, params.window_size+1))))
        confidence_intervals.write_csv(output, 
                                       separator=",", 
                                       include_header=True, 
                                       float_precision=self.float_precision)
        print(colored(f"Bootstrap has succesfully been completed for taxonomic rank {taxonomic_rank}=`{rank}`.", "green"))
        return confidence_intervals

    def average_phylums(self, domain: str, output: str, join_templates: int = 0) -> None:
        self.load_table()
        enrichment_df: pd.DataFrame = self.enrichment_df.to_pandas()
        if "domain" not in enrichment_df:
            raise KeyError(f"Invalid specified domain `{domain}`.")
        enrichment_df = enrichment_df.query(f"domain == '{domain}'")
        if "phylum" not in enrichment_df: 
            raise KeyError(f"No phylums detected for domain `{domain}`.")
        ## Phylum calculation
        ## keep assembly accessions with a present phylum
        enrichment_df_phylum = enrichment_df.dropna(subset=['phylum'])
        phylum_to_kingdom = dict(zip(enrichment_df_phylum["phylum"], enrichment_df_phylum["kingdom"]))

        if "strand_polarity" in enrichment_df_phylum and not join_templates:
            print(colored("Template & Non-template partition detected.", "blue"))
            enrichment_df_phylum = enrichment_df_phylum.groupby(["phylum", "strand_polarity", "biotype"])
        else:
            enrichment_df_phylum = enrichment_df_phylum.groupby(["phylum", "biotype"])
        # enrichment_df_phylum = enrichment_df_phylum.grouby(["species_taxid"])\
        #                                           .agg({str(i): "sum" for i in range(-params.window_size, params.window_size+1)})

        enrichment_df_phylum = enrichment_df_phylum.agg({str(i): "sum" for i in range(-params.window_size, params.window_size+1)})
        enrichment_df_phylum = enrichment_df_phylum.apply(lambda row: row / row.mean(), axis=1)
        for col in range(-params.window_size, params.window_size+1):
            enrichment_df_phylum[str(col)] = enrichment_df_phylum[str(col)].round(2)
        enrichment_df_phylum.reset_index(inplace=True)
        enrichment_df_phylum.loc[:, "kingdom"] = enrichment_df_phylum["phylum"].map(phylum_to_kingdom)
        enrichment_df_phylum.loc[:, "domain"] = domain
        enrichment_df_phylum.to_csv(output, sep=",", mode="w", header=True, index=True)
        print(colored(f"Phylum averaging has succesfully been completed for domain=`{domain}`.", "green"))
        return

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="""Utility""")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--N", type=int, default=1_000)
    parser.add_argument("--design", type=str, default="design.csv")
    parser.add_argument("--enrichment", type=str, default="")
    parser.add_argument("--polarity", type=int, choices=[0, 1], default=1)
    parser.add_argument("--partition_col", type=str, default=None)
    parser.add_argument("--partition_list", nargs="+", type=list, default=None)
    parser.add_argument("--float_precision", type=int, default=3)
    parser.add_argument("--taxonomic_rank", type=str, default="domain", choices=["domain", "kingdom", "phylum"])
    parser.add_argument("--output", type=str, default="bootstrap.txt")
    parser.add_argument("--window_size", type=int, default=500)
    parser.add_argument("--rank", type=str, default="Bacteria", choices=["Eukaryota", "Archaea", "Viruses", "Bacteria"])

    args = parser.parse_args()
    N = args.N 
    alpha = args.alpha
    design = args.design 
    rank = args.rank
    window_size = int(args.window_size)
    taxonomic_rank = args.taxonomic_rank
    enrichment_file = args.enrichment
    output = args.output
    polarity = args.polarity
    partition_col = args.partition_col
    partition_list = args.partition_list
    float_precision = args.float_precision

    param = params(alpha=alpha, N=N, window_size=window_size)
    bootstrapper = Bootstrapper(params=param, 
                                design=design, 
                                enrichment_file=enrichment_file,
                                float_precision=float_precision)
    bootstrapper.bootstrap_enrichment(taxonomic_rank=taxonomic_rank,
                                      rank=rank,
                                      polarity=polarity,
                                      partition_col=partition_col,
                                      partition_list=partition_list,
                                      output=output)
