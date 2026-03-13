from pathlib import Path 
from termcolor import colored 
from typing import Optional
import gzip
import numpy as np
import pandas as pd
from tqdm import tqdm
# from nonbdna_pipeline.tss_tes_processing import merge_with_summary

def bootstrap_sample_density(density_df: pd.DataFrame, 
                             rank: str, 
                             taxonomy: str,
                             site: str,
                             assembly_summary: str,
                             bootstrap_taxonomic_level: str = "species",
                             biotype: str = ".",
                             n_samples: int = 1000, 
                             window_size: int = 500,
                             polarity: str = ".",
                             save_output: bool = True,
                             seed: int = 42,
                             density_file: Optional[str | Path] = None,
                             alpha: float = 0.05) -> dict[str, pd.Series]:
    """
    Perform bootstrap sampling on density dataframe to calculate enrichment confidence intervals.
    """
    for col in ["biotype", "polarity"]:
        if col not in density_df:
            raise KeyError(f"`{col}` column not found in density dataframe columns.")
    assembly_summary = Path(assembly_summary).resolve()
    if not assembly_summary.is_file():
        raise FileNotFoundError(f"Invalid assembly summary file `{assembly_summary}`. File does not exist!")
    summary_df = pd.read_table(assembly_summary, usecols=["#assembly_accession", rank, bootstrap_taxonomic_level])
    density_df = density_df.merge(summary_df, 
                                  on="#assembly_accession", 
                                  how="inner")
    if rank not in density_df:
        raise KeyError(f"Taxonomic rank `{rank}` not found in density dataframe columns.")
    if bootstrap_taxonomic_level not in density_df:
        raise KeyError(f"Taxonomic level `{bootstrap_taxonomic_level}` not found in density dataframe columns.")
    # When biotype is "all" or ".", select the "." rows which represent the
    # unfiltered total.  This avoids double-counting with per-biotype rows
    # (protein_coding, non_coding) that coexist when use_biotype is enabled.
    biotype_filter_value = "." if biotype in ("all", ".") else biotype
    density_df = density_df[
        (density_df[rank] == taxonomy)
        & (density_df["site"] == site)
        & (density_df["polarity"] == polarity)
        & (density_df["biotype"] == biotype_filter_value)
    ].reset_index(drop=True)
    if density_df.shape[0] == 0:
        raise ValueError(f"Empty dataframe!")
    RANGE = list(map(str, range(-window_size, window_size+1)))
    # Calculate enrichment across window
    window_occurrences = density_df[RANGE].sum(axis=0)
    window_average = window_occurrences.mean()
    enrichment_across_window = (window_occurrences / window_average).round(3)

    # bootstrap on taxonomic level 
    density_df = density_df.groupby(bootstrap_taxonomic_level).agg({loci: "sum" for loci in RANGE})
    bootstrap_results = dict()
    bootstrap_df = []
    rng = np.random.RandomState(seed)
    for _ in tqdm(range(n_samples), leave=True):
        iter_seed = rng.randint(0, 2**32 - 1)
        sample_df = density_df.sample(frac=1.0, replace=True, random_state=iter_seed)
        sample_window_occurrences = sample_df[RANGE].sum(axis=0)
        window_average = sample_window_occurrences.mean()
        sample_window_enrichment = pd.Series((sample_window_occurrences / window_average).to_numpy())
        bootstrap_df.append(sample_window_enrichment)

    # Compute Confidence Interval (quantile-based)
    bootstrap_df = pd.concat(bootstrap_df, axis=1).T
    bootstrap_df.columns = RANGE
    enrichment_ci_lower = bootstrap_df.quantile(alpha/2, axis=0).round(3)
    enrichment_ci_upper = bootstrap_df.quantile(1-alpha/2, axis=0).round(3)
    bootstrap_results = {
        "average": enrichment_across_window,
        "ci_lower": enrichment_ci_lower,
        "ci_upper": enrichment_ci_upper
    }
    if save_output and (isinstance(density_file, str) or isinstance(density_file, Path)):
        density_file = Path(density_file)
        parent_dir = density_file.parent.joinpath("bootstrap_results")
        parent_dir.mkdir(parents=True, exist_ok=True)
        # Map 'all' and '.' biotypes to 'all' in filename for consistency
        biotype_tag = "all" if biotype in ("all", ".") else biotype
        outfile = parent_dir.joinpath(f"{density_file.name.split('_all_')[0]}_{rank}_{taxonomy}_{bootstrap_taxonomic_level}_bootstrap_summary.{biotype_tag}.{site}.tsv.gz")
        print(colored(f"Saving results at ... {outfile}", "green"))
        with gzip.open(outfile, "wt") as fout:
            fout.write(f"site\t{rank}\tbiotype\tpolarity\tmethod\t" + "\t".join(RANGE) + "\n")
            for (method, df) in bootstrap_results.items():
                fout.write(f"{site}\t{taxonomy}\t{biotype_tag}\t{polarity}\t{method}\t" + "\t".join(str(round(val, 3)) for val in df[RANGE].values) + "\n")
        print(colored("DONE!", "green"))
    return bootstrap_results

def main():
    import argparse 
    parser = argparse.ArgumentParser()
    parser.add_argument("density_file", type=str)
    parser.add_argument("--site", "-s", type=str, choices=["TSS", "TES"], default="TSS")
    parser.add_argument("--alpha", "-a", type=float, default=0.05)
    parser.add_argument("--n_samples", "-n", type=int, default=1000)
    parser.add_argument("--window_size", "-w", type=int, default=500)
    parser.add_argument("--rank", "-r", type=str, default="domain")
    parser.add_argument("--taxonomy", "-t", type=str, default="Bacteria")
    parser.add_argument("--polarity", "-p", type=str, default=".")
    parser.add_argument("--assembly_summary", "-asm", type=str, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--save_output", "-o", action="store_true")
    parser.add_argument("--biotype", "-b", type=str, default="all", 
                        choices=["pseudogene", 
                                 "protein_coding", 
                                 "rRNA", 
                                 "non_coding", 
                                 "tRNA",
                                 "all", 
                                 "snRNA"]
                                 )
    parser.add_argument("--bootstrap_taxonomic_level", "-l", type=str, default="family")
    args = parser.parse_args()

    density_file = Path(args.density_file)
    if not density_file.is_file():
        raise FileNotFoundError(f"Density file `{density_file}` not found.")
    density_df = pd.read_table(density_file)
    bootstrap_results = bootstrap_sample_density(
        density_df=density_df,
        rank=args.rank,
        site=args.site,
        taxonomy=args.taxonomy,
        bootstrap_taxonomic_level=args.bootstrap_taxonomic_level,
        biotype=args.biotype,
        polarity=args.polarity,
        n_samples=args.n_samples,
        window_size=args.window_size,
        seed=args.seed,
        alpha=args.alpha,
        density_file=density_file,
        save_output=args.save_output,
        assembly_summary=args.assembly_summary,
    )

if __name__ == "__main__": main()