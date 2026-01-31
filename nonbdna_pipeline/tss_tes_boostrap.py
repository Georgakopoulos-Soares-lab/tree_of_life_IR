import pandas as pd
from pathlib import Path 
from termcolor import colored 
import gzip

def bootstrap_sample_density(density_df: pd.DataFrame, 
                             rank: str, 
                             taxonomy: str,
                             bootstrap_taxonomic_level: str = "species",
                             n_samples: int = 1000, 
                             window_size: int = 500,
                             alpha: float = 0.05) -> dict[str, pd.DataFrame]:
    """
    Perform bootstrap sampling on density dataframe to calculate enrichment confidence intervals.
    """
    if rank not in density_df:
        raise KeyError(f"Taxonomic rank `{rank}` not found in density dataframe columns.")
    if bootstrap_taxonomic_level not in density_df:
        raise KeyError(f"Taxonomic level `{bootstrap_taxonomic_level}` not found in density dataframe columns.")
    density_df = density_df[density_df[rank] == taxonomy].reset_index(drop=True)
    RANGE = list(map(str, range(-window_size, window_size+1)))
    # Calculate enrichment across window
    window_occurrences = density_df[RANGE].sum(axis=1)
    window_average = window_occurrences.mean()
    enrichment_across_window = (window_occurrences / window_average).round(3)

    # bootstrap on taxonomic level 
    density_df = density_df.groupby(bootstrap_taxonomic_level).agg({loci: "sum" for loci in RANGE})
    bootstrap_results = dict()
    bootstrap_df = []
    for _ in range(n_samples):
        sample_df = density_df.sample(frac=1.0, replace=True)
        sample_window_occurrences = sample_df[RANGE].sum(axis=1)
        window_average = window_occurrences.mean()
        sample_window_enrichment = sample_window_occurrences / window_average
        bootstrap_df.append(sample_window_enrichment)

    # Compute Confidence Interval (quantile-based)
    bootstrap_df = pd.concat(bootstrap_df)
    enrichment_ci_lower = bootstrap_df.quantile(alpha/2).round(3)
    enrichment_ci_upper = bootstrap_df.quantile(1-alpha/2).round(3)
    bootstrap_results = {
        "average": enrichment_across_window,
        "ci_lower": enrichment_ci_lower,
        "ci_upper": enrichment_ci_upper
    }
    return bootstrap_results

def main():
    import argparse 
    parser = argparse.ArgumentParser()
    parser.add_argument("density_file", type=str)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--n_samples", type=int, default=1000)
    parser.add_argument("--window_size", type=int, default=500)
    parser.add_argument("--rank", type=str, default="domain")
    parser.add_argument("--taxonomy", type=str, default="Bacteria")
    parser.add_argument("--bootstrap_taxonomic_level", type=str, default="family")
    args = parser.parse_args()

    density_file = Path(args.density_file)
    if not density_file.is_file():
        raise FileNotFoundError(f"Density file `{density_file}` not found.")
    density_df = pd.read_table(density_file)
    parent = density_file.parent.joinpath("bootstrap_results")
    parent.mkdir(exist_ok=True)
    bootstrap_results = bootstrap_sample_density(
        density_df=density_df,
        rank=args.rank,
        taxonomy=args.taxonomy,
        bootstrap_taxonomic_level=args.bootstrap_taxonomic_level,
        n_samples=args.n_samples,
        window_size=args.window_size,
        alpha=args.alpha
    )
    outfile = parent.joinpath(f"{density_file.stem}_{args.rank}_{args.taxonomy}_bootstrap_summary.tsv")
    RANGE = list(map(str, range(-args.window_size, args.window_size+1)))
    with gzip.open(outfile, "rt") as fout:
        for key, df in bootstrap_results.items():
            occurrences = df[RANGE].to_csv(fout, mode="a", sep="\t", header=True, index=False)

if __name__ == "__main__": main()