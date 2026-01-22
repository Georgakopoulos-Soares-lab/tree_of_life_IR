import pyranges as pr
import pandas as pd
from pathlib import Path
from termcolor import colored
from tqdm import tqdm
from sanitize import sanitize_df

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate total bases using pyranges for different GC thresholds, partitioned by spacer length")
    parser.add_argument("--indir", "-i", type=str, required=True, help="Input directory containing TSV files")
    parser.add_argument("--pattern", "-p", type=str, default="IR", choices=["IR", "STR", "MR"])
    parser.add_argument("--min_partition", "-m", type=int, help="Minimum partition value")
    parser.add_argument("--max_partition", "-x", type=int, help="Maximum partition value")

    args = parser.parse_args()
    min_partition = args.min_partition
    max_partition = args.max_partition
    if min_partition is None or max_partition is None:
        if args.pattern == "IR":
            partition_list = list(range(0, 9))
        elif args.pattern == "MR":
            partition_list = list(range(0, 8))
        elif args.pattern == "STR":
            partition_list = list(range(1, 10))
    else:
        partition_list = list(range(min_partition, max_partition + 1))
    partition_list += ["."]

    if args.pattern == "IR" or args.pattern == "MR":
        partition_col = "spacer_length"
    else:
        partition_col = "sru"

    indir = Path(args.indir)
    outdir = indir / f"summary_{args.pattern}"
    outdir.mkdir(exist_ok=True, parents=True)

    infiles = [infile for infile in indir.glob(f"*_genomic_{args.pattern}.processed.tsv")]
    print(colored(f"Total files: {len(infiles)}", "green"))
    results = []

    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    for infile in tqdm(infiles, desc="Processing files"):
        df = pd.read_csv(infile, sep='\t')
        accession_id = extract_id(infile)
        df = sanitize_df(df).to_pandas()
        # Filter by maximum spacer length if specified
        for partition in partition_list:
            if partition == ".":
                df_filtered = df.copy()
            else:
                df_filtered = df[df[partition_col] == partition].copy()

            if len(df_filtered) == 0:
                results.append({
                    '#assembly_accession': accession_id,
                    partition_col: partition,
                    'num_regions_before_merge': 0,
                    'num_regions_after_merge': 0,
                    f'total_bases_{args.pattern}': 0,
                    'avg_region_length': 0
                })
                continue

            # Overall statistics (all spacer lengths combined)
            gr_all = pr.PyRanges(
                chromosomes=df_filtered["seqID"],
                starts=df_filtered["start"],
                ends=df_filtered["end"]
            )
            gr_all_merged = gr_all.merge()
            total_bases_all = (gr_all_merged.df['End'] - gr_all_merged.df['Start']).sum()
            num_regions_all = len(gr_all_merged)

            results.append({
                '#assembly_accession': accession_id,
                 partition_col: partition,
                'num_regions_before_merge': len(df_filtered),
                'num_regions_after_merge': num_regions_all,
                f'total_bases_{args.pattern}': total_bases_all,
                'avg_region_length': total_bases_all / num_regions_all if num_regions_all > 0 else 0
            })

    # Save detailed results to CSV
    results_df = pd.DataFrame(results)
    assembly_df = pd.read_table(f"data/assembly_summary_with_tree.csv.gz")
    results_df = results_df.merge(
        assembly_df[["#assembly_accession",
                        "species_taxid",
                        "organism_name",
                        "total_gene_count",
                        "non_coding_gene_count",
                        "protein_coding_gene_count",
                        "genome_size",
                        "genome_size_ungapped",
                        "gc_percent",
                        "genus",
                        "family",
                        "order",
                        "class",
                        "phylum",
                        "kingdom",
                        "domain"]],
        on="#assembly_accession",
        how="left"
    )
    results_df.loc[:, f"density_{args.pattern}"] = 1e3 * results_df[f"total_bases_{args.pattern}"] / results_df["genome_size"]
    output_file = outdir / f"extractions_{args.pattern}_merged.tsv"
    results_df.to_csv(output_file, index=False, sep="\t")
    print(f"\nDetailed partition matrix saved to: {output_file}.")

if __name__ == "__main__":
    main()
