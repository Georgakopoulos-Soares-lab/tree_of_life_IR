import pyranges as pr
import pandas as pd
from pathlib import Path
from tqdm import tqdm

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate total bases using pyranges for different GC thresholds")
    parser.add_argument("--indir", "-i", type=str, required=True, help="Input directory containing TSV files")
    parser.add_argument("--outdir", "-o", type=str, default="merged_output", help="Output directory")
    parser.add_argument("--pattern", "-p", type=str, default="IR")

    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    infiles = [infile for infile in indir.glob(f"*_genomic_{args.pattern}.processed.tsv")]
    print(f"Total files loaded: {len(infiles)}.")
    gc_thresholds = [0.0, 0.1, 0.2, 0.3]
    results = []
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    for infile in tqdm(infiles, desc="Processing files"):
        df = pd.read_csv(infile, sep='\t')

        # Calculate GC content of arms
        df["gc_arm"] = df["sequence_of_arm"].str.lower().str.count("[gc]") / df["arm_length"]

        for threshold in gc_thresholds:
            df_filtered = df[df["gc_arm"] >= threshold].copy()

            if len(df_filtered) == 0:
                total_bases = 0
                num_regions = 0
            else:
                gr = pr.PyRanges(
                    chromosomes=df_filtered["seqID"],
                    starts=df_filtered["start"],
                    ends=df_filtered["end"]
                )
                gr_merged = gr.merge()
                total_bases = (gr_merged.ends - gr_merged.starts).sum()
                num_regions = len(gr_merged)

            results.append({
                '#assembly_accession': extract_id(infile.name),
                'gc_threshold': threshold,
                'num_regions_before_merge': len(df_filtered),
                'num_regions_after_merge': num_regions,
                'total_bases': total_bases,
                'avg_region_length': total_bases / num_regions if num_regions > 0 else 0
            })

            print(f"{infile.name} - GC≥{threshold}: {num_regions} regions, {total_bases:,} total bases")

    # Save results to CSV
    results_df = pd.DataFrame(results)
    output_file = outdir / "merged_coords_summary.csv"
    results_df.to_csv(output_file, index=False)

    print(f"\nResults saved to: {output_file}")
if __name__ == "__main__":
    main()
