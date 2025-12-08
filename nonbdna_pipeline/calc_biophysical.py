import polars as pl

if __name__ == "__main__":

    # Script that calculates densities from taxonomic ranks
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description="""""")
    parser.add_argument("--dest", "-d", type=str, default="density")
    parser.add_argument("--source", "-s", type=str, required=True)
    parser.add_argument("--mode", "-m", type=str, required=True, 
                        choices=["IR", "MR", "DR", "STR", "Z", "G4"]
                        )
    # DO NOT CHANGE THIS!
    parser.add_argument("--MAX_ARM_LENGTH", type=int, default=31)

    args = parser.parse_args()
    source = args.source
    mode = args.mode
    MAX_ARM_LENGTH = args.MAX_ARM_LENGTH
    dest = Path(args.dest).resolve().joinpath("biophysical")
    dest.mkdir(exist_ok=True, parents=True)
    print("CREATED DIRECTORY `biophysical`...")

    # load not-empty assemblies
    df = pl.read_parquet(source)
    df = df.with_columns(
                pl.col("sequenceOfArm")
                    .str.to_lowercase()
                    .str.replace_all(r"[^gc]", "")
                    .str.len_chars()
                    .alias("gc_content")
                    )\
            .with_columns(
                ( (pl.col("gc_content") / pl.col("armLength")) / pl.col("gc_percent"))\
                        .alias("gc_enrichment")
            )

    # calculate gc enrichment in comparison to background
    # partition on arm length & domain
    df_grouped = df.filter(df["armLength"] < MAX_ARM_LENGTH)\
                .group_by(["armLength", "superkingdom"], maintain_order=True)\
                .agg([
                        pl.col("gc_enrichment").mean().round(2),
                        pl.col("gc_content").sum(),
                        pl.col("start").count().alias("total_occurrences")
                    ])\
                .with_columns(
                        (1e2 * pl.col("gc_content") / (pl.col("total_occurrences") * pl.col("gc_content")))\
                                .alias("gc_total_percentage")
                        )
    df_grouped.write_csv(f"{dest}/biophysical_arm_{mode}_gc_enrichment.csv", 
                         separator="\t")

    # calculate gc enrichment in comparison to background
    # partition on spacer length & domain
    df_grouped = df.group_by(["spacerLength", "superkingdom"], maintain_order=True)\
                    .agg([
                            pl.col("gc_enrichment").mean().round(2)
                        ])
    df_grouped.write_csv(f"{dest}/biophysical_spacer_{mode}_gc_enrichment.csv", separator="\t")

    # calculate arm length % per domain
    domain_total = df.group_by("superkingdom", maintain_order=True)\
                            .agg([
                                pl.col("start").count().alias("total_occurrences")
                                ])

    df_grouped = df.group_by(["armLength", "superkingdom"], maintain_order=True)\
                        .agg([
                            pl.col("start").count().alias("occurrences")
                            ])\
                        .join(domain_total, on="superkingdom", how="left")\
                        .with_columns(
                                (pl.col("occurrences") * 1e2 / pl.col("total_occurrences"))\
                                            .round(2)\
                                            .alias("arm_percentage")
                            )
    df_grouped.write_csv(f"{dest}/biophysical_arm_occurrences_{mode}.csv", separator="\t")

    # calculate spacer length % per domain
    df_grouped = df.group_by(["spacerLength", "superkingdom"], maintain_order=True)\
                        .agg([
                            pl.col("start").count().alias("occurrences")
                            ])\
                        .join(domain_total, on="superkingdom", how="left")\
                        .with_columns(
                                (pl.col("occurrences") * 1e2 / pl.col("total_occurrences"))\
                                        .round(2)\
                                        .alias("spacer_percentage")
                            )
    df_grouped.write_csv(f"{dest}/biophysical_spacer_occurrences_{mode}.csv", separator="\t")

    # calculate spacer length & arm length population per domain
    domains = ["Eukaryota", "Bacteria", "Viruses", "Archaea"]
    df_grouped = df.filter(pl.col("armLength") < MAX_ARM_LENGTH)\
                        .group_by(["armLength", "spacerLength", "superkingdom"], maintain_order=True)\
                        .agg([
                                pl.col("start").count().alias("occurrences")
                            ])\
                        .join(domain_total, on="superkingdom", how="left")\
                        .with_columns(
                                (pl.col("occurrences") * 1e2 / pl.col("total_occurrences"))\
                                        .round(2)\
                                        .alias("percentage")
                        )\
                        .sort(by=["superkingdom", "percentage"])
    df_grouped.write_csv(f"{dest}/biophysical_arm_spacer_prevalence_{mode}.csv", separator="\t")