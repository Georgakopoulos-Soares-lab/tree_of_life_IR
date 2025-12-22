def main():
    import argparse
    import pandas as pd
    from pathlib import Path
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("--design_file", type=str)
    parser.add_argument("--success", type=str)
    args = parser.parse_args()
    success = set()
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    with open(args.success, mode="r") as f:
        for line in f:
            line = line.strip()
            success.add(extract_id(line))
    design_df = pd.read_table(args.design_file)
    design_df = design_df[design_df["accession_id"].isin(success)]
    design_renamed = Path(args.design_file.replace(".csv", ".filtered.csv"))
    design_df.to_csv(design_renamed, mode="w", sep="\t", index=False, header=True)

if __name__ == "__main__": main()

