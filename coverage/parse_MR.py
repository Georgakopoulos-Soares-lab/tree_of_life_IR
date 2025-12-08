from pathlib import Path
import pandas as pd
from utils import load_bucket
from tqdm import tqdm
import re

def max_STR_coverage(seq: str) -> float:
    """
    Returns the STR that coverages the most of sequence.
    """
    seq = seq.upper()
    str_lengths = range(1, 9)
    f = lambda n: 2 if n >= 3 else (3 if n == 2 else 7)
    patterns = [r"(([ATGC]{%s}?)\2{%s,})" % (n, f(n)) for n in str_lengths]
    patterns.append(r"(([ATGC]{9,}?)\2{2,})")  # {9,}
    coverage = [0 for _ in range(len(seq))]
    for pattern in patterns:
        matches = re.finditer(pattern, seq)
        for match in matches:
            start = match.start()
            end = match.end()
            motif = match.group(1)
            full_match = match.group(2)
            for pos in range(start, end):
                coverage[pos] = 1
    return min(4, (1e2 * sum(coverage) / len(seq)) // 20)

nucleotides = {"a", "g", "c", "t"}

def parse_HDNA(infile: str, outdir: str = "extractions_MR", threshold: float = 0.9) -> pd.DataFrame:
    df = pd.read_table(infile)
    HDNA_outfile = outdir.joinpath("HDNA", infile.name.replace("_MR.processed.tsv", "_HDNA.processed.tsv"))
    GT_outfile = outdir.joinpath("GT", infile.name.replace("_MR.processed.tsv", "_GT.processed.tsv"))
    MR_outfile = outdir.joinpath("MR", infile.name)
    if df.shape[0] > 0:
        df.loc[:, "sequence_of_arm"] = df["sequence_of_arm"].str.lower()
        df.loc[:, "sequence_of_spacer"] = df["sequence_of_spacer"].str.lower()
        df.loc[:, "sequence"] = df["sequence"].str.lower()
        df = df[df["sequence"].apply(lambda seq: all(c in nucleotides for c in seq))].reset_index(drop=True)
        ###
        df.loc[:, "AT_content"] = df["sequence_of_arm"].str.count("a|t").div(df["arm_length"])
        df.loc[:, "GA_content"] = df["sequence_of_arm"].str.count("g|a").div(df["arm_length"])
        df.loc[:, "GT_content"] = df["sequence_of_arm"].str.count("g|t").div(df["arm_length"])
        df.loc[:, "is_HDNA"] = ((df["AT_content"] < 0.8) & ((df["GA_content"] >= threshold) | (df["GA_content"] < 1 - threshold))).astype(int)
        df.loc[:, "is_GT"] = ((df["AT_content"] < 0.8) & ((df["GT_content"] >= threshold) | (df["GT_content"] < 1 - threshold))).astype(int)
        df.loc[:, "STR_coverage_sequence"] = df["sequence"].apply(max_STR_coverage).astype(int)
        # df.loc[:, "coverage_quantile"] = (df["STR_coverage_sequence"] // 20).astype(int)
        # df.loc[:, "STR_coverage_spacer"] = df["sequence_of_spacer"].apply(detect_STR_coverage)
        # df.loc[:, "STR_coverage_arm"] = df["sequence_of_arm"].apply(detect_STR_coverage)

        HDNA_df = df[df["is_HDNA"] == 1].copy()
        HDNA_df.loc[:, "strand"] = HDNA_df["GA_content"].apply(lambda x: "+" if x >= threshold else "-")
        GT_df = df[df["is_GT"] == 1].copy()
        GT_df.loc[:, "strand"] = GT_df["GT_content"].apply(lambda x: "+" if x >= threshold else "-")
    else:
        df.loc[:, "STR_coverage_sequence"] = []
        # df.loc[:, "coverage_quantile"] = []
        HDNA_df = df
        GT_df = df
    HDNA_df.to_csv(HDNA_outfile, sep="\t", index=False, header=True)
    GT_df.to_csv(GT_outfile, sep="\t", index=False, header=True)
    df.to_csv(MR_outfile, sep="\t", index=False, header=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("schedule", type=str)
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--outdir", type=str, default="/scratch/nmc6088/extractions_HDNA")

    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)
    outdir.joinpath("HDNA").mkdir(exist_ok=True)
    outdir.joinpath("GT").mkdir(exist_ok=True)
    outdir.joinpath("MR").mkdir(exist_ok=True)
    outdir.mkdir(exist_ok=True)
    schedule = args.schedule
    threshold = args.threshold
    data = load_bucket(schedule, bucket_id=args.bucket_id)
    for infile in tqdm(data):
        infile = Path(infile).resolve()
        parse_HDNA(infile, outdir=outdir, threshold=args.threshold)
