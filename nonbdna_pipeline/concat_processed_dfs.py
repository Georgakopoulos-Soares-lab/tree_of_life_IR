def main():
    import argparse
    import pandas as pd
    from pathlib import Path
    import os
    from tqdm import tqdm
    parser = argparse.ArgumentParser()
    parser.add_argument("indir", type=str)
    parser.add_argument("--pattern", type=str, default="IR", choices=["IR", "STR"])
    args = parser.parse_args()
    indir = Path(args.indir).resolve()
    pattern = args.pattern
    if not indir.is_dir():
        raise ValueError(f"Directory `{indir}` does not exist.")
    infiles = [infile for infile in indir.glob(f"*_{pattern}.processed.tsv")]
    result = indir / f"{pattern}_biophysical.processed.all.tsv"
    if Path(result).is_file():
        os.remove(result)
    with result.open("a+") as fout:
        for i, infile in tqdm(enumerate(infiles), total=len(infiles)):
            df = pd.read_table(infile)
            if df.shape[0] == 0:
                continue
            df.to_csv(fout, header=i==0, sep="\t", index=False)
if __name__ == "__main__": main()
