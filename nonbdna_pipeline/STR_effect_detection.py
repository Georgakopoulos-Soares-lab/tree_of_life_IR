from pybedtools import BedTool 
from tqdm import tqdm 
from pathlib import Path
import pandas as pd
from utils import load_bucket
import csv

COVERAGE_FIELDS = ["total_hits", "overlapping_bp", "compartment_length", "coverage"]

if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("schedule", type=str, default="schedule.json")
    parser.add_argument("--pattern", type=str, default="IR", choices=["IR"])
    parser.add_argument("--outdir", type=str)
    parser.add_argument("--bucket_id", type=int)
    parser.add_argument("--partition", type=str, default="spacer_length")
    
    # # # # # # #
    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)

    outfile = outdir.joinpath(f"bucket_{args.bucket_id}_STR_effects.tsv")
    f = open(outfile, mode="w", encoding="utf-8")
    writer = csv.DictWriter(f, fieldnames=["#assembly_accession", 
                                                "total_IR_bp", 
                                                "total_STR_bp",
                                                "overlap_ratio"], 
                            delimiter="\t")
    pattern = args.pattern
    schedule = args.schedule
    bucket_id = args.bucket_id 
    partition_col = args.partition

    IR_parent = Path(f"extractions_{pattern}")
    STR_parent = Path("extractions_STR")
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    extract_name = lambda x: Path(x).name.split(".fna")[0]
    
    # Process Files
    files = load_bucket(schedule, bucket_id=bucket_id)
    for file in tqdm(files):
        file = Path(file)
        name = extract_name(file)
        IR_file = IR_parent / f"{name}_{pattern}.tsv.gz"
        STR_file = STR_parent / f"{name}_STR.tsv.gz"
        if not STR_file.is_file() or not IR_file.is_file():
            print(file)
            continue
        accession_id = extract_id(file)
        motif_df = pd.read_table(file, usecols=["seqID", "start", "end", partition_col])
        motif_bed = BedTool.from_dataframe(motif_df).sort().merge()
        STR_df = pd.read_table(STR_file, usecols=["seqID", "start", "end"])
        STR_bed = BedTool.from_dataframe(STR_df).sort().merge()
        # # # #  total
        motif_df = pd.read_table(
                            motif_bed.coverage(STR_bed).fn,
                            header=None,
                            names=["seqID", "start", "end"] + COVERAGE_FIELDS
                        )
        total_STR_bp = motif_df["overlapping_bp"].sum()
        total_bp = motif_df["compartment_length"].sum()
        overlap_ratio = round(total_STR_bp / total_bp, 2)
        writer.writerow({
                        "#assembly_accession": accession_id,
                        "total_bp": total_bp,
                        "total_STR_bp": total_STR_bp,
                        "overlap_ratio": overlap_ratio
                        })
    f.close()
        

