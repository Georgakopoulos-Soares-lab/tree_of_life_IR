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
    parser.add_argument("--indir", type=str)
    parser.add_argument("--outdir", type=str)
    parser.add_argument("--bucket_id", type=int)
    parser.add_argument("--partition", type=str, default="spacer_length")
    
    # # # # # # #
    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)
    indir = Path(args.indir).resolve()
    level = 2

    outfile = outdir.joinpath(f"bucket_{args.bucket_id}_STR_effects.tsv")
    f = open(outfile, mode="w", encoding="utf-8")
    writer = csv.DictWriter(f, fieldnames=["#assembly_accession", 
                                           "partition",
                                            "total_IR_bp", 
                                            "total_STR_bp",
                                            "overlap_ratio"], 
                            delimiter="\t")
    pattern = args.pattern
    schedule = args.schedule
    bucket_id = args.bucket_id 
    partition_col = args.partition
    IR_parent = indir / Path("extractions_shuffled_level_12-12-2025_2_IR")
    STR_parent = indir / Path("extractions_shuffled_level_12-13-2025_2_STR")
    assert IR_parent.is_dir()
    assert STR_parent.is_dir()
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    extract_name = lambda x: Path(x).name.split("_level")[0]

    # Process Files
    files = load_bucket(schedule, bucket_id=bucket_id)
    partition_list = list(range(0, 9)) + ["."]
    for file in tqdm(files):
        file = Path(file)
        name = extract_name(file)
        # example GCA_000007325.1_ASM732v1_genomic_level_2.shuffled_IR.processed.tsv
        IR_file = IR_parent / f"{name}_level_{level}.shuffled_IR.processed.tsv"
        STR_file = STR_parent / f"{name}_level_{level}.shuffled_STR.processed.tsv"
        if not STR_file.is_file() or not IR_file.is_file():
            print(file)
            continue
        ####
        accession_id = extract_id(file)
        motif_df = pd.read_table(IR_file, usecols=["seqID", "start", "end", partition_col])
        motif_bed = BedTool.from_dataframe(motif_df).sort() # .merge()
        STR_df = pd.read_table(STR_file, usecols=["seqID", "start", "end"])
        STR_bed = BedTool.from_dataframe(STR_df).sort().merge()
        # # # #  total
        motif_df = pd.read_table(
                            motif_bed.coverage(STR_bed).fn,
                            header=None,
                            names=["seqID", "start", "end", partition_col] + COVERAGE_FIELDS
                        )
        for partition in partition_list:
            if partition != ".":
                temp_df = motif_df[motif_df[partition_col] == partition].copy()
            else:
                temp_df = motif_df.copy()

            total_STR_bp = temp_df["overlapping_bp"].sum()
            total_bp = temp_df["compartment_length"].sum()
            if total_bp == 0:
                overlap_ratio = None
            else:
                overlap_ratio = round(1e2 * total_STR_bp / total_bp, 2)
            writer.writerow({
                        "#assembly_accession": accession_id,
                        "partition": partition,
                        "total_IR_bp": total_bp,
                        "total_STR_bp": total_STR_bp,
                        "overlap_ratio": overlap_ratio
                        })
    f.close()
        

