# example
# accession_id    gff     extraction      group   phylum  kingdom superkingdom    genome_size
# GCA_000002515.1 /storage/group/izg5139/default/external/gff/GCA_000002515.1_ASM251v1_genomic.gff        /storage/group/izg5139/default/external/quadrupia_database/g4/gquadruplex/GCA_000002515.1_ASM251v1_genomic.csv     fungi   Ascomycota      Fungi   Eukaryota       10689156

if __name__ == "__main__":
    import sys
    from pathlib import Path
    import argparse
    from termcolor import colored
    from tqdm import tqdm
    import pandas as pd
    from collections import defaultdict
    # # # # #
    extract_id = lambda accession: "_".join(Path(accession).name.split("_")[:2])
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("indir", type=str)
    parser.add_argument("--gff_files", type=str)
    parser.add_argument("--tree_of_life", type=str, default="../nonbdna_pipeline/assembly_summary_with_tree.csv.gz")
    parser.add_argument("--pattern", type=str, default="STR", choices=["STR", "MR", "HDNA", "GT", "IR", "DR", "GQ"])
    # # # # 
    args = parser.parse_args()
    indir = Path(args.indir).resolve()
    assert indir.is_dir()
    pattern = args.pattern

    # Load Tree of Life
    TREE_OF_LIFE = ["#assembly_accession", "species_taxid", "family", "order", "class", "group", "phylum", "kingdom", "domain", "genome_size"]
    assembly_df = pd.read_table(args.tree_of_life, usecols=TREE_OF_LIFE)

    # Load GFF Files
    GFF_PATH = Path(args.gff_files).resolve()
    GFF_files = defaultdict(list)
    total_gff = 0
    with open(GFF_PATH, mode="r", encoding="UTF-8") as f:
        for line in f:
            total_gff += 1
    with open(GFF_PATH, mode="r", encoding="UTF-8") as f:
        for line in tqdm(f, total=total_gff, leave=True):
            line = line.replace("\"", "").replace(",", "").strip()
            if Path(line).is_file():
                accession_id = extract_id(line)
                GFF_files["accession_id"].append(accession_id)
                GFF_files["gff"].append(line)
    print(f"Total GFF files loaded: {len(GFF_files)}.")
    GFF_files = pd.DataFrame(GFF_files)\
                    .merge(
                                assembly_df,
                                how="left",
                                left_on="accession_id",
                                right_on="#assembly_accession"
                            )
    print(f"Sourcing from {indir}...")
    files = {extract_id(file): file for file in indir.glob(f"*_{pattern}.processed.tsv")}
    print(f"Total files found: {len(files)}.")
    print("Creating preliminary schedule...")
    with open(f"design_{pattern}.csv", mode="w", encoding="UTF-8") as f:
        f.write("accession_id\textraction\n")
        for file_id, file in files.items():
            f.write(f"{file_id}\t{file}\n")
    print("Creating enhanced design file...")
    design_df = pd.read_table(f"design_{pattern}.csv")\
                    .merge(GFF_files,
                            how="inner",
                            on="accession_id")
    design_df.to_csv(f"design_{pattern}_with_tree.csv", mode="w", sep="\t", index=False, header=True)
    print(colored("DONE!", "green"))
