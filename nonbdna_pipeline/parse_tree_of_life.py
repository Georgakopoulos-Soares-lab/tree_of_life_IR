import pandas as pd

TREE_COLS = ["tax_id", "species", "genus", "family", "order", "class", "phylum", "subkingdom", "kingdom", "domain"]
ASSEMBLY_COLS = ["#assembly_accession", "species_taxid", "taxid", 
                 "group", "genome_size", "genome_size_ungapped", 
                 "organism_name", "isolate", "infraspecific_name",
                 "gc_percent", 
                 "total_gene_count", "protein_coding_gene_count", "non_coding_gene_count",
                 "assembly_level", "assembly_type", "ftp_path",
                 ]
if __name__ == "__main__":

    import argparse 
    parser = argparse.ArgumentParser(description="Parse a tree of life file and print the parsed data.")
    parser.add_argument("--assembly_summary", type=str, default="assembly_summary.txt.gz", help="Path to the assembly summary file.")
    parser.add_argument("--tree", type=str, default="ncbi_lineages_2025-07-12.csv.gz", help="Path to the tree of life file.")

    args = parser.parse_args()

    tree_df = pd.read_csv(args.tree, usecols=TREE_COLS) 
    headers = pd.read_table("headers.txt").columns.tolist()
    assembly_df = pd.read_table(args.assembly_summary, header=None, names=headers, usecols=ASSEMBLY_COLS)\
                    .query("assembly_level == 'Complete Genome'")

    assembly_df_tree = assembly_df\
                    .merge(
                                tree_df,
                                left_on="species_taxid",
                                right_on="tax_id",
                                how="left"
                    )\
                    .drop_duplicates(subset=["#assembly_accession"])\
                    .reset_index(drop=True)

    null_ids = set(assembly_df_tree["taxid"][assembly_df_tree["tax_id"].isna()].unique())
    if len(null_ids) >  0:
        print(f"Warning: {len(null_ids)} taxids in assembly summary not found in tree of life. These will be set to 'unknown'.")
    if False:
        assembly_df_tree_n = assembly_df\
                        .merge(
                                    tree_df[tree_df["tax_id"].isin(null_ids)],
                                    left_on="taxid",
                                    right_on="tax_id",
                                    how="inner"
                            )\
                    .drop_duplicates(subset=["#assembly_accession"])\
                    .reset_index(drop=True)
        assembly_df = pd.concat([assembly_df_tree, assembly_df_tree_n], ignore_index=True)\
                    .drop_duplicates(subset=["#assembly_accession"])\
                    .reset_index(drop=True)
        null_ids = set(assembly_df["taxid"][assembly_df["tax_id"].isna()].unique())
        if len(null_ids) > 0:
            print(f"Warning: {len(null_ids)} taxids in assembly summary not found in tree of life. These will be set to 'unknown'.")
    assembly_df = assembly_df_tree

    # Replace values based on groups
    assembly_df.loc[assembly_df["group"] == "vertebrate_other", "kingdom"] = "Metazoa"
    assembly_df.loc[assembly_df["group"] == "vertebrate_other", "domain"] = "Eukaryota"
    assembly_df.loc[assembly_df["group"] == "vertebrate_mammalian", "kingdom"] = "Metazoa"
    assembly_df.loc[assembly_df["group"] == "vertebrate_mammalian", "domain"] = "Eukaryota"
    assembly_df.loc[assembly_df["group"] == "protozoa", "kingdom"] = "Protista"
    assembly_df.loc[assembly_df["group"] == "protozoa", "domain"] = "Eukaryota"
    assembly_df.loc[assembly_df["group"] == "plant", "kingdom"] = "Plantae"
    assembly_df.loc[assembly_df["group"] == "plant", "domain"] = "Eukaryota"
    assembly_df.loc[assembly_df["group"] == "fungi", "kingdom"] = "Fungi"
    assembly_df.loc[assembly_df["group"] == "fungi", "domain"] = "Eukaryota"
    assembly_df.loc[assembly_df["group"] == "bacteria", "domain"] = "Bacteria"
    assembly_df.loc[assembly_df["group"] == "archaea", "domain"] = "Archaea"
    assembly_df.loc[assembly_df["group"] == "viral", "domain"] = "Viruses"
    assembly_df.loc[:, "kingdom"] = assembly_df["kingdom"].replace("Viridiplantae", "Plantae")
    # Print total NaN
    print(assembly_df[assembly_df["domain"].isna()]["group"].value_counts())
    print(assembly_df[assembly_df["kingdom"].isna()]["group"].value_counts())
    assembly_df.to_csv("assembly_summary_with_tree.csv", index=False, sep="\t")
