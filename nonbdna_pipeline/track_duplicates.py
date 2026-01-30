# NCBI Database statistics
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import json
from tqdm import tqdm
from termcolor import colored
import logging

def load_groups():
    assembly_groups = {}
    with open("assembly_groups.txt", mode="r", encoding="UTF-8") as f:
        for line in f:
            line = line.strip()
            accession_id, group = line.split("\t")
            assembly_groups.update({accession_id: group})
    return assembly_groups

def filter_duplicated_assemblies():
    logging.basicConfig(level=logging.INFO, 
                        format="%(levelname)s:%(asctime)s:%(message)s", 
                        filename="assembly_NIKOL.log")
    assemblies = [accession for accession in Path("files").resolve().glob("*.fna.gz")]

    extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])
    extract_unique_id = lambda accession_id: accession_id.split('.')[0].split("_")[1]
    extract_version = lambda accession_id: int(accession_id.split('.')[1])

    assembly_groups = load_groups()

    gff_files = {extract_id(gff): gff for gff in Path("files").resolve().glob("*.gff.gz")}
    seen_accessions = {}
    total_assemblies = len(assemblies)
    total_assemblies_replaced = 0

    for accession in tqdm(assemblies, total=total_assemblies, leave=True, position=0):

        accession_id = extract_id(accession)
        unique_id = extract_unique_id(accession_id)
        current_version = extract_version(accession_id)
        # replace GCA with GCF accession if both exist (reason: more enriched GFF features)

        if unique_id not in seen_accessions:
            seen_accessions.update({unique_id: str(accession)})
        else:
            # id already exists
            logging.info(f"Accession {accession_id} already exists.")
            existing_accession = seen_accessions[unique_id]
            existing_id = extract_id(existing_accession)
            old_version = extract_version(existing_id)

            # # # # # # # # # # # # # # # # # # # # # 
            # replace the existing accession with the new accession in one of the following scenarions:
            # - the existing accession is from genbank and the new accession is from refseq | Reason: refseq has more enriched GFF features
            # - the existing accession is from genbank and the new accession (either from genbank or from refseq) has a more latest version and the existing

            if existing_id.startswith("GCA_") and accession_id.startswith("GCF_"):
                logging.info(f"Replacing assembly {existing_id} with {accession_id} (Reason: GCA ---> GCF).")
                seen_accessions.update({unique_id: str(accession)})
                total_assemblies_replaced += 1
            elif existing_id.startswith("GCA_") and old_version < current_version:
                logging.info(f"Replacing assembly {existing_id} with {accession_id} (Reason: replacing latest version {old_version} ---> {current_version}).")
                seen_accessions.update({unique_id: str(accession)})
                total_assemblies_replaced += 1
            elif existing_id.startswith("GCF_") and accession_id.startswith("GCF_") and old_version < current_version:
                logging.info(f"Replacing assembly {existing_id} with {accession_id} (Reason: replacing latest version {old_version} ---> {current_version}).")
                seen_accessions.update({unique_id: str(accession)})
                total_assemblies_replaced += 1
            else:
                logging.info(f"Maintaning original assembly {accession_id}.")

    logging.info(f"Total assemblies: {total_assemblies}.")
    logging.info(f"Total filtered assemblies: {len(seen_accessions)}.")
    logging.info(f"Total replaced assemblies: {total_assemblies_replaced}.") 
    logging.info("Assocating assemblies with gff files...")

    assemblies_with_gffs = set()
    for unique_id in seen_accessions:
        full_path = seen_accessions[unique_id]
        accession_id = extract_id(full_path)
        if accession_id in gff_files:
            assemblies_with_gffs.add(unique_id)

    logging.info(f"Total {len(assemblies_with_gffs)} assemblies correspond to a GFF file.")
    now = datetime.now()
    year = now.year
    month = now.month
    day = now.day

    group_counts = defaultdict(int)
    with open(f"filtered_assemblies_{year}_{month}_{day}.txt", mode="w", encoding="UTF-8") as f:
        for unique_id, accession in seen_accessions.items():
            has_gff = unique_id in assemblies_with_gffs
            accession_id = extract_id(accession)
            group = assembly_groups[accession_id]
            group_counts[group] += 1
            f.write(f"{accession}\t{group}\t{int(has_gff)}\n")

    with open(f"db_statistics_{year}_{month}_{day}.txt", mode="w", encoding="UTF-8") as f:
        for group, counts in group_counts.items():
            f.write(f"{group}\t{counts}\n")

if __name__ == "__main__":
    filter_duplicated_assemblies()