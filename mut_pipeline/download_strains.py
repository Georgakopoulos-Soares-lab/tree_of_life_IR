def main():
    import argparse
    import csv
    from pathlib import Path
    import tempfile
    from subprocess import run
    from typing import Optional

    parser = argparse.ArgumentParser(
        description="Download bacterial strains from NCBI RefSeq."
    )
    parser.add_argument(
        "species_taxid",
        type=int,
        help="NCBI Taxonomy ID of the species.",
    )
    parser.add_argument(
        "--output_dir", type=str, help="Directory to save downloaded strains."
    )
    parser.add_argument(
        "--variant_file",
        type=str,
        help="File containing list of variant Taxonomy IDs to download.",
    )
    args = parser.parse_args()
    collection: Optional[list[str]] = None
    with open(args.variant_file, mode="r", encoding="UTF-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["species_taxid"] == str(args.species_taxid):
                collection = row["group"].split("&")
                break
    if collection is None:
        print(
            f"No variants found for species taxid {args.species_taxid} in {args.variant_file}."
        )
        return
    with tempfile.NamedTemporaryFile("w", delete=False) as tmp_file:
        for accession in collection:
            tmp_file.write(accession + "\n")
        tmp_file.flush()
        cmd = f"datasets download genome accession --inputfile {tmp_file.name} --include gff3,genome"
        run(cmd, shell=True, check=True)
    run(f"unzip -n ncbi_dataset.zip", shell=True, check=True)
    target = Path(f"genomes_{args.species_taxid}")
    target.mkdir(exist_ok=True)
    run(f"mv ncbi_dataset/data/**/*.fna {target}", shell=True, check=True)
    


if __name__ == "__main__":
    main()
