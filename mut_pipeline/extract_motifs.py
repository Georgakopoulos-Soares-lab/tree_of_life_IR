def main():
    import argparse
    from pathlib import Path

    from tqdm import tqdm

    from minditool import MindiTool

    parser = argparse.ArgumentParser(""".""")
    parser.add_argument("--genome", type=str, default="genomes")
    parser.add_argument("--pattern", type=str, default="STR", choices=["STR", "IR"])
    parser.add_argument("--species_taxid", type=int, default=None)
    args = parser.parse_args()
    genomes = Path(args.genome).resolve()
    fasta = [file for file in genomes.glob("*.fna")]
    if args.species_taxid:
        name = f"motifs_{args.pattern}_{args.species_taxid}"
    else:
        name = f"motifs_{args.pattern}"
    tempdir = Path(name)
    for accession in tqdm(fasta):
        t = MindiTool(tempdir=tempdir)
        t.extract(accession=accession, pattern=args.pattern)
        if args.pattern in t.fn:
            t.sanitize(accession=accession, mode=args.pattern)



if __name__ == "__main__":
    main()
