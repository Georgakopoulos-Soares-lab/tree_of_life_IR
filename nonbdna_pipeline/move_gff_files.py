def main():
    from pathlib import Path 
    import argparse
    import shutil
    parser = argparse.ArgumentParser()
    parser.add_argument("indir")
    parser.add_argument("-o", "--outdir", default="./gff_db")
    args = parser.parse_args()
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True, parents=False)
    genome_infiles = {infile.parent.name: infile for infile in indir.glob("**/*.fna")}
    infiles = {infile.parent.name: infile for infile in indir.glob("**/*.gff")}
    for accession_id, infile in infiles.items():
        genome_name = genome_infiles[accession_id].name.split(".fna")[0]
        gff_destination = outdir.joinpath(genome_name + ".gff")
        shutil.copy(infile, gff_destination)
if __name__ == "__main__": main()