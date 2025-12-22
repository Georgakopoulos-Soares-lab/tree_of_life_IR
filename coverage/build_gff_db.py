def main():
    import shutil
    import sys
    from pathlib import Path

    outdir = Path(sys.argv[2])
    outdir.mkdir(exist_ok=True)
    files = {file.parent.name: file for file in Path(sys.argv[1]).glob("**/*.gff")}
    fasta = {file.parent.name: file for file in Path(sys.argv[1]).glob("**/*.fna")}
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    for file_id, file in files.items():
        outfile = outdir / fasta[file_id].name.split(".fna")[0] + ".gff"
        shutil.move(file, outfile)
