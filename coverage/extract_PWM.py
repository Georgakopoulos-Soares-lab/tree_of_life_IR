import pandas as pd
from pathlib import Path
import subprocess
import os

extract_name = lambda accession: Path(accession).name.split('.filtered')[0]
extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])

def extract_intersections(accession: os.PathLike[str], out: os.PathLike[str], compartment: str, gff_parent: os.PathLike[str], window_size: int = 500):
    accession = Path(accession).resolve()
    out = Path(out).resolve()
    out.mkdir(exist_ok=True, parents=True)
    
    gff_parent = Path(gff_parent).resolve()
    accession_id = extract_id(accession)
    gff_loc = gff_parent.joinpath(extract_name(accession) + '.gff.gz')
    
    # filter GFF file
    # GFF Headers:
    # chromosome, db, compartment, start, end, '_', 'strand', '_', ID

    gff_df = pd.read_csv(gff_loc, 
                         delimiter="\t", 
                         comment="#", 
                         header=None,
                         dtype={"start": int, "end": int},
                         usecols=[0, 2, 3, 4, 6, 8],
                         names=["chromosome", 
                                "compartment", 
                                "start", 
                                "end", 
                                "strand", 
                                "metadata"
                                ]
                    )
    gff_df = gff_df[gff_df['compartment'] == compartment]\
                    .drop(columns=['compartment'])

    gff_df["biotype"] = gff_df["metadata"].apply(lambda x: x.split("gene_biotype=")[1].split(";")[0] if "gene_biotype=" in x else "")
    gff_df["start"] = gff_df["start"] - 1
    gff_df["end"] = gff_df["end"] - 1

    positions = ["start", "end"]
    gff_df.drop(columns=["metadata"], inplace=True)
    
    # read extractions
    # tandem headers
    # chromosome   start     end  length  sru  consensus_repeats consensus              sequence

    extraction_bed_headers = ["chromosome", "start", "end", "consensus", "sequence"]
    df = pd.read_parquet(accession, engine="fastparquet").reset_index(drop=True)
    df["start"] = df["start"] - 1
   
    ext_bed_out = f"{out}/{accession_id}.ext.bed"
    df[extraction_bed_headers].to_csv(ext_bed_out, index=False, mode="w", sep="\t", header=None)

    for pos in positions:
        gff_df[f"{pos}_inf"] = (gff_df[pos] - window_size).apply(lambda x: max(x, 0))
        gff_df[f"{pos}_sup"] = gff_df[pos] + window_size + 1
    
        gff_bed_headers = ["chromosome", f"{pos}_inf", f"{pos}_sup", "strand"]
        gff_bed_out = f"{out}/{accession_id}.{pos}.bed"
        gff_df[gff_bed_headers].to_csv(gff_bed_out, mode="w", sep="\t", index=False, header=None)
    
        intersection_out = f"{out}/{accession_id}.intersect.{pos}.bed"
        
        gff_sorted_bed_out = f"{out}/{accession_id}.{pos}.sorted.bed"
        command = f"sort -k1,1 -k2,2n {gff_bed_out} > {gff_sorted_bed_out}"
        subprocess.run(command, shell=True, check=True)

        command = f"bedtools intersect -a {gff_sorted_bed_out} -b {ext_bed_out} -wo > {intersection_out}"
        subprocess.run(command, shell=True, check=True)
        
        os.remove(gff_bed_out)
        os.remove(gff_sorted_bed_out)
        
        # post-process intersection file
        intersectdf = pd.read_csv(intersection_out, 
                                        header=None, 
                                        delimiter="\t",
                                        usecols=[0, 2, 3, 5, 6, 7, 8, 9],
                                        dtype={"comp_end": int, "start": int, "end": int},
                                        names=[
                                               "chromosome", 
                                               "comp_end",
                                               "strand",
                                               "start",
                                               "end",
                                               "consensus",
                                               "sequence",
                                               "overlap",
                                              ]
                                )

        intersectdf["comp_end"] = intersectdf["comp_end"] - window_size - 1
        intersectdf.to_csv(intersection_out, mode="w", sep="\t", index=False, header=None)

    os.remove(ext_bed_out)


if __name__ == "__main__":
    
    accession = "/home/dollzeta/emergencyroom/localdb/GCF_000002725.2_ASM272v2_genomic.filtered.parquet.snappy"
    bed_out = "bed_out"
    extract_intersections(accession, out=bed_out)









