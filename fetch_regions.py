import sys
from pathlib import Path
from utils import parse_fasta
from Bio.Seq import Seq

to_keep = {"GCF_900639715.1": ("NZ_LR135490.1", "+"),
	   "GCF_900639415.1": ("NZ_LR135245.1", "+"),
	   "GCF_900635415.1": ("NZ_LR134107.1", "+"),
}
extract_id = lambda x: '_'.join(Path(x).name.split('_')[:2])
to_use = "NZ_LR135490.1"
accessions = {extract_id(file): file for file in Path().cwd().glob("*.fna") if not file.name.endswith("1.fna")}
print(accessions)
f = open("genomes.txt", "w")
f.write("#file\tname\ttags\n")
id = ["A", "B", "C"]
j = 0
for i, (accession_id, (seq_id, strand)) in enumerate(to_keep.items()):
    accession_df = accessions[accession_id]
    for seqID, seq in parse_fasta(accession_df):
        name = f"{accession_id}_{seqID}.fna"
        if seqID == seq_id:
            print(seqID, accession_id, j, i)
            print(accession_df, seqID, strand)
            if strand == "-":
                print("CHANGED")
                seq = str(Seq(seq).reverse_complement())
            with open(name, "w") as g:
                g.write(f">{to_use}\n{seq}\n")
            f.write(f"{name}\t{id[j]}\tlw:1.5\n")
            j += 1
            break
        
f.close()

