from utils import parse_fasta
import json
from tqdm import tqdm
from pathlib import Path
primates = {
            "mPonPyg2": "P. pygmaeus",
            "mPanTro3": "P. troglodytes",
            "mGorGor1": "G. gorilla",
            "chm13v2": "H. sapiens",
            "mPanPan1": "P. paniscus",
            "mPonAbe1": "P. abelii",
            "mSymSyn1": "S. syndactylus",
    }
if __name__ == "__main__":
    import sys
    with open(sys.argv[1], mode="r") as f:
        data = json.load(f)
    genome_sizes = dict()
    for bid, file in tqdm(data.items()):
        file = file[0]
        species = primates.get(Path(file).name.split('.')[0], None)
        if species is None:
            continue
        genome_size = 0
        for seqID, seq in parse_fasta(file):
            genome_size += len(seq)
        genome_sizes[species] = genome_size
    with open("primate_genome_sizes.tsv", mode="w") as f:
        for species, gsize in genome_sizes.items():
            f.write(f"{species}\t{gsize}\n")
