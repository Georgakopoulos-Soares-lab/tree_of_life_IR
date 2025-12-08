from pathlib import Path
from utils import parse_fasta
import sys
from tqdm import tqdm
import pandas as pd
genomes = [file for file in Path(sys.argv[1]).resolve().glob("*.fna")]
for infile in genomes:
    df = infile.parent.joinpath(infile.name.replace(".fna", "_STR.processed.tsv"))
    print(infile, df)
    data = pd.read_table(df)
    for seqID, seq in parse_fasta(infile):
        temp = data[data["seqID"] == seqID]
        for _, row in tqdm(temp.iterrows(), total=temp.shape[0]):
            start = int(row["start"])
            end = int(row["end"]) 
            sequence = row["sequence"]
            length = int(row["sequence_length"])
            original = seq[start: end]
            arm = row["sequence_of_arm"]
            sru = int(row["sru"])
            repeating = int(row["consensus_repeats"])
            assert length == len(original) == len(sequence)
            assert original == sequence == arm * repeating
        
    
