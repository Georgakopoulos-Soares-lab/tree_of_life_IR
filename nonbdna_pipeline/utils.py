from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path
import gzip
from typing import Iterator
import json

def parse_fasta(file: str) -> Iterator[tuple[str, str]]:
    if Path(file).name.endswith(".gz"):
        f = gzip.open(file, "rt")
    else:
        f = open(file, mode="r", encoding="utf-8")
    for record in SimpleFastaParser(f):
        seqID = record[0]
        if " " in seqID:
            seqID = seqID.split(" ")[0]
        seq = record[1]
        yield seqID, seq.strip().lower()
    f.close()

def load_schedule(schedule: str) -> dict[int, list[str]]:
    with open(schedule, mode="r", encoding="utf-8") as f:
        return json.load(f)
    
def load_bucket_ids(SCHEDULE: str) -> list[int]:
    bucket_ids = []
    with open(SCHEDULE, mode="r", encoding="UTF-8") as f:
        schedule = json.load(f)
        for bucket_id in schedule:
            bucket_ids.append(int(bucket_id))
    return bucket_ids

def load_bucket(schedule: str, bucket_id: int) -> list[str]:
    with open(schedule,  mode="r", encoding="utf-8") as f:
        return json.load(f)[str(bucket_id)]
