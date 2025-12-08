import json
from pathlib import Path
import numpy as np
import pandas as pd
import subprocess
from typing import Optional
import os

extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
extract_idv2 = lambda accession: '.'.join(Path(accession).name.split('.')[:2])
correct_strands = {"+", "-"}

def process_enrichment(s_sites: pd.DataFrame, e_sites: pd.DataFrame, assembly, window_size: int) -> list[dict[str, int]]:
    global correct_strands
    counts = {}
    sequences = {
                 ("start", "+"): s_sites.query("strand=='+'"), 
                 ("start", "-"): e_sites.query("strand=='-'"),
                 ("end", "+"): e_sites.query("strand=='+'"),
                 ("end", "-"): s_sites.query("strand=='-'"),
            }

    from_overlap = s_sites[s_sites['strand'].isin(correct_strands)]['overlap'].sum() + e_sites[e_sites['strand'].isin(correct_strands)]['overlap'].sum()
    total_overlap = 0
    bins = {}
    for (ts, strand) in sequences:
        annotations = sequences[ts, strand]
        overlap = annotations['overlap'].sum()
        total_overlap += overlap

        bin_counts = construct_bin_counts(annotations, window_size, strand)
        assert overlap == sum(sum(counts for cell in array_counts for counts in cell.values()) for array_counts in bin_counts.values())
        bins.update({(ts, strand): bin_counts})

    for position in ["start", "end"]:
        counts[position] = {sru: 
                        [
                          {
                            'A': strpos['A'] + strneg['A'],
                            'T': strpos['T'] + strneg['T'],
                            'G': strpos['G'] + strneg['G'],
                            'C': strpos['C'] + strneg['C'],
                            } 
                          for strpos, strneg in zip(bins[position, "+"][sru], bins[position, "-"][sru][::-1])
                        ] \
                            for sru in range(1, 10, 1)}


    total_start_sum = sum(sum(val for cell in array_counts for val in cell.values()) for array_counts in counts["start"].values())
    total_end_sum = sum(sum(val for cell in array_counts for val in cell.values()) for array_counts in counts["end"].values())

    assert total_overlap == from_overlap == total_start_sum + total_end_sum, f"Invalid overlap detected in {assembly}. Expected: {total_overlap}, {from_overlap}, vs. {total_start_sum+total_end_sum}, s_sites: {s_sites.shape[0]}, e_sites: {e_sites.shape[0]}."
    
    return counts

def reduce_enrichment(files: list[os.PathLike[str]], 
                                window_size: int, 
                                compartment: str, 
                                out: os.PathLike[str], 
                                bucket_id: Optional[str]=None) -> None:
    out = Path(out).resolve()
    start_sites = {extract_idv2(file): file for file in files if file.name.split(".")[-2] == "start" and file.name.endswith(".bed")}
    end_sites = {extract_idv2(file): file for file in files if file.name.split(".")[-2] == "end" and file.name.endswith(".bed")}

    # overlap row
    # NC_001905.3	3203	4204	gene	protein_coding	NC_001905.3	3254	3270	ac	acacacacacacacac	16
    headers = [
               "chromosome", 
               "comp_start",
               "strand",
               "start", 
               "end", 
               "consensus", 
               "sequence",
               "overlap",
            ]
    
    assembly_ids = start_sites.keys()
    for assembly in assembly_ids:

        start_site = start_sites[assembly]
        end_site = end_sites[assembly]

        s_sites = pd.read_csv(start_site, delimiter="\t", header=None, names=headers)
        e_sites = pd.read_csv(end_site, delimiter="\t", header=None, names=headers)

        # biotypes = set(s_sites['biotype'])
        # split for different biotypes perphaps # TODO
        counts = process_enrichment(s_sites, e_sites, assembly, window_size)
        for transcription_site in ["start", "end"]:
            target = out.joinpath(f"density_{transcription_site}_{bucket_id}_{compartment}_{window_size}").with_suffix(".enrichment.csv")
            distribution = counts[transcription_site]
            if not target.is_file():
                with target.open("w", encoding="UTF-8") as f:
                    f.write("#assembly_accession,sru,nucleotide," + ",".join(str(i) for i in range(-window_size, window_size+1)) + "\n")

            with target.open("a", encoding="UTF-8") as f:
                assembly_id = extract_idv2(start_site)
                for sru, density in distribution.items():
                    for nucleotide in ['A', 'G', 'C', 'T']:
                        nucleotide_stream = [stats[nucleotide] for stats in density]
                        f.write(f"{assembly_id},{sru},{nucleotide}," + ",".join(str(occurrences) for occurrences in nucleotide_stream) + "\n")


def invert(nucleotide: str) -> str:
    match nucleotide:
        case 'A':
            return 'T'
        case 'G':
            return 'C'
        case 'C':
            return 'G'
        case 'T':
            return 'A'
        case 'N':
            return 'N'
        case _:
            raise ValueError(f"Unknown nucleotide {nucleotide}.")

invert_sequence = lambda sequence: ''.join(map(invert, sequence))

def construct_bin_counts(annotations, window_size: int, strand: str) -> list[dict]:
    assert strand == "+" or strand == "-", f"Invalid strand {strand}."
    stranded_counts = {sru: [{'A':0, 'G':0, 'C':0, 'T':0} for _ in range(-window_size, window_size+1)] for sru in range(1,10,1)}

    # construct bin counts
    annotations['sequence'] = annotations['sequence'].str.upper()
    for _, pattern in annotations.iterrows():

        sequence = pattern['sequence']
        strand = pattern['strand']
        sru = len(pattern['consensus'])
        
        if strand == "-":
            sequence = invert_sequence(sequence)

        LOWER = window_size - (pattern['comp_start'] - pattern['start'])
        UPPER = window_size - (pattern['comp_start'] - (pattern['end']-1))

        minimum_seq = max(pattern['comp_start']-window_size, pattern['start']) - pattern['start']
        maximum_seq = min(pattern['comp_start']+window_size, pattern['end']-1) - pattern['start']
        overlapping_chunk = sequence[minimum_seq: maximum_seq+1]
        assert minimum_seq <= maximum_seq, f"Invalid sequence length {len(sequence)}. UPPER: {UPPER}, LOWER: {LOWER}, Overlap: {overlapping_chunk}, LOWER: {minimum_seq}, UPPER: {maximum_seq}."

        if UPPER < 0 or LOWER > window_size * 2:
            raise ValueError("Out of bounds.")
        
        UPPER = min(UPPER, 2 * window_size)
        LOWER = max(LOWER, 0)

        try:
            assert LOWER <= UPPER, f"Invalid sequence length {len(sequence)}. UPPER: {UPPER}, LOWER: {LOWER}, Overlap: {overlapping_chunk}, LOWER: {minimum_seq}, UPPER: {maximum_seq}."
            assert len(sequence) >= UPPER - LOWER + 1, f"Invalid sequence length {len(sequence)}. UPPER: {UPPER}, LOWER: {LOWER}, Overlap: {overlapping_chunk}, LOWER: {minimum_seq}, UPPER: {maximum_seq}."
            assert len(overlapping_chunk) == UPPER - LOWER + 1, f"Invalid sequence length {len(sequence)}. UPPER: {UPPER}, LOWER: {LOWER}, Overlap: {overlapping_chunk}, LOWER: {minimum_seq}, UPPER: {maximum_seq}."
        except AssertionError:
            breakpoint()
       
        for idx, i in enumerate(range(LOWER, UPPER+1), 0):
            stranded_counts[sru][i][overlapping_chunk[idx]] += 1

    return stranded_counts


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--window_size", type=int, default=500)
    parser.add_argument("--inn", type=str, default="bed_out/enrichment")
    parser.add_argument("--out", type=str, default="bed_out/enrichment/out")

    files = []

    args = parser.parse_args()
    window_size = args.window_size
    inn = Path(args.inn)
    out = Path(args.out)
    out.mkdir(exist_ok=True, parents=True)
    # bucket_id = args.bucket_id
    
    files = [file for file in inn.glob("*.bed")]
    reduce_enrichment(files, window_size, out, 'gene')
