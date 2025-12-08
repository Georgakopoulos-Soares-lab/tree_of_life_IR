import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator

import attr
import pandas as pd
from attr import field
from utils import parse_fasta


@attr.s(slots=True, auto_attribs=True)
class STRDetector:

    min_sru: int = field(
        default=2, metadata={"help": "Minimum SRU length to detect"}, converter=int
    )
    max_sru: int = field(
        default=9, metadata={"help": "Maximum SRU length to detect"}, converter=int
    )
    outdir: Path = field(
        default=Path("."),
        metadata={"help": "Output directory for results"},
        converter=Path,
    )
    motif: str = field(init=False)
    motifs: dict[int, str] = field(init=False)

    def __attrs_post_init__(self) -> None:
        self.outdir.mkdir(exist_ok=True)
        print(f"Output directory --> {self.outdir}")
        self.max_sru = max(self.min_sru, self.max_sru)
        self.motif: str = "[agct]"
        self.motifs: dict[int, str] = {
            1: r"(%s)\1{10,}" % self.motif,
            2: r"(%s)\1{4,}" % (self.motif * 2),
            3: r"(%s)\1{3,}" % (self.motif * 3),
            4: r"(%s)\1{2,}" % (self.motif * 4),
            5: r"(%s)\1{2,}" % (self.motif * 5),
            6: r"(%s)\1{2,}" % (self.motif * 6),
            7: r"(%s)\1{2,}" % (self.motif * 7),
            8: r"(%s)\1{2,}" % (self.motif * 8),
            9: r"(%s)\1{2,}" % (self.motif * 9),
        }

    def _detect(self, sequence: str) -> Iterator[dict]:
        data = defaultdict(list)
        for SRU in range(self.min_sru, self.max_sru + 1):
            motif = self.motifs.get(SRU, None)
            if motif is None:
                raise ValueError(
                    f"SRU {SRU} is not supported. Supported SRUs are from {self.min_sru} to {self.max_sru}."
                )
            matches = re.finditer(motif, sequence, re.IGNORECASE)
            for match in matches:
                start = match.start()
                end = match.end()
                microsatellite = match.group(0)
                sequence_of_arm = match.group(1)
                sru = len(sequence_of_arm)
                assert sru == SRU, f"SRU {sru} does not match expected {SRU}"
                length = end - start
                consensus_repeats = length // sru
                data["start"].append(start)
                data["end"].append(end)
                data["sequence_of_arm"].append(sequence_of_arm)
                data["arm_length"].append(len(sequence_of_arm))
                data["sequence_of_spacer"].append(".")
                data["spacer_length"].append(0)
                data["sequence"].append(microsatellite)
                data["sequence_length"].append(length)
                data["sru"].append(sru)
                data["consensus_repeats"].append(consensus_repeats)
                data["type"].append("STR")
                data["method"].append("Consensus Motif")
            yield data

    @staticmethod
    def extract_name(accession: str) -> str:
        accession = Path(str(accession).replace(".gz", ""))
        if accession.name.endswith(".fasta"):
            accession = accession.name.split(".fasta")[0]
            return accession
        elif accession.name.endswith(".fa"):
            accession = accession.name.split(".fa")[0]
            return accession
        elif accession.name.endswith(".fna"):
            accession = accession.name.split(".fna")[0]
            return accession
        raise ValueError(f"Accession {accession} is not a valid fasta file.")

    def detect(self, accession: str, save_output: bool = True) -> pd.DataFrame:
        df = []
        for seqID, seq in parse_fasta(accession):
            result = []
            for _result in self._detect(seq):
                if _result:
                    _result = pd.DataFrame(_result)
                    result.append(_result)
            if not result:
                continue
            data = (
                pd.concat(result, ignore_index=True)
                .drop_duplicates(subset=["start", "end"])
                .sort_values(by=["start"], ascending=True)
                .reset_index(drop=True)
            )
            data.loc[:, "seqID"] = seqID
            data = data[
                [
                    "seqID",
                    "start",
                    "end",
                    "sequence_of_arm",
                    "arm_length",
                    "sequence_of_spacer",
                    "spacer_length",
                    "sequence",
                    "sequence_length",
                    "sru",
                    "consensus_repeats",
                    "type",
                    "method",
                ]
            ]
            df.append(data)
        if not df:
            print(f"No STRs found in {accession}")
            return df
        df = pd.concat(df, ignore_index=True)
        outfile = (
            self.outdir / f"{STRDetector.extract_name(accession)}_STR_consensus.tsv"
        )
        if save_output:
            df.to_csv(outfile, sep="\t", index=False, header=True, mode="w")
        return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Detect STRs in a sequence.")
    parser.add_argument("accession", type=str, help="The sequence to detect STRs in")
    parser.add_argument("--outdir", type=str, default=".")
    parser.add_argument(
        "--min_sru", type=int, default=2, help="Minimum SRU length to detect"
    )
    parser.add_argument(
        "--max_sru", type=int, default=9, help="Maximum SRU length to detect"
    )
    args = parser.parse_args()
    detector = STRDetector(
        min_sru=args.min_sru, max_sru=args.max_sru, outdir=args.outdir
    )
    STR_df = detector.detect(args.accession, save_output=True)
