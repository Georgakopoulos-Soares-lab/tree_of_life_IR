from pathlib import Path
import logging 
import json 
import gzip
import os
from attr import field
import attr
import pandas as pd
from nonbdna_pipeline.minditool import MindiTool
from nonbdna_pipeline.logger import Logger

@attr.s
class Extractor:
    outdir: Path = field(init=True, factory=lambda: Path(".").resolve())
    subdir: Path = field(init=False)
    logdir: Path = field(init=False)
    status: Path = field(init=False)
    def __attrs_post_init__(self) -> None:
        self.outdir = Path(self.outdir)
        self.outdir.mkdir(exist_ok=True)
        self.extractions_dir = self.outdir.joinpath("extractions")
        self.extractions_dir.mkdir(exist_ok=True)

        self.status = self.outdir.joinpath("status")
        self.status.mkdir(exist_ok=True)

        self.logdir = self.outdir.joinpath("log_debug_nonbdna")
        self.logdir.mkdir(exist_ok=True)

        self.subdir = self.outdir.joinpath("merged_bulk")
        self.subdir.mkdir(exist_ok=True)
        return

    @staticmethod
    def load_bucket(schedule: str, bucket_id: int) -> list[str]:
        """
        Load a bucket of accessions from a JSON file.
        """
        with open(schedule, mode="r", encoding="UTF-8") as f:
            buckets = json.load(f)
        return buckets.get(str(bucket_id), [])

    def _setup_logging(self, bucket_id: int) -> None:
        logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s', 
                        filename=f'{self.logdir}/mindi_tool_{bucket_id}.log', 
                        filemode='w'
                        )
        return

    def process_bucket(self, schedule: str, bucket_id: int, pattern: str) -> None:
        self._setup_logging(bucket_id=bucket_id)
        infiles = Extractor.load_bucket(schedule=schedule, bucket_id=bucket_id)
        total_accessions = len(infiles)
        if total_accessions == 0:
            logging.info(f"No accessions found for bucket `{bucket_id}` in schedule `{schedule}`.")
            return
        else:
            logging.info(f"Loaded `{total_accessions}` accessions for bucket `{bucket_id}` in schedule `{schedule}`.")
        fout = {}
        # Initializing processing
        for m in pattern.split(","):
            extractions_outfile = self.subdir.joinpath(f"bucket_{bucket_id}_mode_{m}.tsv.gz")
            empty_outfile = self.subdir.joinpath(f"bucket_{bucket_id}_mode_{m}_empty.tsv")
            fout[m] = gzip.open(extractions_outfile, "wt"), open(empty_outfile, "w", encoding="UTF-8")
        # if isinstance(pattern, str):
        #    pattern = [pattern]
        # Processing begins
        for i, accession in enumerate(infiles, 1):
            # logging.info(f"Processing accession: {accession}. Progress: {(i-1) * 1e2 / total_accessions:.2f}% (Bucket ID: {args.bucket_id}).")
            hunter = MindiTool(tempdir=self.outdir)
            result = hunter.extract(accession, 
                                    pattern=pattern)
            if result is None:
                logging.error(f"Extraction failed for accession: {accession}.")
                continue
            for m in pattern.split(","):
                try:
                    hunter.sanitize(accession, mode=m)
                except AssertionError:
                    logging.error(f"Extraction failed for accession: {accession} (mode {m}).")
                    os.remove(hunter.fn[m])
                    continue
                try:
                    result_df = pd.read_table(hunter.fnp[m])
                    if result_df.shape[0] == 0:
                        fout[m][1].write(accession + "\n")
                    # result_df.loc[:, "#assembly_accession"] = MindiTool.extract_id(accession)
                    result_df.to_csv(fout[m][0], 
                                     header=i==1, 
                                     index=False, 
                                     sep="\t")
                except pd.errors.EmptyDataError:
                    fout[m][1].write(accession + "\n")
                os.remove(hunter.fn[m])
        for m in pattern.split(","):
            for i in range(2):
                fout[m][i].close()
        logging.info(f"Completed processing for bucket `{bucket_id}` in schedule `{schedule}`.")
        # Bucket Complete Message
        status_outfile = self.status.joinpath(f"bucket_{bucket_id}_status.completed")
        with status_outfile.open("w", encoding="UTF-8") as fout:
            fout.write(f"Bucket {bucket_id} has been completed\n")

        return infiles

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="MindiTool Extraction and Validation")
    parser.add_argument("--bucket_id", type=int, default=0, help="Bucket ID for parallel processing")
    parser.add_argument("--schedule", type=str, default="all", help="Schedule for processing (e.g., all, daily, weekly)")
    parser.add_argument("--pattern", type=str, default="STR")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory for extracted files")
    args = parser.parse_args()
    if args.outdir == ".":
        pattern__repr__ = "_".join(args.pattern.split(","))
        args.outdir = f"extractions_{pattern__repr__}"

    Extractor(outdir=args.outdir)\
            .process_bucket(schedule=args.schedule, 
                            bucket_id=args.bucket_id,
                            pattern=args.pattern
                        )
