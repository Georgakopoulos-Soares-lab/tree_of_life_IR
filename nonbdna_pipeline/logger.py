import time 
import logging
import attr
from attr import field
from pathlib import Path

@attr.s
class Logger:

    LOG_INTERVAL: int = 240
    files_processed: int = field(init=False, default=0)
    total_files: int = field(init=True, converter=int)

    def _setup_logging(self, bucket_id: int) -> None:
        DATE = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            filename=self.log_dir.joinpath(f"tss_tes_processing_{DATE}_{bucket_id}.log"),
            # handlers=[
            #     logging.FileHandler(self.log_dir.joinpath(f"tss_tes_processing_{DATE}.log")),
            #     logging.StreamHandler()
            # ]
        )
    def _log_progress(self) -> None:
        while True:
            prc = self.files_processed / self.total_files * 1e2 if self.total_files > 0 else 0.0
            logging.info(f"Progress: {prc:.2f}% files processed.")
            time.sleep(Logger.LOG_INTERVAL)