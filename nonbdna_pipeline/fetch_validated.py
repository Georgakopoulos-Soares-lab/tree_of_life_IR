class PipelineIncomplete(Exception):
    pass

class Validated:

    def __init__(self, indir: str, pattern: str):
        from pathlib import Path
        from subprocess import run
        from termcolor import colored
        self.indir = Path(indir).resolve()
        if not self.indir.is_dir():
            raise PipelineIncomplete()
        self.debug_dir = self.indir.joinpath("log_debug_dna")
        if not self.debug_dir.is_dir():
            raise PipelineIncomplete()
        self.pattern = pattern
        self.product = Path(f"processed_succesfully_accessions_{pattern}.txt")

    def fetch_validated(self):
        run(f"cd {self.debug_dir} && cat * | grep -a 'passed all checks' | awk -F ' ' '{{ print $7 }}' > {self.product}", shell=True, check=True)

    def read(self) -> set[str]:
        validated = set()
        extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
        if not self.product.is_file():
            raise ValueError(f"Missing validated file `{self.product}`.")
        with open(self.product, mode="r", encoding="UTF-8") as fin:
            for line in fin:
                line = line.strip()
                accession_id = extract_id(line)
                validated.add(accession_id)
        total_validated = len(validated)
        color = "red" if not total_validated else "green"
        print(colored(f"Total validated assemblies: {total_validated}.", color))
        return validated
