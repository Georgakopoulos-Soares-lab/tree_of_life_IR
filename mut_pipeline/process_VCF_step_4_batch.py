import pickle
from process_VCF_step_3 import *

class PolymorphismDetector(BatchProcessor):
    def extract_baseline(self):
        files = self.fetch_files()
        models = dict()
        outdir = Path("models")
        outdir.mkdir(exist_ok=True)
        print(colored(f"Total files detected {len(files)}.", "green"))
        for species_taxid, VCF_in in tqdm(files.items()):
            model = TrinucleotideModel()
            reference = self.fasta_files[species_taxid]
            model.fit(VCF_in, reference)
            models[species_taxid] = model
            outfile = outdir.joinpath(f"trinucleotide_model_event_rate_{model.offset}_{species_taxid}.pkl")
            # model.saveas(outdir=self.outdir, species_taxid=species_taxid)
            with open(outfile, "wb") as f:
                pickle.dump(model, f)
        return models

def main():
    import argparse
    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("--location", type=str, default="processed_VCF")
    parser.add_argument(
        "--reference", type=str, default="reference/reference_files.txt"
    )
    parser.add_argument("--outdir", type=str, default="trinucleotide_model_out")
    args = parser.parse_args()
    detector = PolymorphismDetector(
        location=args.location, reference=args.reference, outdir=args.outdir
    ).extract_baseline()


if __name__ == "__main__":
    main()