from _process_VCF_step_2 import *

def load_VCF(VCF_in: str, fasta: str = None):
    intervals = []
    vcf = pysam.VariantFile(VCF_in)
    # gs = Fasta(fasta).extract_gs()
    for record in vcf:
        ref = record.ref
        for idx, alt in enumerate(record.alts, start=0):
            if record.info["TYPE"][idx] != "snp":
                continue
            transition = self.spectrum[f"{ref}>{alt}"]
            intervals.append(
                Interval(
                    chrom=record.chrom,
                    start=record.start,
                    end=record.stop,
                    name=transition,
                    score=str(round(record.info["AF"][idx], 4)),
                )
            )
    vcf.close()
    intervals_bed = BedTool(intervals)
    return intervals_bed


def retrieve_motifs(VCF_in: str, motifs: str):
    intervals_bed = load_VCF(VCF_in)
    motifs_df = pd.read_table(motifs)
    motif_bed = BedTool.from_dataframe(motifs_df).sort()

    # intersect_df = pl.read_csv(
    #    intervals_bed.intersect(motifs_bed).fn,
    #    has_header=False,
    #    separator="\t",
    #    new_columns=
    # )


def main():
    import argparse

    parser = argparse.ArgumentParser(description=""".""")
    parser.add_argument("--vcf", type=str)
    parser.add_argument("--motif", type=str)
    parser.add_argument("--fasta")
    args = parser.parse_args()
