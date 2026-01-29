# Snakemake pipeline to run tss_tes_processing.py for all buckets and merge outputs
import os
import json
from pathlib import Path
from nonbdna_pipeline.utils import load_bucket, load_bucket_ids

# Prefer config file; can override with --configfile
configfile: "config_IR.yaml"

SCHEDULE = str(Path(config["schedule"]).resolve())
INDIR = str(Path(config["indir"]).resolve())
GFF_DIR = str(Path(config.get("gff_dir", config.get("gff_indir", ""))).resolve())
ASM = str(Path(config.get("assembly_summary", "")).resolve()) if config.get("assembly_summary") else ""
PATTERN = config.get("pattern", "IR")
WINDOW_SIZE = int(config.get("window_size", 500))
STRAND_POLARITY = int(config.get("strand_polarity", 0))
USE_BIOTYPE = str(config.get("use_biotype", 0))
PSEUDOGENES_TO_GENES = str(config.get("pseudogenes_to_genes", 1))

assert Path(SCHEDULE).is_file(), f"Schedule file {SCHEDULE} not found."
assert Path(INDIR).exists(), f"Input dir {INDIR} not found."
assert GFF_DIR and Path(GFF_DIR).exists(), f"GFF dir {GFF_DIR} not found. Set gff_dir in config."

OUTDIR = Path(INDIR).resolve()
DENSITY_DIR = OUTDIR.joinpath("tss_tes_density")
COVERAGE_DIR = OUTDIR.joinpath("gff_motif_coverage")

def bucket_out(bid: int) -> str:
    return str(DENSITY_DIR.joinpath(f"tss_tes_density_{PATTERN}_bucket_{bid}.tsv.gz"))

PER_BUCKET_OUT = [bucket_out(bid) for bid in BUCKETS]
MERGED_OUT = str(DENSITY_DIR.joinpath(f"tss_tes_density_{PATTERN}_all_buckets.tsv.gz"))

def coverage_out(bid: int) -> str:
    return str(COVERAGE_DIR.joinpath(f"gff_motif_coverage_{PATTERN}_bucket_{bid}.tsv.gz"))

PER_BUCKET_COV = [coverage_out(bid) for bid in BUCKETS]
MERGED_COV_OUT = str(COVERAGE_DIR.joinpath(f"gff_motif_coverage_{PATTERN}_all_buckets.tsv.gz"))

rule all:
    input:
        MERGED_OUT,
        MERGED_COV_OUT

rule tss_tes_bucket:
    output:
        lambda wildcards: bucket_out(int(wildcards.bid))
    params:
        schedule=SCHEDULE,
        indir=INDIR,
        gff=GFF_DIR,
        pattern=PATTERN,
        asm=ASM,
        window=WINDOW_SIZE,
        sp=STRAND_POLARITY
    wildcard_constraints:
        bid="\d+"
    shell:
        r"""
        python -u nonbdna_pipeline/tss_tes_processing.py {params.schedule} \
            -i {params.indir} \
            --gff_indir {params.gff} \
            -p {params.pattern} \
            --window_size {params.window} \
            --strand_polarity {params.sp} \
            --bucket_id {wildcards.bid} \
            --assembly_summary {params.asm}
        """


rule motif_coverage_bucket:
    output:
        lambda wildcards: coverage_out(int(wildcards.bid))
    params:
        schedule=SCHEDULE,
        indir=INDIR,
        gff=GFF_DIR,
        pattern=PATTERN,
        use_biotype_flag=("--use_biotype" if USE_BIOTYPE in ("1", "true", "True", "yes") else ""),
        pseudo_flag=("--pseudogenes_to_genes" if PSEUDOGENES_TO_GENES in ("1", "true", "True", "yes") else "")
    wildcard_constraints:
        bid="\d+"
    shell:
        r"""
        python -u nonbdna_pipeline/gff_motif_coverage.py {params.schedule} \
            --bucket_id {wildcards.bid} \
            -p {params.pattern} \
            -i {params.indir} \
            -g {params.gff} \
            {params.use_biotype_flag} \
            {params.pseudo_flag}
        """


def _reduce_outputs(input, output) -> None:
    import gzip
    wrote_header = False
    with gzip.open(output[0], "wt") as fout:
        for infile in input:
            with gzip.open(infile, "rt") as fin:
                for i, line in enumerate(fin):
                    if i == 0:
                        if not wrote_header:
                            fout.write(line)
                            wrote_header = True
                    else:
                        fout.write(line)

rule reduce_tss_tes_buckets:
    input:
        PER_BUCKET_OUT
    output:
        MERGED_OUT
    run:
        _reduce_outputs(input, output)


rule reduce_gff_motif_coverage_buckets:
    input:
        PER_BUCKET_COV
    output:
        MERGED_COV_OUT
    run:
        _reduce_outputs(input, output)