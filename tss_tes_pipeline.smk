# Snakemake pipeline to run tss_tes_processing.py for all buckets and merge outputs
import os
import json
from pathlib import Path
import tarfile
import pybedtools

# Prefer config file; can override with --configfile
configfile: "nonbdna_pipeline/config_IR.yaml"

SCHEDULE = str(Path(config["schedule"]).resolve())
INDIR = str(Path(config["indir"]).resolve())
GFF_DIR = str(Path(config.get("gff_dir", config.get("gff_indir", ""))).resolve())
GFF_TSS_TES_DIR = str(Path(config.get("gff_indir_tss_tes")).resolve())
ASM = str(Path(config.get("assembly_summary", "")).resolve()) if config.get("assembly_summary") else ""
PATTERN = config.get("pattern", "IR")
WINDOW_SIZE = int(config.get("window_size", 500))
STRAND_POLARITY = int(config.get("strand_polarity", 0))
USE_BIOTYPE = str(config.get("use_biotype", 0))
MIN_PARTITION = config.get("min_partition")
MAX_PARTITION = config.get("max_partition")
PARTITION = config.get("partition_col")
GFF_SUFFIX = config.get("gff_suffix", ".agat.gff")
PSEUDOGENES_TO_GENES = str(config.get("pseudogenes_to_genes", 1))

assert Path(SCHEDULE).is_file(), f"Schedule file {SCHEDULE} not found."
assert Path(INDIR).exists(), f"Input dir {INDIR} not found."
assert GFF_DIR and Path(GFF_DIR).exists(), f"GFF dir {GFF_DIR} not found. Set gff_dir in config."
assert GFF_TSS_TES_DIR and Path(GFF_TSS_TES_DIR).exists(), f"GFF dir {GFF_TSS_TES_DIR} not found. Set gff_tss_tes_dir in config."

print(f"USING BIOTYPE: {USE_BIOTYPE}")

OUTDIR = Path(INDIR).resolve()
DENSITY_DIR = OUTDIR.joinpath("tss_tes_density")
COVERAGE_DIR = OUTDIR.joinpath("gff_motif_coverage")

Path("garbage").mkdir(exist_ok=True)
pybedtools.helpers.set_tempdir("garbage")

def _load_bucket_ids(schedule_path: str) -> list[int]:
    with open(schedule_path, "r", encoding="UTF-8") as f:
        buckets = json.load(f)
    try:
        return sorted(int(k) for k in buckets.keys())
    except Exception:
        return list(range(len(buckets)))


BUCKETS = _load_bucket_ids(SCHEDULE)


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
        f"{DENSITY_DIR}/tss_tes_density_{PATTERN}_bucket_{{bid}}.tsv.gz"
    params:
        schedule=SCHEDULE,
        indir=INDIR,
        gff=GFF_TSS_TES_DIR,
        pattern=PATTERN,
        asm=ASM,
        window=WINDOW_SIZE,
        sp=STRAND_POLARITY,
        use_biotype_flag=("--use_biotype" if USE_BIOTYPE in ("1", "true", "True", "yes") else "")
    wildcard_constraints:
        bid=r"\d+"
    shell:
        r"""
        tss-tes-processor {params.schedule} \
            -i {params.indir} \
            --gff_indir {params.gff} \
            -p {params.pattern} \
            --window_size {params.window} \
            --strand_polarity {params.sp} \
            --bucket_id {wildcards.bid} \
            --assembly_summary {params.asm} \
            {params.use_biotype_flag}
        """


rule motif_coverage_bucket:
    output:
        f"{COVERAGE_DIR}/gff_motif_coverage_{PATTERN}_bucket_{{bid}}.tsv.gz"
    params:
        schedule=SCHEDULE,
        indir=INDIR,
        gff=GFF_DIR,
        pattern=PATTERN,
        gff_suffix=GFF_SUFFIX,
        partition_col=PARTITION,
        min_partition=MIN_PARTITION,
        max_partition=MAX_PARTITION,
        use_biotype_flag=("--use_biotype" if USE_BIOTYPE in ("1", "true", "True", "yes") else ""),
        pseudo_flag=("--pseudogenes_to_genes" if PSEUDOGENES_TO_GENES in ("1", "true", "True", "yes") else "")
    wildcard_constraints:
        bid=r"\d+"
    shell:
        """
        gff-motif-coverage {params.schedule} \
            --bucket_id {wildcards.bid} \
            -p {params.pattern} \
            -i {params.indir} \
            -g {params.gff} \
            --partition_col {params.partition_col} \
            --min_partition {params.min_partition} \
	    --max_partition {params.max_partition} \
            --gff_suffix {params.gff_suffix}
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
