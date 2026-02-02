# Tree of Life Inverted Repeats

This repository contains tools and pipelines for detecting non-B DNA motifs and computing coverage/enrichment (IR, MR, STR, DR, etc.). The code is packaged under the `nonbdna_pipeline` Python package and includes a small extraction utility that wraps the internal `MindiTool` extraction logic.

This README shows a minimal workflow to install the package locally with pip and run the extraction/merge step using the included CLI entrypoint (or module). It also shows how to run the bootstrap analysis script that computes taxonomic bootstrap confidence intervals.

## Prerequisites
 - Python 3.10+ (the project was developed on modern Python 3.11/3.12 builds)
 - git
 - system build tools for any compiled Python dependencies (if needed)
 - recommended: create an isolated virtual environment (venv or conda/micromamba)

## Install

```bash
# from the repository root
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .
```

## Requirements (optional)

If you prefer installing a requirements file, create a `requirements.txt` with the dependencies below (example). Then install with:

```bash
pip install -r requirements.txt
```

Example minimal `requirements.txt` (adjust versions as needed):

```
pandas
polars
pyranges
pybedtools
termcolor
attrs
pytest
pyyaml
python-dotenv
```

## Run nonbdna extraction (Stream & Merge)

The main extraction/merge utility is implemented in `nonbdna_pipeline/stream_and_merge_bucket.py` and provides a command-line `main()` function.

Example usage:

```bash
# run bucket 0 for STR pattern (uses defaults for schedule/indir unless you override)
python -m nonbdna_pipeline.stream_and_merge_bucket 0 \
    --pattern STR \
    --schedule schedule_tandem_extractions.json \
    --indir /path/to/extractions_STR_MR \
    --partition_col sru --min_partition 1 --max_partition 10
```

The script will read the schedule JSON file (default `schedule_tandem_extractions.json`), consult the log files under the configured `indir` to find validated processed extraction files, and write merged outputs under `indir/merged/<pattern>/`.

### Extract a single accession using MindiTool

If you want to extract non-B DNA from a single FASTA accession (for quick tests), use the internal `MindiTool` wrapper. Example Python usage:

```python
from nonbdna_pipeline.minditool import MindiTool

# Either set env var `nonBDNA` to the path of the non-B-DNA binary
# or pass it directly to the constructor
mt = MindiTool(nonBDNA="/path/to/nonbdna_executable", tempdir="./tmp_mindi")

# Run extraction for one accession (FASTA/FA/.gz). Pattern can be 'IR', 'MR', 'STR', etc.
mt.extract("example/GC0000001.1_reference.fna", pattern=["IR"])  # returns the MindiTool instance

# Get processed dataframe for that mode (if available)
df_ir = mt.to_dataframe("IR")
print(df_ir.head())

# Clean up temporary processed files when done
mt.cleanup()
```

This runs the underlying non-B-DNA binary on a single accession and writes processed `.processed.tsv` files into `tempdir` (or the current working directory if not set). Use this for quick checks before running the full per-bucket pipeline.

## Run taxonomic bootstrap (tss_tes_bootstrap)

The bootstrap helper script is `nonbdna_pipeline/tss_tes_boostrap.py`. It expects a density TSV (for example the merged tss/t es density output) and will compute bootstrap confidence intervals by taxonomic grouping.

Example to run for a particular taxonomy/biotype:

```bash
python -m nonbdna_pipeline.tss_tes_boostrap \
    /path/to/tss_tes_density_IR_all_buckets.tsv.gz \
    --rank domain --taxonomy Bacteria \
    --biotype protein_coding \
    --bootstrap_taxonomic_level family \
    --n_samples 1000 --window_size 500
```

## Outputs
- Extraction/merge outputs: written under `indir/merged/<pattern>/` (see `StreamAndMerge.process_bucket` for exact filenames).
- Bootstrap outputs: the script writes results under `$(density_file.parent)/bootstrap_results/` and the filename is derived from the input stem plus the rank/taxonomy.

## Further help
- See `nonbdna_pipeline/` for the available modules. The main scripts with CLI entrypoints are:
  - `stream_and_merge_bucket.py` (extraction/merge runner)
  - `gff_motif_coverage.py` (motif coverage processor)
  - `tss_tes_boostrap.py` (bootstrap analysis)

If you want, I can add example Snakemake rules that call these scripts for common datasets or create a short tutorial notebook demonstrating a full run on a small sample dataset.