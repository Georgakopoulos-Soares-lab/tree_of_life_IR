PYTEST?=micromamba run -n nonbdna pytest
TESTS_DIR=tests
REPORT_DIR=.reports
REPORT_FILE=$(REPORT_DIR)/pytest_report.txt
# 
GFF_INDIR=nonbdna_pipeline/gff_db  
MIN_PARTITION=0
MAX_PARTITION=9
PARTITION="spacer_length"
BUCKET=0
PATTERN=IR

.PHONY: test test-verbose ci-clean

WINDOW_SIZE=10

run:
	@python nonbdna_pipeline/tss_tes_processing.py nonbdna_pipeline/new_schedule_4_cf186650022b49b7b8c27d123dccc3ca.json --indir nonbdna_pipeline/extractions_IR --gff_indir nonbdna_pipeline/gff_db -p IR -s 0 --ignore_errors --window_size $(WINDOW_SIZE)

ci-clean:
	@mkdir -p $(REPORT_DIR)
	@rm -f $(REPORT_FILE)

# Default test run (quiet)
test: ci-clean
	$(PYTEST) -q $(TESTS_DIR) | tee $(REPORT_FILE)

test-verbose: ci-clean
	$(PYTEST) -vv -s --disable-warnings --durations=10 --maxfail=1 $(TESTS_DIR) | tee $(REPORT_FILE)

dag:
	snakemake --dag -s tss_tes_pipeline.smk \
		--configfile nonbdna_pipeline/config_IR.yaml | dot -Tpng -o tss_tes_dag.png

snake:
	bash submit_tss_tes.sh


motif_run:
	gff-motif-coverage nonbdna_pipeline/new_schedule_4_cf186650022b49b7b8c27d123dccc3ca.json \
			--indir nonbdna_pipeline/extractions_IR \
			--pattern $(PATTERN) \
			--gff_indir $(GFF_INDIR) \
			--bucket_id $(BUCKET) \
			--partition_col $(PARTITION) \
			--min_partition $(MIN_PARTITION) \
			--max_partition $(MAX_PARTITION) \
 			--use_biotype \
			--gff_suffix .gff

tss_tes_run:
	echo "Hi"

inspect_motif_run:
	gzcat nonbdna_pipeline/extractions_IR/gff_motif_coverage/gff_motif_coverage_IR_bucket_0.tsv.gz


# MUTATIONS
process_vcf:
	python mut_pipeline/process_VCF_step_1.py \
	--vcf_in mut_pipeline/test_data/sample.vcf.gz \
	--reference mut_pipeline/test_data/reference.fasta \
	--annotate
def main():
