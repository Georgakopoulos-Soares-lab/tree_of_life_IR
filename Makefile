PYTEST?=micromamba run -n nonbdna pytest
TESTS_DIR=tests
REPORT_DIR=.reports
REPORT_FILE=$(REPORT_DIR)/pytest_report.txt

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

# Verbose with timings and warnings disabled
# Includes -vv, -s (no capture), durations, and detailed summary
# Add --maxfail=1 to stop on first failure for quick iteration
# You can set PYTEST to plain 'pytest' if not using micromamba
# Example: make PYTEST=pytest test-verbose

test-verbose: ci-clean
	$(PYTEST) -vv -s --disable-warnings --durations=10 --maxfail=1 $(TESTS_DIR) | tee $(REPORT_FILE)

