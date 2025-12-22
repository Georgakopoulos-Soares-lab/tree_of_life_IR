
GENOMES = "../nonbdna_pipeline/ncbi_dataset/data"
ACCESSIONS_LISTING = "accessions_listing.txt"
EXTRACTIONS="../nonbdna_pipeline/nonbdna_data_extractions"
PATTERN ?="IR"
EXTRACTIONS_DIR = "$(EXTRACTIONS)/extractions_$(PATTERN)"
GFF = "GFF_FILES.txt"
TREE_OF_LIFE = "../nonbdna_pipeline/data/assembly_summary_with_tree.csv.gz"
BUCKETS ?= 2
TEMPDIR ?= "tempdir"
DESIGN ?= "design_IR_with_tree.normal.csv"
SCHEDULE ?= "schedule_test.json"

# GFF Files
build_gff_db:
	@python build_gff_db.py ncbi_dataset/data gff_files

# Fetch Validated
fetch_validated:
	cat $(EXTRACTIONS_DIR)/log_debug_nonbdna/* | grep -a 'passed all checks' | awk -F ' ' '{ print $$7 }' > validated_accessions.txt

# Create Design
create_design:
	@python create_des.py $(EXTRACTIONS_DIR) --gff_files $(GFF) --pattern $(PATTERN) --tree_of_life $(TREE_OF_LIFE)
	# filter Design
	# TODO


# Create Schedule for GFF
# use GFF listing here
schedule:
	# @find $(EXTRACTIONS_DIR) -type f -name "*.tsv" -exec readlink -f {} \;  > $(ACCESSIONS_LISTING)
	@python scheduling.py --total_buckets $(BUCKETS) --accessions_path $(GFF)

run_density:
	@bash submit_density_IR.sh $(SCHEDULE) $(DESIGN) $(BUCKETS) $(TEMPDIR)
