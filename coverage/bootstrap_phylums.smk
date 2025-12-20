from termcolor import colored
import csv
from gff_utils import CoverageExtractor
from bootstrap_enrichment import Bootstrapper

out = Path(config['out']).resolve()
out.mkdir(exist_ok=True)
mode = config['mode']
alpha = round(float(config['alpha']), 2)
DESIGN = config['DESIGN']
partition_col = config['partition_col']
compartment = config['compartment']
polarity_type = config['polarity_type']

# load phylums
PHYLUMS = set()
delimiter = CoverageExtractor._sniff_delimiter(DESIGN)
with open(DESIGN, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter=delimiter)
    for row in reader:
        PHYLUMS.add(row['phylum'].replace(' ', '-'))
PHYLUMS = list(filter(len, PHYLUMS))
total_phylums = len(PHYLUMS)
info_color = "red" if total_phylums == 0 else "blue"
print(colored(f"Total phylums detected: {total_phylums}.", info_color))
# <<

# constants
BIOTYPES = config['BIOTYPES']
if not BIOTYPES:
    BIOTYPES = ["protein_coding", "non_coding", "."]
DOMAINS = ["Archaea", "Bacteria", "Eukaryota", "Viruses"]
# DOMAINS = ["Bacteria"]
SITES = ["TSS", "TES"]
# <<

# create directories
dest_dir_phylum = Path(f"{out}/mode_{mode}_partition_{partition_col}/phylum_bootstrap")
dest_dir_phylum.mkdir(exist_ok=True)

print(f"CHOSEN MODE: `{mode}`.")
print(f"Biotypes: `{BIOTYPES}`.")
print(f"Redirecting phylum level outputs to --> `{dest_dir_phylum}`.")
# <<

rule all:
    input:
        expand(['%s/mode_%s_partition_%s/phylum_bootstrap/enrichment_bootstrap_alpha_%s_GC.{biotype}.{site}.%s.phylum.{phylum}.csv' % (out, mode, partition_col, alpha, mode)], 
biotype=BIOTYPES,
											site=SITES, 
											phylum=PHYLUMS),

rule taxonomy_phylum_bootstrap:
    input:
        DESIGN,
        '%s/mode_%s_partition_%s/enrichment_%s_%s.{biotype}.{site}.parquet' % (out, mode, partition_col, compartment, mode),
    output:
        '%s/mode_%s_partition_%s/phylum_bootstrap/enrichment_bootstrap_alpha_%s.{biotype}.{site}.%s.phylum.{phylum}.csv' % (out, mode, partition_col, alpha, mode),
    params:
        window_size=int(config['window_size']),
        alpha=round(float(config['alpha']), 2),
        N=int(config['N']),
        partition_col=config["partition_col"],
	polarity=int(config["polarity"]),
        partition_list=config["partition_list"],
        mode=config['mode'],
    run:
        bootstrapper = Bootstrapper(design=input[0], enrichment_file=input[1], params=params)
        bootstrapper.bootstrap_enrichment(taxonomic_rank="phylum", 
                                          rank=wildcards.phylum.replace('-', ' '), 
                                          polarity=params.polarity,
                                          partition_col=params.partition_col,
                                          partition_list=params.partition_list,
                                          output=output[0])
