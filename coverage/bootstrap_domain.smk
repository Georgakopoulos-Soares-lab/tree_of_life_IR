from termcolor import colored
import csv
from tqdm import tqdm
from gff_utils import CoverageExtractor
from bootstrap_enrichment import Bootstrapper
import shutil
import pandas as pd
import tempfile

out = Path(config['out']).resolve()
out.mkdir(exist_ok=True)
mode = config['mode']
alpha = round(float(config['alpha']), 2)
design = config['DESIGN']
partition_col = config['partition_col']
compartment = config['compartment']
polarity_type = config['polarity_type']
total_buckets = int(config['total_buckets'])
window_size = int(config['window_size'])

# constants
BIOTYPES = config['BIOTYPES']
if not BIOTYPES:
    BIOTYPES = ["protein_coding", "non_coding", "."]
DOMAINS = ["Archaea", "Bacteria", "Eukaryota", "Viruses"]
# DOMAINS = ["Bacteria"]
SITES = ["TSS", "TES"]
# <<

valid_bucket_ids = []
for bucket in range(total_buckets):
    for site in SITES:
        for biotype in BIOTYPES:
            if Path(f'{out}/mode_{mode}_partition_{partition_col}/enrichment_bucket_{bucket}_{window_size}_{compartment}_{mode}_GC.{biotype}.{site}.txt').is_file():
                valid_bucket_ids.append(bucket)

print(f"Total valid bucket ids: {len(valid_bucket_ids)}.")
# load phylums
PHYLUMS = set()
delimiter = CoverageExtractor._sniff_delimiter(design)
with open(design, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter=delimiter)
    for row in reader:
        PHYLUMS.add(row['phylum'].replace(' ', '-'))
PHYLUMS = list(filter(len, PHYLUMS))
total_phylums = len(PHYLUMS)
info_color = "red" if total_phylums == 0 else "blue"
print(colored(f"Total phylums detected: {total_phylums}.", info_color))
# <<


# create directories
dest_dir_domain = Path(f"{out}/mode_{mode}_partition_{partition_col}/domain_bootstrap")
dest_dir_domain.mkdir(exist_ok=True)
dest_dir_phylum = Path(f"{out}/mode_{mode}_partition_{partition_col}/phylum_bootstrap")
dest_dir_phylum.mkdir(exist_ok=True)

print(f"CHOSEN MODE: `{mode}`.")
print(f"Biotypes: `{BIOTYPES}`.")
print(f"Redirecting domain level outputs to --> `{dest_dir_domain}`.")
print(f"Redirecting phylum level outputs to --> `{dest_dir_phylum}`.")
# <<

rule all:
    input:
        expand(['%s/mode_%s_partition_%s/domain_bootstrap/enrichment_bootstrap_alpha_%s.{biotype}.{site}.%s.domain.{domain}.csv' % (out, mode, partition_col, alpha, mode),
                '%s/mode_%s_partition_%s/domain_bootstrap/enrichment_phylums.{biotype}.{site}.%s.{domain}.csv' % (out, mode, partition_col, mode)],
                                 biotype=BIOTYPES,
                                 site=SITES,
                                 domain=DOMAINS)
rule concat:
    input:
        expand('%s/mode_%s_partition_%s/enrichment_bucket_{bucket}_%s_%s_%s_GC.{{biotype}}.{{site}}.txt' % (out, mode, partition_col, window_size, compartment, mode),
                        bucket=valid_bucket_ids
                ),
    output:
        "%s/mode_%s_partition_%s/enrichment_%s_%s.{biotype}.{site}.txt" % (out, mode, partition_col, compartment, mode),
        "%s/mode_%s_partition_%s/enrichment_%s_%s.{biotype}.{site}.parquet" % (out, mode, partition_col, compartment, mode)
    run:
        df_all = []
        with tempfile.NamedTemporaryFile(delete=False, mode="a") as tmpfile:
            for i, file in tqdm(enumerate(input), total=len(input)):
                df = pd.read_table(file)
                df.to_csv(tmpfile, sep="\t", index=False, header=i==0, mode="a+")
            infile = tmpfile.name
            # df_all.append(df)
        # df_all = pd.concat(df_all, ignore_index=True)
        # df_all.to_csv(output[0], sep="\t", index=False, header=True, mode="w")
        df_all = pd.read_table(infile)
        shutil.move(infile, output[0])
        # df_all.to_csv(output[0], sep="\t", index=False, header=True, mode="w")
        df_all.to_parquet(output[1], compression="snappy")

rule taxonomy_domain_bootstrap:
    input:
        '%s/mode_%s_partition_%s/enrichment_%s_%s.{biotype}.{site}.txt' % (out, mode, partition_col, compartment, mode),
        '%s/mode_%s_partition_%s/enrichment_%s_%s.{biotype}.{site}.parquet' % (out, mode, partition_col, compartment, mode),
        design,
    output:
        '%s/mode_%s_partition_%s/domain_bootstrap/enrichment_bootstrap_alpha_%s.{biotype}.{site}.%s.domain.{domain}.csv' % (out, mode, partition_col, alpha, mode),
        '%s/mode_%s_partition_%s/domain_bootstrap/enrichment_phylums.{biotype}.{site}.%s.{domain}.csv' % (out, mode, partition_col, mode)
    params:
        window_size=int(config['window_size']),
        alpha=round(float(config['alpha']), 2),
        polarity=int(config["polarity"]),
        mode=config['mode'],
        N=int(config['N']),
        partition_col=config['partition_col'],
        partition_list=config['partition_list'],
        join_templates=config['join_templates']
    run:
        bootstrapper = Bootstrapper(design=input[2], enrichment_file=input[1], params=params)
        bootstrapper.bootstrap_enrichment(taxonomic_rank="domain",
                                          rank=wildcards.domain,
                                          partition_col=params.partition_col,
                                          polarity=params.polarity,
					  partition_list=params.partition_list,
                                          output=output[0])
        bootstrapper.average_phylums(domain=wildcards.domain,
                                     output=output[1],
                                     join_templates=params.join_templates)
