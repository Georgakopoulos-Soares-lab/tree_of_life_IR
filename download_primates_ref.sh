#!/bin/bash

mkdir -p primates_ref
cd primates_ref
# -----------------------------------------------------------------------------------------------------
# Genomes
# >>

mkdir -p fasta 
cd fasta/

# Primates
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.pri.cur.20231205.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/mSymSyn1.pri.cur.20240514.fasta.gz
# Homo sapiens
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
# GRCh38
# hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Gunzip
gunzip * 

# -----------------------------------------------------------------------------------------------------
# Satellite Regions
# >>
cd ..
mkdir -p centromeres
cd centromeres/

# Homo Sapiens CenSat
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed

# Primates 
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/mGorGor1_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/mPonPyg2_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/repeats/mSymSyn1_v2.1_CenSat_v2.0.bb

# convert to BED format
bigBedToBed mGorGor1_v2.0_CenSat_v2.0.bb mGorGor1_v2.0_CenSat_v2.0.bed
bigBedToBed mPanPan1_v2.0_CenSat_v2.0.bb mPanPan1_v2.0_CenSat_v2.0.bed
bigBedToBed mPanTro3_v2.0_CenSat_v2.0.bb mPanTro3_v2.0_CenSat_v2.0.bed
bigBedToBed mPonAbe1_v2.0_CenSat_v2.0.bb mPonAbe1_v2.0_CenSat_v2.0.bed
bigBedToBed mPonPyg2_v2.0_CenSat_v2.0.bb mPonPyg2_v2.0_CenSat_v2.0.bed
bigBedToBed mSymSyn1_v2.1_CenSat_v2.0.bb mSymSyn1_v2.1_CenSat_v2.0.bed

# -----------------------------------------------------------------------------------------------------
echo "DONE!"
cd ../..
