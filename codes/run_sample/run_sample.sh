#!/bin/bash

# Usage: Run from within the sample directory containing the WES-generated VCF files.
# bash run_sample.sh
# Software requirements:
	# 1. BCFtools
	# 2. BEDOPS
	# 3. BEDTOOLS


set -e # Exit immediately if a command fails

# ==========================================
# CONFIGURATION
# ==========================================
REF_DIR="../../ref/"             # Path containing prep_ref.sh outputs
PYTHON_UTILS="../../codes/run_sample/mut_generation.py"
SAMPLE=$(basename $(pwd))

echo "Running sample $SAMPLE"

echo "1. Grepping PASS somatic calls..."
bcftools view -f PASS *vcf.gz -Oz -o ${SAMPLE}_pass.vcf.gz
zcat ${SAMPLE}_pass.vcf.gz | vcf2bed > filtered.bed

echo "2. Hard filter artifact indels..."
python3 $PYTHON_UTILS filter_indels filtered.bed hardfiltered.bed

echo "3. Call exons which have a somatic mutation..."
bedtools intersect -a $REF_DIR/new_myfmt.gff -b hardfiltered.bed -wa -u > mut_exons_regions.gff

echo "4. Mutate the reference exons to have the somatic callings..."
python3 $PYTHON_UTILS mutate_exons mut_exons_regions.gff hardfiltered.bed $REF_DIR/exons.final.fasta mutated_exons.fa

echo "5. Revcom..."
python3 $PYTHON_UTILS revcom mutated_exons.fa mutated_exons_final.fa

echo "6. Build isoforms from mutant exons..."
python3 $PYTHON_UTILS build_isoforms $REF_DIR/tx_exon_map.tsv $REF_DIR/exons.final.fasta mutated_exons_final.fa mutated_isoforms.fa

echo "7. Copying reference transcriptome to working directory..."
cp $REF_DIR/ref_MANE_transcriptome.fa ./${SAMPLE}-expanded_MANE_transcriptome.fa

echo "8. Adding mutant isoforms to reference..."
cat mutated_isoforms.fa >> ./${SAMPLE}-expanded_MANE_transcriptome.fa

echo "Removing tmps..."
rm filtered.bed mut_exons_regions.gff mutated_exons.fa

echo "Done. The full expanded transcriptome can be found in the file '${SAMPLE}-expanded_MANE_transcriptome.fa'"
echo "FINISHED SCRIPT"
