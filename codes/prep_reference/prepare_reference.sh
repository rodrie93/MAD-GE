#!/bin/bash

# Preparation pipeline for reference transcriptome and mutant-aware maps
# Software requirements for this script:
    # 1. Python3
    # 2. Bedtools
    # 3. gffread
    # 4. Salmon

# Usage: bash annot_prep_edit.sh <genome.fa> <annot.gff3>

set -e # Exit immediately if a command exits with a non-zero status.

GENOME=$1
ANNOT=$2
PYTHON_UTILS="../codes/prep_reference/prep_reference.py" 
THREADS=8

if [ -z "$GENOME" ] || [ -z "$ANNOT" ]; then
    echo "Usage: bash annot_prep_edit.sh <genome.fa> <annot.gff3>"
    exit 1
fi

echo "1. Filtering for MANE Select transcripts..."
# Preserves GFF header for gffread compatibility
(grep "^#" "$ANNOT" || true; grep "MANE_Select" "$ANNOT") > mane_filtered.gff3

echo "2. Generating custom annotations and transcript maps..."
python3 $PYTHON_UTILS parse_gff mane_filtered.gff3 new_myfmt.gff tx_exon_map.tsv

echo "3. Extracting raw exon sequences from genome..."
# Converts the custom 10-col format to standard BED to guarantee bedtools compatibility
awk -F'\t' '{OFS="\t"; print $1, $4-1, $5, $3, ".", $7}' new_myfmt.gff > tmp_exons.bed
bedtools getfasta -fi "$GENOME" -bed tmp_exons.bed -nameOnly -fo tmp_exons.fa

echo "4. Generating final stranded fasta..."
python3 $PYTHON_UTILS format_fasta tmp_exons.fa new_myfmt.gff exons.final.fasta

echo "5. Building wild-type Reference Transcriptome..."
gffread -w ref_MANE_transcriptome.fa -g "$GENOME" mane_filtered.gff3

echo "6. Indexing wild-type Transcriptome with Salmon..."
salmon index -t ref_MANE_transcriptome.fa -i MANE-index -p $THREADS -k 31

echo "Cleaning up temporary files..."
rm mane_filtered.gff3 tmp_exons.bed tmp_exons.fa

echo "=========================================="
echo "Pipeline finished successfully! Outputs:"
echo " 1. new_myfmt.gff"
echo " 2. tx_exon_map.tsv"
echo " 3. exons.final.fasta"
echo " 4. ref_MANE_transcriptome.fa"
echo " 5. MANE-index/"
echo "=========================================="
