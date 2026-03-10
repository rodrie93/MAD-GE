#!/bin/bash

# Usage: Run from within the sample directory containing both the RNA FASTQ files 
# and the newly generated expanded transcriptome.
# bash salmon_mut_quant.sh

# Required Softwares:
# 1. Salmon

set -e # Exit immediately if a command fails

start_time=$(date +%s.%N)

# ==========================================
# CONFIGURATION
# ==========================================
REF_DIR="../../ref/"             # Path containing prep_ref.sh outputs
SALMON_BIN="salmon"                        # Command or path to salmon
THREADS=8
SAMPLE=$(basename $(pwd))

# Automatically grab FASTQ files (assuming they are already trimmed/QC'd)
R1=$(ls *R1*fastq.gz | head -n 1) 
R2=$(ls *R2*fastq.gz | head -n 1)

echo "=========================================="
echo "Starting Salmon Quantification for Sample: $SAMPLE"
echo "Start time: $(date)"
echo "=========================================="

# 1. Traditional Quantifications (Wild-type only)
echo "--> Running traditional quantification against wild-type reference..."
$SALMON_BIN quant -i $REF_DIR/MANE-index -l A -1 $R1 -2 $R2 -p $THREADS \
    --validateMappings --dumpEq --writeUnmappedNames -o MANE_trad/

# 2. Index the new expanded transcriptome
# This links directly to the output naming convention of our mut_pipe.sh script
EXPANDED_TX="${SAMPLE}-expanded_MANE_transcriptome.fa"

echo "--> Indexing expanded patient-specific transcriptome ($EXPANDED_TX)..."
$SALMON_BIN index -t $EXPANDED_TX -i ${SAMPLE}-MANE-index -p $THREADS -k 31

# 3. Mutant Quantifications (Expanded)
echo "--> Running mutant-aware quantification..."
$SALMON_BIN quant -i ${SAMPLE}-MANE-index -l A -1 $R1 -2 $R2 -p $THREADS \
    --validateMappings --dumpEq --writeUnmappedNames -o quants_muts/

end_time=$(date +%s.%N)
elapsed_time=$(echo "($end_time - $start_time) / 60" | bc)

echo "=========================================="
echo "Finished Sample: $SAMPLE"
echo "Elapsed time: $elapsed_time minutes"
echo "Outputs are ready in MANE_trad/ and quants_muts/"
echo "=========================================="
