# MAD-GE: Mutant Aware Deconvolution of Gene Expression

MAD-GE builds personalized transcriptomes to quantify mutant alleles. Originally coded pre-LLM by a biologist without a strong CS background, this pipeline was recently refactored and tidied with Gemini to ensure transparency and reproducibility for publication.

## Quick Start & Installation

MAD-GE relies on standard bioinformatics tools:
- Python3
- BEDTOOLS
- BCFtools
- BEDOPS
- gffread
- Salmon

**Option A: HPC Modules (Manual Setup)**
If you are running this on an HPC cluster, you can simply clone the repository and load your local modules before executing the pipeline. This is how the pipeline was originally developed.

```bash
git clone https://github.com/rodrie93/MAD-GE.git
cd MAD-GE

# Example module loading (adjust to your cluster's naming conventions)
module load BEDTOOLS
module load BCFtools

# Make sure all required tools are are in your PATH
```

**Option B: Conda/Mamba**
```bash
# Clone the repository
git clone https://github.com/rodrie93/MAD-GE.git
cd MAD-GE

# Create and activate the environment
conda env create -f environment.yml
conda activate mut_pipe
```

💻 OS Compatibility
MAD-GE is written using Bash and Python. It runs natively on Linux and macOS. Windows users must use the Windows Subsystem for Linux (WSL) (e.g., Ubuntu on WSL2) to execute the .sh wrappers.


📁 Directory Structure
MAD-GE expects the following organizational structure to run smoothly across multiple patient samples. Ensure your working directory looks like this before executing the pipeline:

```text
MAD-GE/
├── environment.yml         # Conda environment configuration
├── codes/                  # Pipeline scripts
│   ├── prep_reference/
│   │   ├── prepare_reference.sh
│   │   └── prep_reference.py
│   └── run_sample/
│       ├── mut_generation.py
│       ├── run_sample.sh
│       └── salmon_mut_quant.sh
├── ref/                    # Reference data
│   ├── genome.fa
│   ├── annotations.gff
│   └── ... (other generated ref files will appear here)
└── samples/                # Patient data, tumor tissue RNA fastq's and WES-derived VCF
    ├── P001/
    │   ├── P001-T-R1.fastq.gz
    │   ├── P001-T-R2.fastq.gz
    │   └── P001-WES.vcf.gz
    └── P002/
        └── ...
```

## Step 1: Prepare the Reference (Run Once)
Before processing individual samples, you must build the customized reference and the base Salmon index.
The transcriptome I used was the MANE subset (Matched Annotation from NCBI and EMBL-EBI).

Navigate to your reference folder:

```bash
cd ref/
```

Run the preparation script. Ensure your paths to the genome and annotations are correct, or are in the same folder.

```bash
bash ../codes/prep_reference/prepare_reference.sh path/to/genome.fa path/to/annotations.gff
```
This should take less than 10 minutes

## Output (1):
The "prepare_reference.sh" script will produce the following output files:
1. "new_myfmt.gff": as in "new my_format gff", is a file that I use downstream to build the mutant isoforms. Each row corresponds to an exon, in order (with the :# suffix indicating its sequential order in the isoform), with its genomic coordinates, strand, and exon/transcript/gene Ensembl IDs.

2. "tx_exon_map.tsv": as in "transcript-exon map", is a similar file as above but it is "isoform-centered". Each row is an isoform, where:
    2.1 The first column is its ENST ID, 
    2.2 Second column is the sequential order of exons conforming it, and 
    2.3 Third column is the corresponding ENSG (gene) ID

3. "exons.final.fasta": is a fasta file with all reference (or wild-type) exon sequences. This is used downstream in conjunction with the previous files to sequentially concatenate exons and build an isoform.

4. "ref_MANE_transcriptome.fa": is the transcriptome fasta, with only MANE-subset genes (should be around ~19k transcripts).

5. MANE-index/: is the Salmon-generated index of the MANE transcriptome fasta from above. This is used to quantify "traditional" gene expression counts, which disregards reads mutations.

## Step 2: Run Sample Pipeline

For each sample, MAD-GE will filter somatic calls, mutate the reference exons, build personalized isoforms, and quantify allele expression.

You can loop this across all your sample directories:

```bash
cd ../samples/ 
for sample in *; do
    cd $sample
    # 1. Generate the personalized mutant transcriptome
    bash ../../codes/run_sample/run_sample.sh

    # 2. Run quantification using Salmon (make sure run_sample.sh finishes first)
    # Make sure to setup your SBATCH options in the header of this script if using an HPC
    # Otherwise, just "bash" it
    sbatch ../../codes/run_sample/salmon_mut_quant.sh

    cd ..

```

## Output (2)

For each sample directory, after running "run_sample.sh", MAD-GE generates:
- Expanded Transcriptome FASTA: The personalized reference containing both wild-type and mutated isoforms.
- Mutation Annotations (.tsv): Detailed genomic and transcriptomic coordinates of all applied mutations.

After running the "salmon_mut_quant.sh" script, MAD-GE will output:
- Salmon Quantifications: Folders containing the traditional (wild-type only) and mutant-aware expression estimates (quant.sf) for downstream analysis. The traditional quantification will be found in the "MANE_trad" folder, while the mutant-aware counts will be stored at the "quants_muts" folder. These are outputs from Salmon, so it is advised to get used to their tool outputs and documentation as well.