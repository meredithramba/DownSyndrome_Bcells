#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=igblast
#SBATCH --output=../../logs_updated/igblast_%A.out
#SBATCH --error=../../logs_updated/igblast_%A.err
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --partition=cpu
#SBATCH --mail-user=mr4527@cumc.columbia.edu
#SBATCH --mail-type=END,FAIL

# Load conda environment
source /groups/db3700_gp/mr4527/miniconda/etc/profile.d/conda.sh
conda activate complete_env

# ---- USER SETTINGS ----
DATA="../data"
OUTDIR="IgBLAST"
# -----------------------

# Pull the sample name for this array task
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$DATA/samples.txt")

AssignGenes.py igblast -s $DATA/${SAMPLE}/filtered_contig.fasta \
    -b /groups/db3700_gp/mr4527/igblast --organism human --loci ig \
    --format blast --outdir $OUTDIR --outname $SAMPLE --nproc 4

MakeDb.py igblast -i $OUTDIR/${SAMPLE}_igblast.fmt7 \
    -s $DATA/${SAMPLE}/filtered_contig.fasta \
    -r /groups/db3700_gp/mr4527/germlines/imgt/human/vdj/ \
    --10x $DATA/${SAMPLE}/filtered_contig_annotations.csv --extended \
    --outdir $OUTDIR --outname $SAMPLE --nproc 4

