#!/bin/bash
#SBATCH -J TrimAl
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH -t 24:00:00
#SBATCH -o 02e_seq_trim.out
#SBATCH -e 02e_seq_trim.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load your servers trimal module, make trimal available in your path, or provide the path to trimal as a variable "trimal"
trimal=/Path/to/your/trimal/source/./trimal

${trimal} -in ${input} -out ${output} -resoverlap 0.75 -seqoverlap 80 -automated1

## Param key:
## -resoverlap: Minimum overlap of a positions with other positions in the column to be considered a "good position". Range: [0 - 1]. (see User Guide).
## -seqoverlap: Minimum percentage of "good positions" that a sequence must have in order to be conserved. Range: [0 - 100](see User Guide).

## Run Log ##
## sbatch --export input=fasta.aln,output=fasta.aln /Path/to/GitHub/repo/01b_Sbatch/02e_seq_trim.sbatch