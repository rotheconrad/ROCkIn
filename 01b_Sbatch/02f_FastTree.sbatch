#!/bin/bash
#SBATCH -J FastTree
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH -t 24:00:00
#SBATCH -o 02f_FastTree.out
#SBATCH -e 02f_FastTree.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load your servers fasttree module or make FastTree available in your PATH

FastTree $input > $output

## Run Log ##
## sbatch --export input=fasta.aln,output= /Path/to/GitHub/repo/01b_Sbatch/02f_FastTree.sbatch

