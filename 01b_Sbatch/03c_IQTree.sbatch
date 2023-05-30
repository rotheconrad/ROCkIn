#!/bin/bash
#SBATCH -J IQTree
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=20
#SBATCH --mem=50G
#SBATCH -t 5-00:00:00
#SBATCH -o 03c_IQTree.out
#SBATCH -e 03c_IQTree.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load your servers fasttree module or make FastTree available in your PATH
# IQTREE works well with a conda install

iqtree -s ${input} --prefix ${outpre} -m MFP -alrt 1000 -B 1000 -T THREADS

## Run Log ##
## sbatch --export input=fasta.aln,outpre= /Path/to/GitHub/repo/01b_Sbatch/03c_IQTree.sbatch

