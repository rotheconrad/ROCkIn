#!/bin/bash
#SBATCH -J RAxML-B
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH --mem=10G
#SBATCH -t 4:00:00:00
#SBATCH -o 02f_RAxML-B.out
#SBATCH -e 02f_RAxML-B.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load the modules for your servers RAxML or RAxML must be avialable in your PATH
module load raxml/8.2.11-mva2

# set some veriables.
pre=`pwd`/${output}
out=${name}.RAxML

if [ ! -d ${output} ]; then mkdir ${output}; fi

# Run RAxML for amino acid sequence
raxmlHPC-PTHREADS-AVX -s ${input} -n ${out} -w ${pre} -m PROTGAMMAAUTO -f a -N autoMRE -p 4564821 -T 8 -x 1235

## Run Log ##
## sbatch --export input=fasta.aln,output=directory /Path/to/GitHub/repo/01b_Sbatch/002f_RAxML_AminoAcid-Bootstrap.sbatch

