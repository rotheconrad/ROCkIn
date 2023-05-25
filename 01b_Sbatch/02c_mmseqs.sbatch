#!/bin/bash
#SBATCH -J mmseqs
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem=10G
#SBATCH -t 12:00:00
#SBATCH -o 02c_mmseqs.out
#SBATCH -e 02c_mmseqs.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# make a tmp folder
mkdir ${oDir}/tmp

# create an mmseqs database
DB=02_mmseqsDB
mkdir ${oDir}/$DB
mmseqs createdb $infile ${oDir}/$DB/$DB --dbtype 1

# Cluster input sequences at 90% amino acid sequence identity
C90=02_mmseqs_results
mkdir ${oDir}/$C90
mmseqs cluster ${oDir}/$DB/$DB ${oDir}/$C90/$C90 ${oDir}/tmp --min-seq-id ${ss} --cov-mode 1 -c 0.5 --cluster-mode 2 --cluster-reassign --threads 2
# convert to tsv format
mmseqs createtsv ${oDir}/$DB/$DB ${oDir}/$DB/$DB ${oDir}/$C90/$C90 ${oDir}/$C90.tsv --threads 2
# for amino acid clustering get the representative sequence for annotation.
C90rep=02_mmseqs_reps
mkdir ${oDir}/$C90rep
mmseqs createsubdb ${oDir}/$C90/$C90 ${oDir}/$DB/$DB ${oDir}/$C90rep/$C90rep
mmseqs convert2fasta ${oDir}/$C90rep/$C90rep ${oDir}/${C90rep}.fasta

rm -r ${oDir}/tmp

## Run Log ##
## sbatch --export infile=,oDir=,ss=0.90 /Path/to/GitHub/repo/01b_Sbatch/02c_mmseqs.sbatch
