#!/bin/bash
#SBATCH -J EBI-BLAST
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH -t 24:00:00
#SBATCH -o 00a_ebi_blast.out
#SBATCH -e 00a_ebi_blast.err

cd $SLURM_SUBMIT_DIR

# path to ebi blast script
ebi_blast='/Path/to/GitHub/repo/02_Python/01a_ebi_ncbiblast.py'

# executes the command Run 30 sequences at a time from 1 fasta file.
python ${ebi_blast} --email YOUR_EMAIL --program blastp --stype protein \
--sequence $fasta --database uniprotkb --multifasta --useSeqId --maxJobs 30 --pollFreq 60 \
--outformat ids --exp 10 --alignments 1000

## Run Log:
## sbatch --export fasta=RefSeqs.faa /Path/to/GitHub/repo/01b_Sbatch/01a_ebi_blast.sbatch
