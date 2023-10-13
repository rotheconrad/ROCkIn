#!/bin/bash
#SBATCH -J Blastp
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem=2G
#SBATCH -t 12:00:00
#SBATCH -o 02b_BlastP.out
#SBATCH -e 02b_BlastP.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Set python script directory
scripts=/Path/to/GitHub/repo/02_Python

# Load Modules Needed
# load your systems blast-plus module or have blastp available in your PATH
module load blast-plus/2.10.1

# Make Blast DB
if [ ! -s ${ref}.pdb ]
  then
        makeblastdb -dbtype prot -in $ref
fi

# Make output directory
if [ ! -d 02b_filter ]; then mkdir 02b_filter; fi

# Run Blast
blastp -num_threads 2 -max_target_seqs 10 \
-db ${ref} -query ${qry} -out ${out}.blast  -subject_besthit -outfmt \
'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

# Filter Blast and plot hists
# 1) Removes matches <= 30% identity to verified sequences
# 2) Remove matches with <= 50% sequence alignment (alignment length / query sequence length).

python ${scripts}/02b_Blastp_filter_hist.py -i ${out}.blast -o ${out}

# Get Blast Filtered Fasta Sequence
python ${scripts}/02c_Get_Fasta_from_Filtered_Blast.py -b ${out}.fltrdBstHts.blst -q ${qry} -o ${out}_blast_fltrd.fasta
python ${scripts}/02c_Get_Fasta_from_Filtered_Blast.py -b 02b_filter/${out}_pID100.blst -q ${qry} -o 02b_filter/${out}_pID100.fa

## Run Log ##
## sbatch --export ref=OG.faa,qry=derep_representatives.faa,out=out.blast /Path/to/GitHub/repo/01b_Sbatch/00b_BlastP.sbatch

