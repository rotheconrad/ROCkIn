#!/bin/bash
#SBATCH -J ClustalO
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=12
#SBATCH --mem=10G
#SBATCH -t 4:00:00
#SBATCH -o 02d_ClustalO.out
#SBATCH -e 02d_ClustalO.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load your servers clustalo module or clustalo must be available in your PATH

if [ ! -s ${verified}.aln ]
  then
	clustalo -i $verified -o ${verified}.aln
fi

clustalo -i $newseqs -o ${newseqs}.aln --profile1 ${verified}.aln --threads 12

## Run Log:
## sbatch --export verified=,newseqs= /Path/to/GitHub/repo/01b_Sbatch/02d_seq_alignment.sbatch

