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

# Make output directory
if [ ! -d 03_ROCkIn_Results ]; then mkdir 03_ROCkIn_Results; fi

# Load your servers clustalo module or clustalo must be available in your PATH
module load clustal-omega # change this for your server

if [ ! -s ${verified}.aln ]
  then
        n1=`basename ${verified}`
        clustalo -i $verified -o 03_ROCkIn_Results/03a_${n1}.aln --threads 12
fi

n2=`basename ${newseqs} | cut -d_ -f2-`
clustalo -i $newseqs -o 03_ROCkIn_Results/03a_${n2}.aln --profile1 03_ROCkIn_Results/03a_${n1}.aln --threads 12

## Run Log:
## sbatch --export verified=,newseqs= /Path/to/GitHub/repo/01b_Sbatch/02d_seq_alignment.sbatch

