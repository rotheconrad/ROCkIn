#!/bin/bash
#SBATCH -J EBI-DBFETCH
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH -t 24:00:00
#SBATCH -o 01b_ebi_dbfetch.out
#SBATCH -e 01b_ebi_dbfetch.err

cd $SLURM_SUBMIT_DIR

dbfetch='/Path/to/GitHub/repo/02_Python/01b_ebi_dbfetch.py'

while read p
	do
	name=`echo $p | cut -d: -f2`
	python ${dbfetch} fetchData $p > ${odir}/${gene}/${name}.dat;
	echo $name
 done < $input

## single examples/test
##  python ${dbfetch} fetchData tr:A0A2P5GGK6_9ENTR > MCR-9.1_Escherichia_coli.ids.A0A2P5GGK6_9ENTR.dat

## Run Log:
## sbatch --export input=,odir=,gene=  /Path/to/GitHub/repo/01b_Sbatch/01b_ebi_dbfetch.sbatch
