#!/usr/bin/env python

''' Get secondary representative from MMSeq2 clusters

MMSeqs2 has a method to retrieve the representative sequence form each cluster.
But, to test our models we want a fresh set of sequences not used in training,
we want a test set. For this, we will select a secondary representative from 
the MMSeqs 2 clusters.

By default no random seed is set so the secondary representative selection is
chosen randomly from each cluster each time the script is rerun (excluding the
primary cluster representative selected by mmseqs2). use the optional
-s parameter to set a fixed seed for reproducible "random" selections.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: May 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random
from collections import defaultdict


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def get_secondary_mmseqs_cluster_reps(infile, outpre):

    '''
    The mmseqs cluster tsv file is two columns.
    column 1 is the cluster representative sequence.
    This sequence is repeated for every sequence in the cluster.
    Column 2 are sequences in the cluster.
    The general idea with this function is to create 2 dictionaries.

    '''

    print('\n\nParsing MMSeqs2 Cluster TSV file ...')
    # initialize variables
    # Store data in two dictionaries.
    data = defaultdict(list)
    secReps = {}
    tmpList = None

    # read through the file and populate the dictionaries
    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            cluster, seq = X[0].split('//')[0], X[1].split('//')[0]
            if seq == cluster:
                if tmpList: secReps[random.choice(tmpList)] = ''
                tmpList = []
            if seq != cluster:
                # store data in dict
                tmpList.append(seq)

    with open(f'{outpre}_secReps_IDlist.txt', 'w') as out:
        for ID, _ in secReps.items():
            out.write(f'{ID}\n')

    return secReps


def get_secRep_fasta(secReps, fastafile, outpre):

    fastaout = f'{outpre}_secReps.fa'
    with open(fastafile, 'r') as fasta, open(fastaout, 'w') as fout:
        for name, seq in read_fasta(fasta):
            n = name[1:].split('//')[0]
            if n in secReps:
                fout.write(f'{name}\n{seq}\n')

    return True


def main():
    #byGenome = defaultdict(lambda: defaultdict(int))

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-c', '--mmseqs_cluster_tsv_input_file',
        help='Please specify the mmseqs cluster tsv input file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f', '--fasta_file_input_to_mmseq',
        help='Please specify the fasta file input to mmseqs!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Prefix for output files of txt id list and fasta file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-s', '--random_seed',
        help='(OPTIONAL) set the random seed to set fixed random behavior!',
        metavar=':',
        type=str,
        required=False,
        default=None
        )
    args=vars(parser.parse_args())

    # define the input parameters
    clusterfile = args['mmseqs_cluster_tsv_input_file']
    fastafile = args['fasta_file_input_to_mmseq']
    outpre = args['output_file_prefix']
    rndsd = args['random_seed']

    random.seed(rndsd)

    # parse the input file
    secReps = get_secondary_mmseqs_cluster_reps(clusterfile, outpre)
    _ = get_secRep_fasta(secReps, fastafile, outpre)

    print('\n\nComplete success space cadet! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

