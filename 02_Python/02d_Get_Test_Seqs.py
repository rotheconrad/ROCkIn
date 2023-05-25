#!/usr/bin/env python

''' Get secondary representative from MMSeq2 clusters

MMSeqs2 has a method to retrieve the representative sequence form each cluster.
But, to test our models we want a fresh set of sequences not used in training,
we want a test set. For this, we will select a secondary representative from 
the MMSeqs 2 clusters.

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

import argparse
from collections import defaultdict

def parse_mmseqs_cluster_tsv(infile):

    '''
    The mmseqs cluster tsv file is two columns.
    column 1 is the cluster representative sequence.
    This sequence is repeated for every sequence in the cluster.
    Column 2 are sequences in the cluster.
    The general idea with this function is to create 2 dicionaries.
    One dictionary stores data by the cluster represenative name as in
    {cluster: genome}. Which genomes are in each cluster.
    The other dicitonary stores data by the genome name.
    {genome: cluster}. Which clusters are in each genome.
    '''

    print('\n\nParsing MMSeqs2 Cluster TSV file ...')
    # initialize variables
    # Store data in two dictionaries.
    byCluster = defaultdict(lambda: defaultdict(int))
    byGenome = defaultdict(lambda: defaultdict(int))

    # read through the file and populate the dictionaries
    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            cluster = X[0]
            genome = X[1].split('_')[0]
            # store data in dict of dict
            byCluster[cluster][genome] += 1
            byGenome[genome][cluster] += 1

    return byCluster, byGenome


def main():
    #byGenome = defaultdict(lambda: defaultdict(int))

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--mmseqs_cluster_tsv_input_file',
        help='Please specify the mmseqs cluster tsv input file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--secondary_representative_test_IDs',
        help='Please specify the name to use for the output file!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define the input parameters
    infile = args['mmseqs_cluster_tsv_input_file']
    outfile = args['secondary_representative_test_IDs']

    # parse the input file
    _ = get_secondary_mmseqs_cluster_reps(infile, outfile)

    print('\n\nComplete success space cadet! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

