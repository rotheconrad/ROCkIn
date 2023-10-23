#!/usr/bin/env python

''' Distance clustering from phylogenetic tree in newick format.

Normal Use:
Input the newick file from IQ-Tree generated in the previous step.
Input the concatenated fasta file of the RepSeqs and the cluster
representatives (these should match the sequences used to build the
multiple alignment the tree is based on)

Experimental/Development Use:
Rerun with the distance matrix to experiment with clustering params.
Re-input the distance matrix on subsequent runs instead of the newick.
Fasta file is optional.

Normal Use:
Converts newick file to a distance matrix and clusters using HDBSCAN.

Output:
Distance matrix as tab separated file (TSV).
Leaf node names with cluster label as TSV file with species and gene
annotations if optional fasta file is provided.

Required Python Packages:
conda install -c conda-forge biopython
conda install -c conda-forge hdbscan
conda install -c intel scikit-learn
* scikit-learn installs with hdbscan

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February, 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, itertools
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import Phylo
from collections import defaultdict
import hdbscan
import sklearn.cluster as cluster
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt


def newick_to_distmat(newick, outpre):
    """Reads in a newick file and returns a distance matrix as pd df"""

    # Read in and parse the newick file with BioPython Phylo package
    nwk = Phylo.read(newick, 'newick')

    # initialize dict of dicts to store distance matrix
    dstmt = defaultdict(lambda: defaultdict(list))
    # loop over all points in the tree and calculate distances
    for x, y in itertools.combinations(nwk.get_terminals(), 2):
        dst = float(nwk.distance(x, y))
        dstmt[x.name][y.name] = dst
        dstmt[y.name][x.name] = dst
    for x in nwk.get_terminals():
        dstmt[x.name][x.name] = 0
    
    # Create dataframe from dictionary
    df = pd.DataFrame(dstmt)
    # rearrange for symmetric matrix move the 1st column to the end.
    # the above loops place in the first column what should be the last column
    lst = df.columns.tolist() # get the column names
    lst.append(lst.pop(0)) # move the 1st column to last column
    df = df[lst] # select columns in the new order
    # write matrix to file
    df.to_csv(f'{outpre}_distmat.tsv', sep='\t')

    return df


def distance_matrix_cluster(df, outpre):
    """Takes a distance matrix and returns clusters"""

    # HDBSCAN
    print('\n\t\t\t\tRunning HDBSCAN algorithm ...')
    min_cluster = int(len(df) * 0.02) # min cluster sized based on input points
    if min_cluster < 2: min_cluster = 2 # min_cluster must be greater than 1
    epsilon = float(np.quantile(df.to_numpy().ravel(), 0.01)) # epsilong to 10% quantile
    hdb = hdbscan.HDBSCAN(
                        metric='precomputed',
                        min_cluster_size=min_cluster,
                        min_samples=1, # places more points in clusters
                        cluster_selection_epsilon=epsilon,
                        cluster_selection_method='leaf' # 'eom' (default) or 'leaf'
                        )
    hdb.fit_predict(df)

    '''
    # We explored these other algorithms but went with hdbscan
    # USE HDBSCAN clusters for n_clusters if not user provided
    if not n_clusters: n_clusters = hdb.labels_.max()

    # DBSCAN
    print('\n\t\t\t\tRunning DBSCAN algorithm ...')
    db = cluster.DBSCAN()
    db.fit_predict(df)

    # AgglomerativeClustering
    print('\n\t\t\t\tRunning Agglomerative algorithm ...')
    agl = cluster.AgglomerativeClustering(
                                        linkage='single',
                                        n_clusters=n_clusters,
                                        affinity='precomputed'
                                        )
    agl.fit_predict(df.to_numpy())

    # SpectralClustering
    print('\n\t\t\t\tRunning Spectral algorithm ...')
    spc = cluster.SpectralClustering(
                                    n_clusters=n_clusters,
                                    affinity='precomputed'
                                    )
    spc.fit_predict(df)

    # MeanShift
    print('\n\t\t\t\tRunning Mean Shift algorithm ...')
    ms = cluster.MeanShift(cluster_all=False)
    ms.fit_predict(df)

    # AffinityPropagation
    print('\n\t\t\t\tRunning Affinity Propagation algorithm ...')
    afp = cluster.AffinityPropagation(affinity='precomputed')
    afp.fit_predict(df)

    # KMeans
    print('\n\t\t\t\tRunning KMeans algorithm ...')
    km = cluster.KMeans(n_clusters=n_clusters)
    km.fit_predict(df)
    '''

    # Add the labels for each algorithm to the df.
    print('\n\t\t\t\tUpdating DataFrame and writing to file ...')
    lbld_data = {
                "Gene_Name": df.index,
                "HDBSCAN": hdb.labels_,
                #"DBSCAN": db.labels_,
                #"Agglomerative": agl.labels_,
                #"Spectral": spc.labels_,
                #"MeanShift": ms.labels_,
                #"Affinity": afp.labels_,
                #"K-Means": km.labels_,
                }
    df_lbld = pd.DataFrame(lbld_data)

    df_lbld.to_csv(f'{outpre}_labeled.tsv', sep='\t', index=False)

    return df


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


def add_gene_info(fasta, outpre):
    """Reads in UniProt Fasta and adds species and gene annotations"""

    # initialize dict for gene info parsed from fasta file
    gene_info = {}
    # initialize dict for fasta file output of each cluster
    byCluster = defaultdict(list)

    # parse gene info out of fasta file deflines
    # uniprot fasta has annotation and species info in the deflines
    with open(fasta, 'r') as file:
        for name, seq in read_fasta(file):
            X = name.rstrip().split('//')
            if len(X) > 1:
                uid, status, func, spec, tax = X[0][1:], X[1], X[2], X[3], X[4]
            else:
                R = 'REFERENCE SEQUENCE'
                uid, status, func, spec, tax = X[0][1:], R, R, R, R

            gene_info[uid] = [status, func, spec, tax, seq]

    # Read the labeled dataframe file and add gene info
    data = f'{outpre}_labeled.tsv'
    outfile = f'{outpre}_annotated.tsv'
    with open(data, 'r') as inf, open(outfile, 'w') as outf:
        header = inf.readline().rstrip()
        outf.write(
            f'{header}\tGene_function\tSpecies_name\tStatus\tTaxonomy\tSeq\n'
            )

        for line in inf:
            X = line.rstrip().split('\t')
            name = X[0]
            clstr = X[1] # select the kmeans row
            info = gene_info[name]
            X.extend(info)
            outf.write('\t'.join(X) + '\n')
            # add fasta sequence for cluster to byCluster dict.
            fa = f'>{name}\n{info[4]}\n'
            byCluster[clstr].append(fa)

    # write fasta file for each cluster
    for cluster, seqs in byCluster.items():
        with open(f'{outpre}_{int(cluster):03}.fa', 'w') as file:
            for entry in seqs:
                file.write(entry)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_newick_file',
        help='(Optional) Please specify the input newick file!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-d', '--input_distance_matrix_file',
        help='(Optional) Please specify the input distance matrix file!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='How do you want to name the output files?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f', '--uniprot_fasta_file',
        help='(Optional) Fasta file from Uniprot.',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...')

    # define input params
    newick = args['input_newick_file']
    distmat = args['input_distance_matrix_file']
    outpre = args['output_prefix']
    fasta = args['uniprot_fasta_file']

    # check for output directory and create if needed
    p = outpre.split('/')
    if len(p) > 1:
        Path('/'.join(p[:-1])).mkdir(parents=True, exist_ok=True)

    if newick:
        # generate distance matrix from newick file
        # returns a pandas dataframe
        print('\n\t\tParsing Newick File and Converting to Distance Matrix ...')
        df = newick_to_distmat(newick, outpre)
    elif distmat:
        # read in distance matrix as pandas dataframe
        print('\n\t\tParsing Newick File and Converting to Distance Matrix ...')
        df = pd.read_csv(distmat, sep='\t', index_col=0)

    else:
        # user needs to input newick or distmat
        print('\nPlease provide a newick or distance matrix file!!')

    # cluster the distance matrix
    print('\n\t\tClustering the Distance Matrix ...')
    df_clstrd = distance_matrix_cluster(df, outpre)

    # Add species and gene annotations to dataframe from Uniprot fasta file
    if fasta:
        print('\n\t\tAdding species and gene annotations to DataFrame ...')
        _ = add_gene_info(fasta, outpre)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
