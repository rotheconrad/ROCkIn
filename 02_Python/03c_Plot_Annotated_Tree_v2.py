#!/usr/bin/env python

'''Plot phylogenetic tree with metadata annotations

This script takes a newick file and the *_annotated.tsv file output from
02g_Tree_Distance_Cluster.py and creates a phylogenetic tree plot with
the clade cluster, genus classification, and functional annotation
information layer added on top of the tree.

Inputs:
    - phylogenetic tree in newick format
    - corresponding output from 03c_Tree_Distance_Cluster.py

Outputs:
    - PDF file of phylogenetic tree with metadata layers

Dependiencies:
    - python 3.7+
    - matplotlib (https://matplotlib.org/)
    - pycircos (https://github.com/ponnhide/pyCircos)
        * pip install python-circos

* there are 58 default colors, but users can also provide custom colors

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Oct 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
import pycircos
import matplotlib.pyplot as plt
from collections import defaultdict


def parse_annotations(infile):
    # parses the annotation file into pycircos friendly format.
    # Returns a list of lists used to build pycircos color dictionaries.

    # intialize lists to store metadata
    # names = leaf node names used to link colors to the tree.
    # cluster = cluster numbers. Used to build color dictionary
    # genes = gene annotation. Used to build color dictionary
    # genus = taxonomic classification. Used to build color dictionary
    # phylum = taxonomic classification. Used to build color dictionary
    # tclass = class taxonomic classification. Used to build color dictionary
    # status = reviewed(swissprot)/unreviewed(TrEMBL). Used to build color dict
    name, clust = [], []
    refs = defaultdict(list)

    # read through annotations.tsv file and populate dictionaries
    print("\n\t\tCluster", "RefSeq")
    with open(infile, 'r') as f:
        header = f.readline()
        for line in f:
            X = line.rstrip().split('\t')
            gname = X[0] # gene name
            name.append(gname)
            cluster = int(X[1]) #hdbscan=1, dbscan=2, kmeans=3
            clust.append(cluster)
            annot = X[2] # gene annotation
            if annot == "REFERENCE SEQUENCE":
                refs[cluster].append(gname)
                print(f'\t\t{cluster}:\t{gname}')

    # throw lists into master list for easy transfer
    data = [name, clust, refs]

    return data


def parse_colors(cfile):
    # input file should be three column tsv file
    # with columns for: gene name, color, refseq
    # where gene name is the terminal leaf node label
    # color is a hex value (ex:#80b1d3) corresponding to the assigned clade
    # and refseq is a 0 for non refseqs and 1 for a reference sequence
    # do not include a header
    color_dict = {}
    legend_dict = defaultdict(list)
    with open(cfile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            gname = X[0]
            color = X[1]
            ref = [2]
            color_dict[gname] = {
                    "color": color, "size": 12, "linewidth": 0.1,
                    "edgecolor": "#303030"
                    }
            if ref == "1":
                legend_dict[color].append(gname)



    return color_dict


def default_colors(data):
    # takes input lists of [leaf_node_names] and [metadata] and returns
    # a dict for pycircos colors and a dict for matplotlib legend.

    # set variables from input
    names = data[0]
    clusters = data[1]
    # number cluster colors by sorted unique values. used for legend
    clust_sort = sorted(list(set(data[1])))
    ref_dict = data[2]

    # default color set
    colors = [
            '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
            '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695', '#7f3b08',
            '#b35806', '#e08214', '#fdb863', '#fee0b6', '#d8daeb', '#b2abd2',
            '#8073ac', '#542788', '#e6f5d0', '#b8e186', '#7fbc41', '#4d9221',
            '#276419', '#8e0152', '#c51b7d', '#de77ae', '#f1b6da', '#fde0ef',
            '#c7eae5', '#80cdc1', '#35978f', '#01665e', '#003c30', '#543005',
            '#8c510a', '#bf812d', '#dfc27d', '#f6e8c3', '#e0e0e0', '#bababa',
            '#878787', '#4d4d4d', '#1a1a1a', '#66c2a5', '#fc8d62', '#8da0cb',
            '#e78ac3', '#a6d854'
            ]

    # generate a key to retrieve the same color for each unique meta value
    #color_key = {}
    # check if enough default colors for number of clusters
    if len(clust_sort) > len(colors):
        print('ERROR: More clusters than default colors.')
        exit(1)
    # build the color key
    #for j, clst in enumerate(clust_sort):
    #    color_key[clst] = colors[j]

    color_key = {clst: colors[j] for j, clst in enumerate(clust_sort)}

    # iterate names/clusters, get default colors, generate pyCircos dict
    color_dict = {}

    for name, cluster in zip(names, clusters):
        color = color_key[cluster]
        color_dict[name] = {
                    "color": color, "size": 12, "linewidth": 0.1,
                    "edgecolor": "#303030"
                    }

    # create dict for legend
    ref_key = {}
    for cl, rf in ref_dict.items():
        ref_seqs = '\n'.join(rf)
        cl_color = colors[clust_sort.index(cl)]
        ref_key[ref_seqs] = cl_color

    legend_dict = [color_key, ref_key]

    return color_dict, legend_dict


def build_plot(newick, colors, legends, outfile):

    # initialize tree figure
    Tarc    = pycircos.Tarc
    Tcircle = pycircos.Tcircle
    tarc = Tarc(tree=newick, format="newick", interspace=1)
    tcircle = Tcircle(figsize=(12,12))

    # build the plot
    tcircle.add_tarc(tarc)
    tcircle.set_tarcs()
    tcircle.plot_tree(
                    tarc.arc_id, rlim=(0,550),
                    cladevisual_dict=colors,
                    linewidth=0.4, linecolor="#606060"
                    )

    # save the plot
    tcircle.figure.savefig(outfile)

    # build legends
    outpre = '.'.join(outfile.split('.')[:-1])
    fnames = ['Clusters', 'RefSeqs']
    for i, legend in enumerate(legends):

        fig, ax = plt.subplots(figsize=(10,10))
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        for label, color in legend.items():
            ax.bar(
                0, 0, color=color, label=label, linewidth=0
                )

        ax.legend(
            title=f"{fnames[i]} Legend", title_fontsize='xx-large', loc="center",
            frameon=False, markerscale=5, fontsize='xx-large', ncol=1
            )

        plt.savefig(f'{outpre}_Legend_{fnames[i]}.pdf')
        plt.close()
    

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-a', '--annotation_input_file',
        help='Please specify the annotation input file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-n', '--newick_input_file',
        help='Please specify the newick input file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='Please specify the output file name!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c1', '--cluster_colors_file',
        help='(OPTIONAL) Please specify the clusters color file!',
        metavar=':',
        type=str,
        required=False,
        default=None
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\nRunning Script...\n')
    
    # define the input params
    annotations = args['annotation_input_file']
    newick = args['newick_input_file']
    outfile = args['output_file_name']
    c1 = args['cluster_colors_file']

    # read the annotations file into dictionaries to be used by pyCircos
    # data is a list of dicts [clust, genes, genus, phylum, tclass, status]
    # each dict contains {leaf_node_label: meta_data_value}
    print('\n\nParsing annotation input file ...')
    data = parse_annotations(annotations)
    # match up the input metadata with colors
    # only default colors for now but add custom colors in the future.

    print('\n\nHandling metadata color annotations ...')

    if c1:
        print(f'\n\nCollecting colors from {c1} ...')
        colors, legends = parse_colors(c1)
    else:
        print(f'\n\nGetting default colors for clusters ...')
        colors, legends = default_colors(data)

    # pass the info to PyCircus and build plot.
    print('\n\nPlotting ... and scheming ... ')

    _ = build_plot(newick, colors, legends, outfile)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()
