#!/usr/bin/env python

'''Plot phylogenetic tree with metadata annotations

This script takes a newick file and the *_annotated.tsv file output from
02g_Tree_Distance_Cluster.py and creates a phylogenetic tree plot with
the clade cluster, genus classification, and functional annotation
information layer added on top of the tree.

Inputs:
    - phylogenetic tree in newick format
    - corresponding output from 02g_Tree_Distance_Cluster.py

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
            '#a6cee3', '#b2df8a', '#fd8d3c', '#ccebc5', '#fee391', '#fe9929',
            '#e31a1c', '#4292c6', '#238b45', '#9ebcda', '#8c96c6', '#67001f',
            '#bc80bd', '#f16913', '#d94801', '#00441b', '#80b1d3', '#fdd0a2',
            '#d9d9d9', '#662506', '#4d004b', '#1f78b4', '#66c2a4', '#fdbf6f',
            '#238b45', '#d4b9da', '#cab2d6', '#fccde5', '#fdb462', '#525252',
            '#33a02c', '#bdbdbd', '#fb6a4a', '#6a3d9a', '#000000', '#2171b5',
            '#9e9ac8', '#00441b', '#810f7c', '#99d8c9', '#bcbddc', '#ff7f00',
            '#88419d', '#fec44f', '#807dba', '#737373', '#a1d99b', '#8dd3c7',
            '#d9d9d9', '#ffffb3', '#df65b0', '#9ecae1', '#54278f', '#ffff99',
            '#252525', '#7f2704', '#41ae76', '#fc9272', '#c6dbef', '#6a51a3',
            '#74c476', '#b3de69', '#67000d', '#980043', '#b15928', '#006d2c',
            '#cc4c02', '#cb181d', '#fb8072', '#8c6bb1', '#fb9a99', '#ec7014',
            '#3f007d', '#6baed6', '#c7e9c0', '#fdae6b', '#a50f15', '#a63603',
            '#bebada', '#c994c7', '#993404', '#ce1256', '#006d2c', '#fcbba1',
            '#bfd3e6', '#ef3b2c', '#e7298a', '#41ab5d', '#08306b', '#ffed6f',
            '#ccece6', '#969696', '#08519c'
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
        help='Please specify the clusters color file!',
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
