#!/usr/bin/env python

'''Filter Blast+ Tabular Output by besthit, match length and percent identity.

This script filters tabular Blast+ output for best hit based on
the bitscore, as well as a user defined percent match length,
and percent identity of the sequence alignment.

Percent match length = alignment length / query sequence length.

This script randomizes the selection of tied matches by default or
removes tied matches all together with the -rtm option.

This script also outputs histograms for:
    1) 'pid': 'Sequence Identity Histogram (%)',
    2) 'alen': 'Alignment Length Histogram (AA)',
    3) 'pml': 'Sequence Alignment Ratio (%)',
    4) 'qlen': 'Query Length Histogram (AA)'

This script reads passing blast matches into memory
RAM usage depends on file size of passing matches.
RAM requirement could be close to the size of the tabular blast file.

This tool takes the following input parameters:

    * Tabular Blast+ output of the following -outfmt:
      '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

    * percent_match_length to filter for as decimal (ex: 0.5 for 50%)
      This is the length of the aligned query sequence to the reference sequence.
      This is used to remove sequences that don't align well or short spurious sequences.

    * percent_identity to filter for as percent (ex: 35 for 35%) 
      this is the sequence similarity of the aligned sequences
      This is used to remove sequences that aren't very similar even though
      they have a long alignment.

    * To filter for best hit only set -pml and -pid to 0

This script returns the following files:

    * input_file.fltrdBstHts.blst
    * 4 pdf files for histograms

This script requires the following packages:

    * argparse
    * random

This file can also be imported as a module and contains the follwing 
functions:

    * tabular_BlastPlus_filter - This function coordinates the filtering.
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random
from collections import defaultdict
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

def best_hits(query, bitscore, d, line, dups):
    """ Filters the besthit based on bitscore """

    if query in d:
        dups += 1
        old_bitscore = float(d[query][0].split('\t')[11])

        if bitscore > old_bitscore:
            d[query] = [line]

        elif bitscore == old_bitscore:
            d[query].append(line)

    else:
        d[query] = [line]

    return d, dups


def tabular_BlastPlus_filter(infile, pml, pid):

    d = {} # initialize dictionary for bitscore besthits
    dups = 0 # counter for number of duplicate matches for reads
    fails = 0 # counter for number of matches failing filters
    passes = 0 # counter for number of matches passing filters
    total = 0 # counter for total blast entries in file

    with open(infile, 'r') as f:

        for l in f:
            total += 1
            X = l.rstrip().split('\t')
            query = X[0] # read identifier
            bitscore = float(X[11]) # bitscore
            pID = float(X[2]) # percent sequence identity
            aLen = int(X[3]) # read alignment length
            qLen = int(X[12]) # full length of read
            pMatch = aLen / qLen # percent match length of read length
            if pID >= pid and pMatch >= pml:
                d, dups = best_hits(query, bitscore, d, l, dups)
                passes += 1
            else:
                fails += 1

    print('Total number of entries in blast file:', total)
    print('Number of entries failing the filters:', fails)
    print('Number of entries passing the filters:', passes)
    print('Number of duplicate blast matches passing filter to remove:', dups)

    return d


def parse_blast(file):
    """ parse blast file for pidents returns list of floats """

    data = {'pid': [], 'alen': [], 'pml': [], 'qlen': []}

    with open(file, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')
            pident = float(X[2])
            alen = int(X[3])
            qlen = int(X[12])
            pml = alen / qlen

            data['pid'].append(pident)
            data['alen'].append(alen)
            data['pml'].append(pml)
            data['qlen'].append(qlen)

    return data


def plot_hist(data, outpre, key):

    # Define the titles
    plot_titles = {
                    'pid': 'Sequence Identity Histogram (%)',
                    'alen': 'Alignment Length Histogram (AA)',
                    'pml': 'Sequence Alignment Ratio (%)',
                    'qlen': 'Query Length Histogram (AA)'
                    }

    print(f'\t\tPlotting {plot_titles[key]}...')

    # Set the colors
    bar_color = '#2171b5'
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    alpha = 0.6

    # Build the plot
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot titles
    ax.set_title(
        plot_titles[key],
        fontsize=20, y=1.02
        )

    # Plot labels
    ax.set_ylabel('Count', fontsize=14)
    ax.set_xlabel('Value', fontsize=14)

    # Set plot/grid style
    ax.minorticks_on()
    ax.tick_params(
        which='minor', axis='both', left=False, bottom=False
        )
    ax.tick_params(
                which='major', axis='both',
                left=True, bottom=True,
                size=6, width=2, tickdir='inout',
                labelsize=12, zorder=10
                )
    ax.yaxis.grid(
        which="minor", color=gridm, linestyle='--',
        linewidth=1, alpha=0.6, zorder=1
        )
    ax.yaxis.grid(
        which="major", color=gridM, linestyle='--',
        linewidth=1.5, alpha=0.4, zorder=1
        )
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_linewidth(2)

    # Plot the data
    ax.hist(
        data[key],
        bins=30,
        #orientation='horizontal',
        rwidth=0.9,
        color=bar_color,
        alpha=alpha,
        )

    # Set plot axis ranges
    #ax.set_xlim(left=0, right=int((max(d['xs'])+min(d['xs']))))

    # adjust layout, save, and close
    #plt.gca().invert_yaxis()
    fig.set_tight_layout(True)
    plt.savefig(f'{outpre}_{key}_histogram.pdf')
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify the tabular magic blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pml', '--percent_match_length',
        help='(Optional) Percent match length to filter for (Default = 0.5).',
        metavar='',
        type=float,
        required=False,
        default=0.5
        )
    parser.add_argument(
        '-pid', '--percent_identity',
        help='(Optional) Percent identity to filter for (Default = 30).',
        metavar='',
        type=float,
        required=False,
        default=30
        )
    args=vars(parser.parse_args())

    # define input parameters
    infile = args['in_file']
    pml = args['percent_match_length']
    pid = args['percent_identity']

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # get the filtered best hits
    filtered_best_hits = tabular_BlastPlus_filter(infile, pml, pid)

    # Write output file
    outfile = infile.split('.')[0] + '_fltrdBstHts.blst'
    outfile100 = infile.split('.')[0] + '_pID100.blst'
    pID100 = defaultdict(list)

    with open(outfile, 'w') as o, open(outfile100, 'w') as o100:

        for k,v in filtered_best_hits.items():
            # write 100% ID's to separate file
            pid = float(v[0].split('\t')[2])
            if pid == 100:
                for line in v:
                    X = line.split('\t')
                    query = X[0].split('//')[0]
                    subject = X[1]
                    pID100[subject].append(query)
                    o100.write(line)

            o.write(random.choice(v))
            
        print(
            'Number of best hit entries written to new file:',
            len(filtered_best_hits), '\n\n'
            )

    pID100out = infile.split('.')[0] + '_pID100_list.txt'
    with open(pID100out, 'w') as out:
        for subject, queries in pID100.items():
            qs = ','.join(queries)
            out.write(f'{subject}\t{qs}\n')

    data = parse_blast(outfile)
    outpre = infile.split('.')[0] + '_fltrdBstHts'
    for key in ['pid', 'alen', 'pml', 'qlen']:
        _ = plot_hist(data, outpre, key)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
