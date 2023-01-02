#!/usr/bin/env python

''' Parse UniProt .dat file to fasta with species and function data

dbfetch messes up the fasta defline annotation details

So, I retrieve the .dat files and parse those to get the defline I want.

input is a directory containing .dat files
output is a single fasta file

deflines are formated as: >1_2_3_4_5_
        1) UniProt ID
        2) Status (reviewed or unreviewed)
        3) functional annotation
        4) species classification
        5) gene name

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April, 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, os


def get_fasta(f, ofile):
    # does the actual parsing of dat to fasta

    # keep track of taxonomy lines
    taxonomy = []
    # keep track of gene sequence
    sequence = []

    # check each line for the info we want
    for line in f:
        # get uniprot id, status and gene length
        if line.startswith('ID '):
            X = line.rstrip().split()
            uid, status, glen = X[1], X[2][:-1], int(X[3])
        # get the gene function
        elif 'SubName:' in line or 'RecName:' in line:
            func = line.rstrip().split('=')[1].split(' {')[0].replace(';', '')
        # get the species name
        elif line.startswith('OS '):
            spec = line.rstrip()[5:-1]
        # get taxonomy
        elif line.startswith('OC '):
            X = line[5:].rstrip().replace('.', '').replace(';', '').split()
            taxonomy.extend(X)
        # get the sequence
        elif line.startswith(' '):
            X = line.rstrip().split()
            sequence.extend(X)

    # put the pieces together
    tax = ';'.join(taxonomy)
    seq = ''.join(sequence)
    lineout = f'>{uid}//{func}//{spec}//{status}//{tax}\n{seq}\n'
    # write line to file
    ofile.write(lineout)

    # check gene length
    if len(seq) != glen: print(f'\n\nGene length discrepency: {uid}')

    return True


def parse_dat_file(idir, outfile):
    """ Opens files. reads. parses. writes. """

    # get a list of files from the input directory
    dfiles = [f for f in os.listdir(idir) if f.split('.')[-1] == 'dat']

    # open outfile
    with open(outfile, 'w') as ofile:
        # read through list of .dat files
        for file in dfiles:
            # open each .dat file
            with open(f'{idir}/{file}', 'r') as f:
                # convert it to a fasta file
                _ = get_fasta(f, ofile)
                
    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_directory',
        help='Please specify the directory with .dat files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')

    # define parameters
    idir = args['input_directory']
    outfile = args['out_file']

    # check idir format
    if idir[-1] == '/': idir = idir[:-1]

    # parse the dat file
    _ = parse_dat_file(idir, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
