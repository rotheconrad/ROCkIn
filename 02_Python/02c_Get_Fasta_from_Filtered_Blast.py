#!/usr/bin/env python

''' Retrieves fasta sequences matching Filtered Blast Output.

For retrieving fasta sequences that map to the reference sequence.

Takes query fasta file and tabular blast output and returns
the sequences matching the blast file. The blast output should be
filtered prior to running this script for best hit and any
desired pIdent or match length cutoffs.

Should work with Blast+ or MagicBlast tabular outputs.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: October 29th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


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


def get_fasta_sequence(blast, query, outfile):
    """Reads files and writes fasta output of matching reads"""

    blastmatch = {}

    with open(blast, 'r') as b:
        for l in b:
            qry = l.split('\t')[0]
            blastmatch[qry] = ''

    with open(query, 'r') as q, open(outfile, 'w') as o:
        for name, seq in read_fasta(q):
            if name[1:].split(' ')[0] in blastmatch:
                o.write(f'{name}\n{seq}\n')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--tabular_blast_file',
        help='Please specify the filtered tabular blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-q', '--query_fasta_file',
        help='Please specify the query fasta file!',
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
    get_fasta_sequence(
                    args['tabular_blast_file'],
                    args['query_fasta_file'],
                    args['out_file']
                    )


if __name__ == "__main__":
    main()
