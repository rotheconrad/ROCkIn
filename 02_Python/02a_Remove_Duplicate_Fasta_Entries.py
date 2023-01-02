#!/usr/bin/env python

''' Removes duplicate fasta entries in a fasta file.

Multiple sequence searches with similar sequences using Blast to NCBI,
UniProt, or any other database can return overlapping results.

To remove duplicates. Concatenate fasta files from the Blast matched
sequences and input to this program.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January, 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
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


def remove_duplicate_fasta(infile, outfile):
    """Reads files and writes fasta output of matching reads"""

    count = 0
    dedup_fasta = {}
    duplicates = 0

    with open(infile, 'r') as file:
        for name, seq in read_fasta(file):
            count += 1
            if name in dedup_fasta:
                duplicates += 1
            else:
                dedup_fasta[name] = seq

    with open(outfile, 'w') as out:
        for name, seq in dedup_fasta.items():
            out.write(f'{name}\n{seq}\n')

    print(f'\n\t\tTotal sequences in file: {count}')
    print(f'\t\tDuplicates Removed: {duplicates}')
    print(f'\t\tUnique sequences retained: {len(dedup_fasta)}')

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f', '--input_fasta_file',
        help='Please specify the input fasta file!',
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
    _ = remove_duplicate_fasta(
                                args['input_fasta_file'],
                                args['out_file']
                                )

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
