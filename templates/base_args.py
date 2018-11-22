#!/usr/bin/env python3
"""

Author: Hidde Bleeker (931202071020)
"""

# Imports
import argparse
from sys import argv


def parse_arguments():
    """Create an argument parser and parse arguments from argv

    :return: Parsed arguments.
    :rtype: dict, with {parameter_name: parameter_value}
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.fa', '.fasta'])
                        else parser.error("Wrong input file extension."),
                        help="Input FASTA file, with .fasta or .fa file "
                             "extension.")
    parser.add_argument('-l', '--linkage', default='complete', # required=True,
                        choices=['single', 'complete', 'average'],
                        help="Specify the hierarchical clustering "
                             "linkage type.")
    if len(argv) <= 1:
        parser.print_help()
        exit(1)
    else:
        return vars(parser.parse_args())


def parse_fasta(lines):
    """Parse a FASTA file. Sequences are stored by ID in a dictionary.

    :param lines: iterator of lines (e.g. in a FASTA file).
    :returns: generator of dictionary of {label:dna_seq}
    """
    label, seq = "", []
    for line in lines:
        if not line.strip() or line.startswith(';'):
            continue
        if line.startswith('>'):
            if label:
                yield label, "".join([subseq for subseq in seq])
            label = line.strip()[1:]
            seq = []
        else:
            seq.append(line.strip())
    yield label, "".join(seq)  # seq is already an iterator


if __name__ == '__main__':

    args = parse_arguments()
    print(type(args), args)

    with open(args['assembled']) as inp_f:
        INPUT_SEQ_DICT = {_id: (_seq, len(_seq)) for _id, _seq in
                          parse_fasta(inp_f)}
