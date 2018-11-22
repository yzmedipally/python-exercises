#!/usr/bin/env python3
"""Script to

Author: Hidde Bleeker (931202071020)
"""

# Imports
import argparse
import subprocess as sbp
from sys import argv
from operator import itemgetter
from os.path import exists



def parse_arguments():
    """Create an argument parser and parse arguments from argv

    :return: Parsed arguments.
    :rtype: dict, with {parameter_name: parameter_value}
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.fa', '.fasta'])
                        else parser.error("Wrong input file extension."),
                        help="Input reference genome as FASTA file, with "
                             ".fasta or .fa extension.")
    parser.add_argument('assembled', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.fa', '.fasta'])
                        else parser.error("Wrong input file extension."),
                        help="Input assembled genome to compare to reference "
                             "as FASTA file, with .fasta or .fa extension.")
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


def find_size_and_n50(seq_dict):
    """Find the N50 length and N50 index

    :param dict seq_dict: ID (str) as keys, tuple of (sequence (str),
    sequence length (int)) as values.
    :return: Assembly size, N50 length, and N50 index.
    :rtype tuple (triplet), of (int, int, int)
    """
    size = sum(_v[1] for _v in seq_dict.values())
    # print(size)
    sorted_seqs = sorted([_val for _val in seq_dict.values()],
                         key=itemgetter(1))
    #
    # print(list(zip(*sorted_seqs))[1])
    # print(list(zip(*sorted_seqs))[1][45])

    size_50 = size // 2
    # print(size_50, sum(list(zip(*sorted_seqs))[1][45:]))

    current_size = 0
    for i, _seq in enumerate(sorted_seqs):
        current_size += _seq[1]
        if current_size >= size_50:
            return size, _seq[1], i + 1


def run_lastz(reference_ip, query_ip, output_format='general',
              output_fn='outlastz.txt'):
    """Run LASTZ, sequence alignment program to compare assembly to reference

    :param str reference_ip: Filename of the reference sequence (FASTA format)
    :param str query_ip: Filename of the query sequence (FASTA format)
    :param str output_format: format type for the lastz utility
    :param str output_fn: Filename to write output to.
    :return:
    """
    if exists(output_fn):
        print("Output file already exists; lastz not run.")
        return
    cmd = "lastz {} {} --format={} --output={}"\
        .format(reference_ip, query_ip, output_format, output_fn)

    print("Running LASTZ, command used:\n{}".format(cmd))
    try:
        sbp.check_call(cmd)
    except sbp.CalledProcessError as err:
        print("LASTZ failed to run. Error:\n{}\nQuitting.".format(err))
    # implicit return





if __name__ == '__main__':

    args = parse_arguments()
    # print(type(args), args)
    with open(args['assembled']) as inp_f:
        ASSEMBLY_DICT = {_id: (_seq, len(_seq)) for _id, _seq in
                         parse_fasta(inp_f)}
    with open(args['reference']) as inp_f:
        REFERENCE_DICT = {_id: (_seq, len(_seq)) for _id, _seq in
                          parse_fasta(inp_f)}

    print(find_size_and_n50(ASSEMBLY_DICT))
    print(find_size_and_n50(REFERENCE_DICT))

    run_lastz(args['reference'], args['assembled'])




    # for _val in ASSEMBLY_DICT.values():
    #     print(_val[1])
    # print("\n\n")
    # for _val in REFERENCE_DICT.values():
    #     print(_val[1])

    # first = next(iter(ASSEMBLY_DICT))
    # print(first, ASSEMBLY_DICT[first][1])
    # first = next(iter(REFERENCE_DICT))
    # print(first, REFERENCE_DICT[first][1])
