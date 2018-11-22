#!/usr/bin/env python3
"""Script to compare assembled genome to reference, with help of LASTZ aligner

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
    parser.add_argument('query', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.fa', '.fasta'])
                        else parser.error("Wrong input file extension."),
                        help="Input query genome to compare to reference "
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
    :return: None, command calling function
    """
    if exists(output_fn):
        print("## Output file already exists; lastz not run.")
        return
    cmd = "lastz {} {} --format={} --output={}"\
        .format(reference_ip, query_ip, output_format, output_fn)

    print("## Running LASTZ, command used:\n## {}".format(cmd))
    try:
        sbp.check_call(cmd, shell=True)
    except sbp.CalledProcessError as err:
        print("## LASTZ failed to run. Error:\n## {}\n## Quitting."
              .format(err))
    # implicit return


def parse_lastz_output(list_of_lines, header=True):
    """Parse the output file of LASTZ for alignments

    :param generator list_of_lines: LASTZ file, as list of lines to parse
    :param bool header: Indicating if input file has a header
    :return: generator of alignments, with coordinates (0-based, open-ended)
    parsed from lastz file
    :rtype generator
    """
    for line in list_of_lines:
        if header:
            header = False
            continue
        if line.strip():
            fields = line.strip().split('\t')
            yield int(fields[4]), int(fields[5])


def find_uncovered(seq_len, covered_coords):
    """Find the positions not present in sequence and return those as indices

    :param int seq_len: Length of the sequence
    :param iterator covered_coords: iterator of pairs of ints; these are the
    ranges which are present (or "covered").
    :return: list of pairs of the ranges that are uncovered by the input
    coordinates; returned as 0-based, half-open indices.
    :rtype generator of pairs of ints
    """
    uncovered = [True, ] * seq_len
    for coords in covered_coords:
        for i in range(coords[0], coords[1]):
            uncovered[i] = False

    curr = []
    for i, not_covered in enumerate(uncovered):
        if not_covered:
            if not curr:
                curr = [i, i + 1]
            else:
                curr[-1] += 1
        else:
            if curr:
                yield tuple(curr)
            curr = []


if __name__ == '__main__':
    # # In the case that argparse is not allowed:
    # if not len(argv) == 3:
    #     print("usage: python3 {} reference.fa query.fa".format(argv[0]))
    # # In this case every args[key] has thus to be replaced with argv[1] or
    # # argv[2].

    # Parse arguments and read in (parse) FASTA files
    args = parse_arguments()
    # print(type(args), args)
    with open(args['query']) as inp_f:
        QUERY_DICT = {_id: (_seq, len(_seq)) for _id, _seq in
                      parse_fasta(inp_f)}
    with open(args['reference']) as inp_f:
        REFERENCE_DICT = {_id: (_seq, len(_seq)) for _id, _seq in
                          parse_fasta(inp_f)}

    # Run LASTZ for sequence alignment
    run_lastz(args['reference'], args['query'], output_fn="outlastz.txt")

    # Find unaligned positions from LASTZ output file.
    with open("outlastz.txt") as lastz_f:
        ALIGNED_POSITIONS = list(parse_lastz_output(lastz_f))
    UNALIGNED_POSITION = find_uncovered(next(iter(REFERENCE_DICT.values()))[1],
                                        ALIGNED_POSITIONS)

    # Print the genome stats and the uncovered regions
    print("{}: TOTAL={}; N50 SIZE={}; N50 INDEX={}"
          .format(args['query'], *find_size_and_n50(QUERY_DICT)))
    print("{}: TOTAL={}; N50 SIZE={}; N50 INDEX={}"
          .format(args['reference'], *find_size_and_n50(REFERENCE_DICT)))
    print("Uncovered regions:")
    GENOME_SEQ = next(iter(REFERENCE_DICT.values()))[0]
    count_total = [0, 0]
    for e in UNALIGNED_POSITION:
        print("{}:{} {:s}".format(e[0], e[1], GENOME_SEQ[e[0]:e[1]]))
        count_total[0] += 1
        count_total[1] += (e[1] - e[0])
    print("Number of uncovered regions: {}\nNumber of uncovered bases: {}"
          .format(*count_total))

    # Testing purposes...
    # for _val in QUERY_DICT.values():
    #     print(_val[1])
    # print("\n\n")
    # for _val in REFERENCE_DICT.values():
    #     print(_val[1])
    # first = next(iter(QUERY_DICT))
    # print(first, QUERY_DICT[first][1])
    # first = next(iter(REFERENCE_DICT))
    # print(first, REFERENCE_DICT[first][1])
