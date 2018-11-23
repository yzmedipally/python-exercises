#!/usr/bin/env python3
"""Script to

Author: Hidde Bleeker (931202071020)
"""

# Imports
import argparse
import subprocess as sbp
import re
from sys import argv
from os.path import exists


def parse_arguments():
    """Create an argument parser and parse arguments from argv

    :return: Parsed arguments; as {argument_name: argument_value}
    :rtype: dict of {str: str}
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('assembly', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.fa', '.fasta'])
                        else parser.error("Wrong input file extension."),
                        help="Input genome assembly FASTA file, with .fasta or"
                             " .fa file extension.")
    parser.add_argument('-o', '--augustus_out', default='aug.gff',
                        type=lambda inp: inp if inp.lower().endswith('.gff')
                        else parser.error("Wrong input file extension."),
                        help="Augustus output filename, with .gff file "
                             "extension.")
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


def run_augustus(in_fn, genemodel='complete', gff3='on', out_fn='aug.gff',
                 species='aspergillus_nidulans'):
    """Run augustus with parameters supplied if output file does not exist

    :param str in_fn: Input (query) filename
    :param str genemodel: Augustus genemodel parameter, default is 'complete'
    :param str gff3: Augustus gff3 parameter, default is 'on'
    :param str out_fn: Output filename, default is 'aug.gff'
    :param str species: Augustus species parameter, default is
    'aspergillus_nidulans'
    :return: None, calls command line utility
    """
    if exists(out_fn):
        print("## Output file already exists; augustus was not run.")
        return
    cmd = "augustus --genemodel={} --gff3={} --outfile={} --species={} {}"\
        .format(genemodel, gff3, out_fn, species, in_fn)

    print("## Running augustus, command used:\n## {}".format(cmd))
    try:
        sbp.check_call(cmd, shell=True)
    except sbp.CalledProcessError as err:
        print("## Augustus failed to run. Error:\n## {}\n## Quitting."
              .format(err))
    # implicit return


def parse_gff3(list_of_lines):
    """Parse a gff3 file for predicted genes

    :param iterator list_of_lines: Input file to read, e.g. an opened gff3 file
    :return: generator of predicted genes, as (id, (start, end), strand)
    :rtype generator of (str, (int, int), str)
    """
    curr = None
    index = None
    strand = ""
    pattern = re.compile(r'^\S+\s+\S+\s+gene\s+(\d+)\s+(\d+)\s+\S+\s+([+-])')
    for line in list_of_lines:
        if not line.strip():
            continue
        if line.startswith('#'):
            start = re.search(r'^# start gene (\S+)', line)
            end = re.search(r'^# end gene (\S+)', line)
            if start:
                curr = start.group(1)
            if end:
                if curr:
                    yield curr, index, strand
            continue
        result = pattern.search(line)
        if result:
            index = tuple(map(int, [result.group(i) for i in [1, 2]]))
            strand = result.group(3)



if __name__ == '__main__':
    # # In the case that argparse is not allowed:
    # if not len(argv) == 3:
    #     print("usage: python3 {} reference.fa query.fa".format(argv[0]))
    # # In this case every ARGS[key] has thus to be replaced with argv[1] or
    # # argv[2].

    # Read command line arguments into dictionary:
    ARGS = parse_arguments()
    # print(type(ARGS), ARGS)

    # Run augustus based on command line parameters.
    run_augustus(ARGS['assembly'], out_fn=ARGS['augustus_out'])

    # Open and read assembly FASTA and Augustus output file:
    with open(ARGS['assembly']) as inp_f:
        ASSEMBLY = {_id: (_seq, len(_seq)) for _id, _seq in
                    parse_fasta(inp_f)}
    with open(ARGS['augustus_out']) as gff_inp:
        PREDICTED = {_id: (_index, _strand) for _id, _index, _strand in
                     parse_gff3(gff_inp)}
        PREDICTED_ORDER = [_id for _id, _index, _strand in
                           parse_gff3(gff_inp)]
    # # Testing:
    # for entry in PREDICTED.items():
    #     print(entry)
    print(PREDICTED_ORDER)

