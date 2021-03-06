#!/usr/bin/env python3
"""Run augustus and print stats of the predicted proteins and their stop codon

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
                        help="Input genome assembly FASTA file, with .fasta "
                             "or .fa file extension.")
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

    :param iterator list_of_lines: Input file to read, e.g. an opened gff file
    :return: generator of predicted genes, as (id, (start, end), strand)
    :rtype generator of (str, (int, int), str)
    """
    curr = None
    index = None
    strand = ""
    seq_id = ""
    patt = re.compile(r'^(\S+)\s+\S+\s+gene\s+(\d+)\s+(\d+)\s+\S+\s+([+-])')
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
                    yield curr, index, strand, seq_id
            continue
        result = patt.search(line)
        if result:
            index = tuple(map(int, [result.group(i) for i in [2, 3]]))
            strand = result.group(4)
            seq_id = result.group(1)


def calc_gc_fraction(seq):
    """Calculates GC content as a fraction for any given sequence

    :param str seq: Sequence to calculate GC fraction for
    :returns float: GC fraction """
    return sum(1 if base in 'GC' else 0 for base in seq.upper()) / len(seq)


def find_stop_codon(seq, indices, strand):
    """Find the stop codon in a sequence based on the index in the sequence

    :param str seq: (DNA) sequence where the protein coding sequence is in
    :param tuple/list indices: Pair of indices where the first element is
    smaller than the second
    :param strand:
    :return:
    """
    if strand not in '+-':
        raise ValueError("Strand information is not in the correct format")
    if indices[0] >= indices[1] or len(indices) != 2:
        raise ValueError("Index information is not in the correct format")
    if not any(isinstance(indices, type_) for type_ in [list, tuple]):
        raise TypeError("Index information is not in the correct type. Should"
                        " be list or tuple.")
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    if strand == '-':
        return "".join(reversed([complement[base] for base in
                                 seq[indices[0]:indices[0]+3]]))
    elif strand == '+':
        return "".join(seq[indices[1]-3:indices[1]])


def calc_gene_and_stop_stats(assembly, predicted):
    """Calculate GC content and stop cod. per gene and GC-content per stopcod.

    :param dict assembly: dictionary containing the assembled genome, as
    {_fragment_id (str): sequence (str)}
    :param predicted: dictionary containing the predicted genes, as
    {gene_id (str): {seq_id (str), (index (int, int)), strand (str)}
    :return: two dictionaries, the first dict with GC content and stop codon
    per gene id; the second with GC content per stop codon.
    :rtype tuple of dicts; ({str: (float, str)}, {str: float})
    """
    gc_stop = {}
    stop_average = {}
    for _id, _val in predicted.items():
        sequence = assembly[_val['seq_id']][0]\
            [_val['index'][0]:_val['index'][1]]
        stop = find_stop_codon(assembly[_val['seq_id']][0], _val['index'],
                               _val['strand'])
        # print(len(sequence), _val['index'], _val['strand'])
        gc_stop[_id] = (calc_gc_fraction(sequence) * 100, stop)
        if stop in stop_average:
            stop_average[stop].append(sequence)
        else:
            stop_average[stop] = [sequence, ]
    stop_average = {_k: calc_gc_fraction("".join(_v))
                    for _k, _v in stop_average.items()}
    return gc_stop, stop_average


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
        ASSEMBLY = {_id.split(" ")[0]: (_seq, len(_seq)) for _id, _seq in
                    parse_fasta(inp_f)}
    with open(ARGS['augustus_out']) as gff_inp:
        PREDICTED, PREDICTED_ORDER = {}, []
        for _id, _index, _strand, _seq_id in parse_gff3(gff_inp):
            PREDICTED[_id] = {'index': (_index[0] - 1, _index[1]),
                              'strand': _strand,
                              'seq_id': _seq_id, }
            PREDICTED_ORDER.append(_id)
    # # Testing:
    # for seq_id in ASSEMBLY:
    #     print(seq_id)
    # for entry in PREDICTED.items():
    #     print(entry)
    # print(PREDICTED_ORDER)

    # Calculate predicted gene and stop codon statistics:
    GENE_STATS, STOP_STATS = calc_gene_and_stop_stats(ASSEMBLY, PREDICTED)

    # Print the output:
    for gene in PREDICTED_ORDER:
        print("{:s}\t{:.3f}\t{:s}".format(gene, *GENE_STATS[gene]))
    print("---")
    for _k, _v in STOP_STATS.items():
        print("{:s}\t{:.3f}".format(_k, _v))
