#!/usr/bin/env python3
"""

Author: Hidde Bleeker (931202071020)

usage: python3 {} reference.fa related.fa
"""

# Imports
from sys import argv, exit
from os.path import exists
import subprocess as sbp


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


def run_needle(ref_seq_fn, rel_seq_fn, out_fn="out.needle", gap_open=8,
               gap_ext=0.5):
    """Run EMBOSS needle on the command line with given given parameters

    :param ref_seq_fn: str, single sequence (reference) input filename
    :param rel_seq_fn: str, input filename of sequences to compare
    :param out_fn: str, output filename of aligned sequences
    :param gap_open: float, penalty for opening a gap in the alignment
    :param gap_ext: float, penalty for extending a gap in the alignment
    :return: None, runs a program.
    """
    if not 3 <= gap_open <= 10:
        raise ValueError("The gap open penalty should be between 3 and 10.")
    if exists(out_fn):
        print("Output filename already exists. Needle is not run.")
    else:
        cmd = "needle -asequence {} -bsequence {} -gapopen {} -gapextend {} " \
              "-outfile {}".format(ref_seq_fn, rel_seq_fn, gap_open, gap_ext,
                                   out_fn)
        try:
            e = sbp.check_call(cmd, shell=True)
        except sbp.CalledProcessError as e:
            print("Needle command failed to run. Error:\n{}".format(e))


def calc_hamm_dist(seq1, seq2):
    """Calculate Hamming distance between two sequences of equal length

    :param seq1: str, sequence of same length as seq2
    :param seq2: str, sequence of same length as seq1
    :return: integer or None, Hamming distance (number of positions at which
    the symbols are different)(None if seqs not of same length)
    """
    if len(seq1) != len(seq2):
        print("Sequences not of same length. Returning without answer.")
        return
    return sum([0 if seq1[i] == seq2[i] else 1 for i in range(len(seq1))])


def calc_perc_id(seq1, seq2):
    """Calculate percentage identity between two aligned sequences

    :param seq1: str, sequence of same length as seq2
    :param seq2: str, sequence of same length as seq1
    :return: float or None, percentage identity (no_identical_pos /
    alignment_length * 100) (None if seqs not of same length)
    """
    if len(seq1) != len(seq2):
        print("Sequences not of same length. Returning without answer.")
        return
    return (len(seq1) - calc_hamm_dist(seq1, seq2)) / len(seq1) * 100


if __name__ == '__main__':
    if not len(argv) == 3:
        print(__doc__.format(argv[0]))
        exit(1)

    # Step 1: Read filenames specified on the command line (using argv)
    ref_fn, rel_fn = argv[1:]
    # ref_fn, rel_fn = "ref.fasta", "related.fasta"

    # 2. Parse FASTA files with multiple protein sequences and determine the
    # lengths of the sequences
    with open(rel_fn) as rel_f:
        rel_seqs = {_id: (_seq, len(_seq)) for _id, _seq in parse_fasta(rel_f)}
    for k, v in rel_seqs.items():
        print("{}:\t{}\t{}".format(k, v[1], v[0]))

    # 3. In your python script, run the program needle to align the protein
    # sequences of related species (related.fasta) to the reference protein
    # from A. thaliana (ref.fasta).
    run_needle(ref_fn, rel_fn)

    # 4. Write a function to calculate the hamming distance between two
    # sequences of equal length
    # 5. Write a function to calculate the percent of identity between two
    # aligned sequences

    # 6. Parse the needle output to extract all pairwise alignments

    # 7. Report a tab-delimited table with one line for each alignment,
