#!/usr/bin/env python3
"""
Author: Hidde Bleeker
"""

# Imports:
import re


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


def parse_gff(lines):
    """ Parse a GFF file for protein sequence.

    :param lines: List of lines (or a file) in the GFF format.
    :return: a dictionary with {gene id: protein sequence}
    """
    begin = re.compile(r"# start gene ([\w]+)")
    end = re.compile(r"# end gene ([\w]+)")
    prot_seq = re.compile(r"protein sequence = \[([\w]+)\]")
    predicted = {}
    reading = False
    curr = ""
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("# "):
            if reading:
                if end.search(line):
                    predicted[end.search(line).group(1)] = \
                        prot_seq.search(curr).group(1)
                    curr = ""
                    reading = False
                else:
                    curr += line.strip()[2:]
            elif begin.search(line):
                reading = True
    return predicted

