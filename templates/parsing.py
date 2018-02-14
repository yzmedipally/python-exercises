#!/usr/bin/env python3
"""
Author: Hidde Bleeker
"""

# Imports:
import re


def parse_fasta(lines):
    """Parse a FASTA file (possibly with multiple sequences). Sequences are stored by ID in a dictionary.

    :param lines: list of lines (or a file) in FASTA format
    :return: dictionary of {label:dna_seq}
    """
    seqs = {}
    for line in lines:
        if not line.strip():
            continue
        if line.startswith('>'):
            has_comma = re.search(r">([\w ]+),[\w ]*", line)
            if has_comma:
                label = has_comma.group(1)
            else:
                label = line.strip()[1:]
            seqs[label] = ""
        else:
            if label in seqs:
                seqs[label] += line.strip()
    return seqs


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

