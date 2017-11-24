#!/usr/bin/env python3
"""
Author: Hidde Bleeker


"""

# imports
from sys import argv
import re


def print_help():
    """ Print a usage manual

    :return: None. is a print function.
    """
    help_str = None

    print("{}".format(help_str))


def next_arg_to_var(argv_index, type_of=str):
    """ Find command line arguments and if necessary convert them

    :param argv_index: index of argument specification (following item in argv is the actual parameter)
    :param type_of: type to convert the argument to
    :return: converted variable, on location argv_index + 1 (1 place further than given index)
    """
    var = None
    try:
        if len(argv) > argv_index + 1:
            var = type_of(argv[argv_index + 1])
        else:
            raise IndexError("No input given for command line argument '{}'".format(argv[argv_index]))
    except NameError:
        print("Not able to convert to given type. Couldn't convert")
    return var


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
    start_pattern = re.compile(r"# start gene ([\w]+)")
    prot_seq_start = re.compile(r"# protein sequence = \[([\w]+)")
    prot_seq_end = re.compile(r"# ([\w]+)\]")
    predicted = {}
    current_gene_id, current_gene_seq = "", ""
    reading_seq = False
    for line in lines:
        if not line.strip():
            continue
        elif line.startswith("#"):
            if reading_seq:
                ending = prot_seq_end.search(line)
                if not ending:
                    seq = re.search(r"# ([\w]+)", line)
                    if seq:
                        current_gene_seq += seq.group(1)
                else:
                    reading_seq = False
                    current_gene_seq += ending.group(1)
                    predicted[current_gene_id] = current_gene_seq
                    current_gene_id, current_gene_seq = "", ""
                continue
            start = start_pattern.search(line)
            seq_start = prot_seq_start.search(line)
            if start:
                current_gene_id = start.group(1)
            if seq_start:
                reading_seq = True
                current_gene_seq += seq_start.group(1)
            else:
                continue
    return predicted


def fasta_dic_to_string(fasta_dic):
    """ Generator function to prevent too large file being loaded into memory when writing fasta files

    :param fasta_dic: Dictionary with {id: sequence}
    :return: generator with formatted string in FASTA format for every entry in the dictionary
    """
    # output_str = ""
    for key, val in fasta_dic.items():
        yield ">{}\n{}\n".format(key, val)


def write_fasta(fasta_dic, output_file="predicted"):
    """ Write FASTA strings from generator function to file

    :param fasta_dic: Dictionary with {id: sequence}
    :param output_file: filename of the output fasta file (without .fa extension)
    :return: None (because it writes a file).
    """
    with open(output_file + ".fa", "w") as file:
        for entry in fasta_dic_to_string(fasta_dic):
            file.write(entry)


def run_blast(input_file1, input_file2, output_file="yeast.blast"):
    pass


if __name__ == '__main__':
    if len(argv) == 1:
        print_help()
    else:
        predicted_gff, reference_fa, output_filename = None, None, None
        try:
            if "-h" in argv or "--help" in argv:
                print_help()
            elif "-p" in argv:
                predicted_gff = next_arg_to_var(argv.index("-p"))
            elif "-r" in argv:
                reference_fa = next_arg_to_var(argv.index("-r"))
            elif "-o" in argv:
                output_filename = next_arg_to_var(argv.index("-o"))
        except IndexError:
            print_help()

        if predicted_gff is not None:
            predicted_fasta = parse_gff(predicted_gff)
            write_fasta(predicted_fasta)

        if all([reference_fa, predicted_gff]) is not None:
            run_blast(reference_fa, predicted_gff)








