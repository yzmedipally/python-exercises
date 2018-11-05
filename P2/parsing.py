#!/usr/bin/env python3
"""Parse a GenBank formatted input file, output FASTA file and tab-delimited statistics file

Author: Hidde Bleeker
931202071020

usage: {} input.gb output_filename
"""

# Imports
from sys import argv


def parse_gb_input(gb_filename):
    """Parse a GenBank formatted file and return dictionary with information per record

    :param gb_filename: str, input filename
    :return: generator, with list entries as [accession_number, (organism_name, sequence)]
    """
    accession, org_name, seq = [], [], []
    reading_seq = False
    with open(gb_filename) as gb_ip:
        for line in gb_ip:
            # print(line.strip(), line.startswith("//"))
            if line.startswith('//'):
                # print("startswith //")
                reading_seq = False
                yield accession, (org_name, "".join(seq))
                accession, org_name, seq = [], [], []
            elif reading_seq:
                for partial_seq in line.strip().upper().split()[1:]:
                    seq.append(partial_seq)
                continue
            elif line.startswith('ACCESSION'):
                # Collect (only) primary accession number
                accession = line.strip().split()[1]
            elif line.startswith('  ORGANISM'):
                org_name = " ".join(line.strip().split()[1:])
            elif line.startswith('ORIGIN'):
                reading_seq = True


def calc_gc_cont(dna_seq):
    """Calculates the GC content of a DNA (or RNA) sequence

    :param dna_seq: str, DNA/RNA sequence
    :return: float, percentage of G+C in the sequence
    """
    count = [1 if _s.upper() in 'GC' else 0 for _s in dna_seq]
    return sum(count) / len(count) * 100


def write_sorted_fasta(seq_dict, output_fn):
    """Print fasta file from given sequence dictionary

    :param seq_dict: dict, with GB accession numbers as keys and (organism, seqeuences,
    GC_percentage) as values
    :param output_fn: Name of the output file to print the sequences to, without extension
    :return: None, it is a function that writes a file output
    """
    seqs_sorted = sorted(seq_dict.items(), key=lambda v: v[1][2], reverse=True)
    with open(output_fn + '.fasta', 'w') as out_f:
        for entry in seqs_sorted:
            out_f.write("".join(['>', str(entry[0]), ' ', str(entry[1][0]), '\n']))
            out_f.write(entry[1][1] + '\n')


def write_sorted_stats(seq_dict, output_fn):
    """Print statistics file from given sequence dictionary

    :param seq_dict: dict, with GB accession numbers as keys and (organism, seqeuences,
    GC_percentage) as values
    :param output_fn: Name of the output file to print the sequences to, without extension
    :return:
    """
    seqs_sorted = sorted(seq_dict.items(), key=lambda v: v[1][2], reverse=True)
    with open(output_fn + '.csv', 'w') as out_f:
        for entry in seqs_sorted:
            out_f.write("{0:s}\t{1:s}\t{2:.2f}\t{3:d}\n"
                        .format(entry[0], entry[1][0], entry[1][2], len(entry[1][1])))


if __name__ == '__main__':
    # Check for right amount of input arguments
    if len(argv) != 3:
        print("Not the right amount of command line arguments given...\nQuitting.")
        print(__doc__.format(argv[0]))
        exit(1)
    # Check if second input argument has the correct file extension
    elif not argv[1].lower().endswith('.gb'):
        print("File in incorrect input format given (should be .gb)...\nQuitting.")
        print(__doc__.format(argv[0]))
        exit(1)
    else:
        # Create dictionary from parsed GB sequences and write FASTA/stats to files.
        sequences = {acc: (name, seq, calc_gc_cont(seq)) for acc, (name, seq) in
                     parse_gb_input(argv[1])}
        write_sorted_fasta(sequences, argv[2])
        write_sorted_stats(sequences, argv[2])
