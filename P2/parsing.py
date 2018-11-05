#!/usr/bin/env python3
"""Parse a GenBank formatted input file, output FASTA and tab-delim. stats files

Author: Hidde Bleeker
931202071020

usage: {} input.gb output_filename
"""

# Imports
from sys import argv
from operator import itemgetter


def parse_gb_input(gb_filename):
    """Parse a GenBank formatted file and return dictionary with info per record

    :param gb_filename: str, input filename
    :return: generator, with elements as accession_nr, organism_name, sequence
    """
    accession, org_name, seq = [], "", []
    reading_seq = False
    with open(gb_filename) as gb_ip:
        for line in gb_ip:
            if line.startswith('//'):
                reading_seq = False
                yield accession, org_name, "".join(seq)
                accession, org_name, seq = [], "", []
            elif reading_seq:
                for partial_seq in line.strip().upper().split()[1:]:
                    seq.append(partial_seq)
                continue
            elif line.startswith('ACCESSION'):
                # Collect (only!) primary accession number
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


def write_sorted_fasta(seqs_sort, output_fn):
    """Print fasta file from given sequence dictionary

    :param seqs_sort: list, with tuples as (GB accession numbers,
    {org_name, seq, GC_percentage}); sorted on order of printing
    :param output_fn: Name of the output file to print the sequences to, without
     extension
    :return: None, it is a function that writes an output file
    """
    with open(output_fn + '.fasta', 'w') as out_f:
        for entry in seqs_sort:
            out_f.write(">{acc:s} {org_name:s}\n{seq:s}\n".format(**entry))


def write_sorted_stats(seqs_sort, output_fn):
    """Print statistics file from given sequence dictionary

    :param seqs_sort: list, with tuples as (GB accession numbers,
    {org_name, seq, GC_percentage}); sorted on order of printing
    :param output_fn: Name of the output file to print the sequences to, without
     extension
    :return: None, it is a function that writes an output file
    """
    with open(output_fn + '.csv', 'w') as out_f:
        for entry in seqs_sort:
            out_f.write("{acc:s}\t{org_name:s}\t{GC_content:.2f}\t{length:d}\n"
                        .format(**entry))


if __name__ == '__main__':
    # Check for right amount of input arguments
    if len(argv) != 3:
        print("Not the right amount of command line arguments given...\n"
              "Quitting.")
        print(__doc__.format(argv[0]))
        exit(1)
    # Check if second input argument has the correct file extension
    elif not argv[1].lower().endswith('.gb'):
        print("File in incorrect input format given (should be .gb)...\n"
              "Quitting.")
        print(__doc__.format(argv[0]))
        exit(1)
    else:
        # Create dictionary from parsed GB sequences and write FASTA/stats to
        # files.
        sequences = [{'acc': acc, 'org_name': name, 'seq': seq, 'GC_content': calc_gc_cont(seq), 'length': len(seq)} for acc, name, seq in parse_gb_input(argv[1])]

        sorted_sequences = sorted(sequences, key=itemgetter('GC_content'),
                                  reverse=True)
        write_sorted_fasta(sorted_sequences, argv[2])
        write_sorted_stats(sorted_sequences, argv[2])
