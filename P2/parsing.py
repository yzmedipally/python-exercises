#!/usr/bin/env python3
"""Parse a Genebank formatted input file, output FASTA file and tab-delimited statistics file

Author: Hidde Bleeker
931202071020
"""

# Imports
from operator import itemgetter
from sys import argv


def parse_gb_input(gb_filename):
    """Parse a GeneBank formatted file and return dictionary with information per record

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


def write_fasta(seq_dict, output_fn):
    """Print fasta file from given sequence dictionary

    :param seq_dict: dict, with GB accession numbers as keys and (organism, seqeuences,
    GC_percentage) as values
    :param output_fn: Name of the FASTA output file to print the sequences to
    :return: None, it is a function that writes a file output
    """
    seqs_sorted = sorted(seq_dict.items(), key=lambda v: v[1][2], reverse=True)
    with open(output_fn, 'w') as out_f:
        for entry in seqs_sorted:
            out_f.write("".join(['>', str(entry[0]), ' ', str(entry[1][0])]))
            out_f.write(entry[1][1])


if __name__ == '__main__':
    if len(argv) != 2:
        print("Not the right amount of command line arguments given.\nQuitting.")
        exit(1)
    elif not argv[1].lower().endswith('.gb'):
        print("File in incorrect input format given (should be .gb).\nQuitting.")
        exit(1)
    else:
        sequences = {acc: (name, seq, calc_gc_cont(seq)) for acc, (name, seq) in
                     parse_gb_input(argv[1])}
        # Create list of sequences sorted by GC percentage
        # seqs_sorted = sorted(sequences.items(), key=lambda v: v[1][2], reverse=True)
        # print([(e[0], e[1][2]) for e in seqs_sorted])
