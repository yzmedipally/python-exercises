#!/usr/bin/env python3
"""Script to Retrieve sequences from GenBank, print stats and shortest in FASTA

Author: Hidde Bleeker (931202071020)

usage: python3 {} GenBank_accession_numbers.txt
"""

# Imports
from sys import argv
from Bio import Entrez


# Rewritten GenBank parse function, as the one from P2 opened the file as well
# (not modular enough...)
def parse_gb_entry(gb_ip):
    """Parse a single GenBank entry and return dictionary with info per record

    :param gb_ip: generator, of GenBank entry
    :return: generator, (accession, org_name, sequence) per entry
    """
    accession, org_name, seq = [], "", []
    reading_seq = False
    for line in gb_ip:
        if line.startswith('//'):
            reading_seq = False
            return accession, org_name, "".join(seq)
            # accession, org_name, seq = [], "", []
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


def parse_acc_nums(ip_filename):
    """Parse a text file with GenBank accession numbers and return as a list
    :param ip_filename: str, input filename to parse
    :return: list of str, GenBank accession numbers
    """
    entries = []
    with open(ip_filename, 'r') as ip_f:
        for line in ip_f:
            entries += line.strip().split()
    return entries


def retrieve_gb_seqs(accession_list):
    """Retrieve the GenBank sequences from a list of accession numbers

    :param accession_list: list of str, GenBank accession numbers
    :return: generator, of Entrez.efetch objects (that is,
    <_io.TextIOWrapper encoding='latin-1'> )
    """
    Entrez.email = "hidde.bleeker@wur.nl"
    for _e in accession_list:
        yield Entrez.efetch(db='nucleotide', id=_e, rettype='gb')


if __name__ == '__main__':
    if not len(argv) == 2:
        print(__doc__.format(argv[0]))

    # Step 1, 2 & 3: Fetch seqs from the file given as argument; put in dict of
    #  dicts with info about org. name, sequence and sequence length per entry
    gb_seqs = retrieve_gb_seqs(parse_acc_nums(argv[1]))
    seqs = {acc_no: {'org_name': name, 'seq': seq, 'seq_len': len(seq)}
            for acc_no, name, seq in [parse_gb_entry(gb_seq)
                                      for gb_seq in gb_seqs]}

    # Step A: Print acc. no., organism name and sequence length for every entry
    for acc, vals in seqs.items():
        print("{}\t{org_name}\t{seq_len}".format(acc, **vals))

    # Step B: Print the label and sequence of shortest seq. in FASTA format
    shortest = min(seqs.items(), key=lambda x: x[1]['seq_len'])
    print(">{:s}\n{:s}".format(shortest[0], shortest[1]['seq']))
