#!/usr/bin/env python3
"""
Parsing Genbank file (.gb) and output FASTA and a statistics table

Author: Hidde Bleeker
Student no: 931202071020

"""

# import statements
from sys import argv
import re


def divide_into_entries(lines):
    """
    Reading the content of a Genbank file into separate entries

    :param lines: a file (e.g. open("filename.txt"))
    :return: Multiple lines (list of lines), spanning an entire entry (up until "//")
    """
    match = re.compile(r"^//")
    curr_entry = []
    for line in lines:
        if not line.split():
            continue
        elif match.search(line):
            yield curr_entry
            curr_entry = []
        else:
            curr_entry.append(line.strip())

    if curr_entry:
        yield curr_entry


class GenBankEntry(object):
    def __init__(self, gb_entry):
        self.access_n = None
        self.organism = None
        self.sequence = ""
        self.__parse_entry(gb_entry)

    def __parse_entry(self, gb_entry):
        """
        Parse one entry in GenBank format

        Only reading accession no, organism and sequence
        """
        for i in range(len(gb_entry)):
            if re.match(r"ACCESSION", gb_entry[i]):
                self.access_n = gb_entry[i].split("   ")[1]
            if re.match(r" {2}ORGANISM", gb_entry[i]):
                self.organism = gb_entry[i].split("  ")[1]
            if re.match(r"ORIGIN", gb_entry[i]):
                for line in gb_entry[i+1:]:
                    for base in line:
                        if base.upper() in "ACTG":
                            self.sequence += base.upper()

    def __calc_gc_content(self):
        """
        Function for calculating percentage of GC in a DNA sequence (type str)
        Not taken into account the occurrences of U, N, Y, K or other FASTA acid codes, or gaps.
        """
        n_gc = self.sequence.count("G") + self.sequence.count("C")
        percentage_gc = n_gc / len(self.sequence)
        if (percentage_gc - n_gc / len(self.sequence)) < 0.001:
            return 100 * percentage_gc

    ###  Get functions
    def get_access_no(self):
        return self.access_n

    def get_organism(self):
        return self.organism

    def get_seq(self):
        return self.sequence

    def get_gc_content(self):
        return self.__calc_gc_content()


if __name__ == '__main__':
    if len(argv) > 1:
        with open(argv[1]) as file:
            # gb_obj_map = []
            entries = []
            for entry in divide_into_entries(file):
                entries += [entry]
            gb_obj_list = [GenBankEntry(e) for e in entries]



