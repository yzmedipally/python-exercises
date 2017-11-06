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
    def __init__(self):
        self.access_n = None
        self.organism = None
        self.sequence = ""

    def parse_entry(self, gb_entry):
        """
        Parse one entry in GenBank format

        Only reading accession no, organism and sequence
        """
        for field in range(len(gb_entry)):
            if re.search(r"^ACCESSION", gb_entry[field]):
                self.access_n = gb_entry[field].split("   ")[1]
            if re.search(r"^  ORGANISM", gb_entry[field]):
                self.organism = gb_entry[field].split("  ")[1]
            if re.search(r"^ORIGIN", gb_entry[field]):
                for line in gb_entry[field+1:]:
                    self.sequence.join([c for c in line if c not in "1234567890 "])









if __name__ == '__main__':
    if len(argv) > 1:
        with open(argv[1]) as file:
            for entry in divide_into_entries(file):
                pass



