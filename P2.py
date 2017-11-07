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
        self.gc_content = self.__calc_gc_content()

    def __parse_entry(self, gb_entry):
        """
        Parse one entry in GenBank format

        Only reading accession no, organism and sequence
        """
        for i in range(len(gb_entry)):
            if re.match(r"ACCESSION", gb_entry[i]):
                self.access_n = gb_entry[i].split("   ")[1]
            if re.search(r"ORGANISM", gb_entry[i]):
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

    def output_fasta_string(self):
        assert self.access_n is not None and self.sequence != ""
        return ">{0:s}\n{1}\n".format(self.access_n, self.sequence)

    def output_tab_delim(self):
        assert [e for e in [self.access_n, self.organism, self.gc_content, self.sequence]] is not None
        return "{0:15s}\t{1:30s}\t{2:.2f}%\t{3:5d}\n".format(self.access_n, self.organism, self.gc_content, len(self.sequence))

    #  Get functions
    def get_gc_content(self):
        if self.gc_content is None:
            self.gc_content = self.__calc_gc_content()
        return self.gc_content


if __name__ == '__main__':
    if len(argv) > 1:
        filename_no_ext = re.search(r'.+\.', argv[1])
        with open(argv[1]) as file:
            # gb_obj_map = []
            entries = []
            for entry in divide_into_entries(file):
                entries += [entry]
            gb_obj_list = [GenBankEntry(e) for e in entries]
            gb_obj_sorted = sorted(gb_obj_list, key=lambda gb_obj: gb_obj.gc_content, reverse=True)

            with open(filename_no_ext.group() + "fasta", 'w') as fasta_dest:
                for entry in gb_obj_sorted:
                    fasta_dest.write(entry.output_fasta_string())

            with open(filename_no_ext.group() + "txt", 'w') as text_dest:
                for entry in gb_obj_sorted:
                    text_dest.write(entry.output_tab_delim())
