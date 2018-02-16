#!/usr/bin/env python3
"""
Usage: {} reference.fa forward.fq reverse.fq


Author: Hidde Bleeker
931202071020
"""

# Imports
from sys import argv
import subprocess as sbp
import os
import re



def print_help():
    """ Print the script description and usage
    :return: None. is a print function.
    """
    print(__doc__.format(argv[0]))


def run_bwa(ref_genome, fw_read, rev_read, out_fn="mapped.sam"):
    """ Run BWA; make index of reference FASTA and map forward and rev. reads
    Only supports mapping paired end reads for now

    :param ref_genome: The fasta file to use as a reference for mapping
    :param fw_read: Forward read
    :param rev_read: Reverse read
    :param out_fn: The output filename; when existing, bwa command not run
    :return: The name of the output sam file
    """
    if not any(os.path.exists(inp_file) for inp_file in
               [ref_genome, fw_read, rev_read]):
        raise FileExistsError("On of the input files does not exist")
    if not os.path.exists(out_fn):
        ind_cmd = "bwa index {}".format(ref_genome)
        is_err = sbp.check_call(ind_cmd, shell=True)
        if not is_err:
            print("Successfully created reference genome index with BWA")
        map_cmd = "bwa mem -t 6 -k 21 {0} {1} {2} > {3}"\
            .format(ref_genome, fw_read, rev_read, out_fn)
        is_err = sbp.check_call(map_cmd, shell=True)
        if not is_err:
            print("Successfully mapped input files to reference genome")
    else:
        print("Already an output with mapped reads exists: {}"
              .format(out_fn))
    return out_fn


def parse_sam(lines):
    """ Parse SAM-files, excluding entries with FLAG >= 2048

    :param lines: List of lines in the gff format
    :return: dictionary with QNAME as keys and list of tuples of FLAG, MAPQ,
    sequence and quality.
    """
    # name_flag_mapq = re.compile(r'(\S+)\s(\d+)\s\S+\s\d+\s(\d+)')
    # seq_qual = re.compile(r'([\S]+)\s([\S]+)\sNM:')
    match = re.compile(r'(\S+)\s(\d+)\s(\S+)\s\d+\s(\d+)(\s\S+)'
                       r'{4}\s+([\S]+)\s([\S]+)\s')
    reads_dic = {}
    for line in lines:
        if line.startswith("@") or not line.strip():
            continue
        else:
            [qname, flag, mapq, rname, sequence, quality] = \
                [match.search(line).group(i) for i in [1, 2, 4, 3, 6, 7]]
            if int(flag) >= 2048:
                continue
            if qname not in reads_dic:
                reads_dic[qname] = [(int(flag), int(mapq), rname, sequence,
                                     quality)]
            else:
                reads_dic[qname].append((int(flag), int(mapq), rname, sequence,
                                         quality))
    return reads_dic


def filter_mappings(mappings):
    """ Filter the mapped reads based on their flags.
    Print statistics to screen.
    Read nuclear mappings to file.

    :param mappings: Dictionary of mapped reads
    :return: Filtered dictionary of remaining (nuclear) mappings and dic of
    statistics of the organellar mappings
    """
    remaining = mappings.copy()
    organellar = {}
    stats = {'chloro': [],
             'mito': [],
             'total': []}
    for key, value in mappings.items():
        if (any(a[0] == 99 for a in value) and any(a[0] == 147 for a in value))\
                or (any(a[0] == 83 for a in value) and any(a[0] == 163
                                                           for a in value)):
            organellar[key] = [value]
            remaining.pop(key)
    for key, value in organellar.items():
        if 'chloro' in any(a[2] for a in value):
            stats['chloro'].append(a[1] for a in value)
        if 'mito' in any(a[2] for a in value):
            stats['mito'].append(a[1] for a in value)
        stats['total'].append(a[1] for a in value)
    return remaining, stats


def output_reads_stats(remaining, stats):
    """ Print out remaining stats (and forward/reverse fq files (no time))

    :param remaining: dic with remaining mappings
    :param stats: dict of stats chloro/mito mappings
    :return: a file
    """
    print("Group\t#proper\tavg MAPQ")
    for key, val in stats.items():
        print("{}\t{}\t{}".format(key, len(val), sum(val) / len(val)))


if __name__ == '__main__':
    # Print help if not used correctly/wrong amount of arguments
    if not len(argv) == 4:
        print_help()
    else:
        # Check if the input files are in the right format
        if any(argv[1].lower().endswith(ext) for ext in ['.fa', ',fasta']) \
            and any(argv[2].lower().endswith(ext) for ext in ['.fq', '.fastq'])\
            and any(argv[3].lower().endswith(ext) for ext in ['.fq', '.fastq']):

            # Read arguments and use them to run bwa as defined in function
            reference, forward_read, reverse_read = argv[1], argv[2], argv[3]
            mapped = run_bwa(reference, forward_read, reverse_read)

            # Open the dictionary with mapped reads and parse
            with open(str(mapped), 'r') as alignments:
                mapped_dic = parse_sam(alignments)
                nuclear, statistics = filter_mappings(mapped_dic)
                output_reads_stats(nuclear, statistics)

