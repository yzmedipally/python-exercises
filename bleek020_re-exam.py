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


if __name__ == '__main__':
    if not len(argv) == 4:
        print_help()
    else:
        if any(argv[1].lower().endswith(ext) for ext in ['.fa', ',fasta']) \
            and any(argv[2].lower().endswith(ext) for ext in ['.fq', '.fastq'])\
            and any(argv[3].lower().endswith(ext) for ext in ['.fq', '.fastq']):
            reference, forward_read, reverse_read = argv[1], argv[2], argv[3]
            print(run_bwa(reference, forward_read, reverse_read))
