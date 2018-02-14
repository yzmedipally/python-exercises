#!/usr/bin/env python3
"""
Author: Hidde Bleeker
"""

# Imports:
import os
import subprocess as sbp


# Script from pfam_domains.py  !!! RETURNS LOG NOT OUTPUT FILENAME !!!
def hmmpress_db(location, database_name):
    """ Use HMMpress on a hmm database file to construct binary HMM files used in the domain database lookup

    :param location: Location/directory (relative to the script run from, or full)
    :param database_name: Name of the database to run HMMpress on
    :return: Error message or 0 (zero) if succeeded
    """
    if os.path.exists(location + "/" + database_name):
        print("Run HMMPress on database")
        try:
            hmmpress_run = sbp.check_call("hmmpress " + database_name, cwd=location, shell=True)
        except sbp.CalledProcessError as err:
            err_msg = "\n### ERROR ###\n" \
                      "HMMpress did not succeed. \nCommand: {}\nReturn code: {}\n"\
                .format(err.cmd, err.returncode)
            print(err_msg)
            return err_msg
        return hmmpress_run
    return "\n### ERROR ###\nPath {} not found.".format(location)


# Example with BLAST db
def make_blast_db(input_fasta, output_db_file="blast_prot_db"):
    """ Create BLAST database from input reference genome

    :param input_fasta: The fasta file to use as a database (string)
    :param output_db_file: The name of the output database
    :return: The name of the output database
    """
    if not os.path.exists(input_fasta):
        raise FileExistsError("Input FASTA file for creating database ('{}') does not exist".format(input_fasta))
    if not os.path.exists(output_db_file):
        cmd = "makeblastdb -in {} -input_type fasta -dbtype prot -out {}".format(input_fasta, output_db_file)
        is_err = sbp.check_call(cmd, shell=True)
        if not is_err:
            print("Successfully created BLAST database")
    else:
        print("Output database already exists.")
    return output_db_file


# Example with running BLAST
def run_blast(query, database, out_file="yeast.blast", option_dic=None,):
    """ Run BLAST with query input to BLAST database

    :param query: Query filename in (FASTA format)
    :param database: The BLAST database
    :param out_file: filename of the BLAST output
    :param option_dic: Dictionary with additional BLAST options
    :return: filename of the BLAST output
    """
    if not all(os.path.exists(fn) for fn in [query, database]):
        raise FileExistsError("Input FASTA file and/or database ('{}' and/or '{}') do not exist"
                              .format(query, database))
    options = ""
    if option_dic:
        print("BLAST options:")
        for key, val in option_dic.items():
            options += " -" + str(key) + " " + str(val)
            print("{}\t{}".format(key, val))
    else:
        print("BLAST run with default options.")

    if not os.path.exists(out_file):
        cmd = "blastp -query {} -db {} -out {}{}"\
            .format(query, database, out_file, options)
        is_err = sbp.check_call(cmd)
        if not is_err:
            print("Successfully run BLAST against database.")
    else:
        print("BLAST output file already exists; BLAST search is not run again.")
    return out_file
