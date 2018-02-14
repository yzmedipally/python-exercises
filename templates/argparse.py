#!/usr/bin/env python3
"""


Author: Hidde Bleeker
"""


# Imports
from sys import argv
import argparse


#  rest goes here


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-databasedir', '-db', required=True, nargs=1, help=".",
                        type=lambda inp: inp if inp.endswith('.fasta') else parser.error("Wrong file extension"))
    parser.add_argument('-input_sequence', '-i', required=True, nargs='+', help=".")

    # Parser testing:
    # args = vars(parser.parse_args("-db asdf -i asdf.fasta asdffdsa.fasta".split(" ")))

    if len(argv) < 2:
        parser.print_help()
    else:
        return vars(parser.parse_args())


if __name__ == '__main__':
    arguments = parse_arguments()
    print(arguments)
