#!/usr/bin/env python3
"""
Usage: {} predicted.gff reference.fa



Author: Hidde Bleeker
931202071020
"""

# Imports
from sys import argv


def print_help():
    """ Print the script description and usage
    :return: None. is a print function.
    """
    print(__doc__.format(argv[0]))


#  rest goes here


if __name__ == '__main__':
    if not len(argv) == 3:
        print_help()
    else:
        if argv[1].endswith('.gff') and argv[2].endswith(".fa"):
            predicted, reference = argv[1], argv[2]
