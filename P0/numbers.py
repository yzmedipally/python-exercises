#!/usr/bin/env python3
"""
Read numbers from command line and print sum of the numbers

Edit for git purposes

Author: Hidde Bleeker
931202071020
"""

# Imports
from sys import argv


if __name__ == '__main__':
    if len(argv) > 1:
        print("{} = {}".format("+".join(argv[1:]),
                               sum(list(map(int, argv[1:])))))
