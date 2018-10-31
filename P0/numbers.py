#!/usr/bin/env python3
"""
Read numbers from command line and print sum of the numbers

Author: Hidde Bleeker
931202071020
"""

# Imports
from sys import argv


def sum_and_print(numbers):
    """Sum and print a list of numbers"""
    print(sum(numbers))


if __name__ == '__main__':
    if len(argv) > 1:
        sum_and_print(list(map(int, argv[1:])))
