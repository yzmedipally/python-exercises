#!/usr/bin/env python3
"""
Author: Hidde Bleeker


"""

# imports
from sys import argv


def print_help():
    help_str = None

    print("{}".format(help_str))


def next_arg_to_var(argv_index, type_of=str):
    var = None
    try:
        if len(argv) > argv_index + 1:
            var = type_of(argv[argv_index + 1])
        else:
            print("No input given for command line argument '{}'".format(argv[argv_index]))
    except NameError:
        print("Not able to convert to given type. Couldn't convert")
    return var


if __name__ == '__main__':
    if len(argv) == 1:
        print_help()
    else:
        if "-h" in argv or "--help" in argv:
            print_help()
        # elif option in argv:

