#!/usr/bin/env python3
"""
Author: Hidde Bleeker
"""

# Imports:


# No need for generators
def write_fasta(fasta_dic, output_file):
    """ Write FASTA strings dictionary to file

    :param fasta_dic: Dictionary with {id: sequence}
    :param output_file: filename of the output fasta file (with .fa extension)
    :return: None (because it writes a file).
    """
    with open(output_file, "w") as out:
        for key, val in fasta_dic.items():
            out.write("{}\n{}\n".format(key, val))





if __name__ == '__main__':
    pass
