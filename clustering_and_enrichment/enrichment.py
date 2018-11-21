#!/usr/bin/env python3
"""

Author: Hidde Bleeker (931202071020)
"""

# Imports
import argparse
from sys import argv
from scipy.stats import fisher_exact


def parse_arguments():
    """Create an argument parser and parse arguments from argv

    :return: Parsed arguments.
    :rtype: dict, with {parameter_name: parameter_value}
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gene_expression', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.csv', '.tsv'])
                        else parser.error("Wrong input file extension."),
                        help="Gene expression data, semicolon separated file, "
                             "with .csv or .tsv file extension.")
    parser.add_argument('ec2kegg', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.csv', '.tsv'])
                        else parser.error("Wrong input file extension."),
                        help="EC2KEGG index, semicolon separated file "
                             "with .csv or .tsv file extension.")
    parser.add_argument('pathways', type=lambda inp: inp
                        if any(inp.lower().endswith(tail) for tail in
                               ['.csv', '.tsv'])
                        else parser.error("Wrong input file extension."),
                        help="KEGG pathways, semicolon separated file "
                             "with .csv or .tsv file extension.")
    if len(argv) <= 1:
        parser.print_help()
        exit(1)
    else:
        return vars(parser.parse_args())


def float_or_int_or_string(entry):
    """Convert a string to float or int if that is possible.

    :param str entry: String to convert to float or int
    :return: converted entry string
    :rtype float, int, str depending on if conversion is possible.
    """
    entry = entry.strip()
    try:
        if '.' in entry:
            return float(entry)
        else:
            return int(entry)
    except ValueError:
        if entry.lower() == 'no':
            return False
        elif entry.lower() == 'yes':
            return True
        else:
            return str(entry)


def parse_comma_delim(list_of_lines):
    """Read in a comma-delimited file with vertical and horizontal headers

    :param generator list_of_lines: generator of lines, tab-delimited with
    floats as data
    :return: data matrix as dictionary in dictionary
    :rtype dict of dicts of {gene_names (str): {descriptives (str): values
    (int, float or str)}
    """
    description = next(list_of_lines).strip().split(';')[1:]
    data_dic = {}
    while True:
        try:
            line = next(list_of_lines).strip().split(';')
            data_dic[line[0].strip()] = {_k.strip(): float_or_int_or_string(_v)
                                         for _k, _v in zip(description,
                                                           line[1:])}
        except StopIteration:
            return data_dic


def parse_ec2kegg_index(list_of_lines):
    """Read in a file describing KEGG pathways per EC reaction number

    :param generator list_of_lines: generator of lines, tab-delimited with
    floats as data
    :return: data matrix as dictionary in dictionary
    :rtype dict of {pathway description (str): {ec_numbers} (set of str)}
    """
    first_time = True
    kegg2ec = {}
    for line in list_of_lines:
        if first_time:
            first_time = False
            continue
        splitted = line.replace('"', '').strip().split(";")
        for pw in splitted[1:]:
            if pw in kegg2ec:
                kegg2ec[pw].append(splitted[0])
            else:
                kegg2ec[pw] = [splitted[0], ]
    return {_key: set(_val) for _key, _val in kegg2ec.items()}


def ec_in_kegg_pw(kegg2ec_index, pathways_index):
    """Build dict describing EC numbers belonging to KEGG pathway from indexes

    :param dict kegg2ec_index: {ec_number (str): pathway description (str)}
    :param dict pathways_index: {kegg_pathway_id (str): pathway name (str)}
    :return: dictionary describing which EC numbers belong to a KEGG pathway
    :rtype dict, of {kegg_pathway_id: ec_number}
    """
    return {_id: kegg2ec_index[_pw] for _id, _pw in pathways_index.items()}



if __name__ == '__main__':

    args = parse_arguments()
    # print(type(args), args)
