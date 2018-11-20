#!/usr/bin/env python3
"""Script for the clustering exercise in BIF-30806

Author: Hidde Bleeker (931202071020)

"""

# Imports
from sys import argv
import argparse
import numpy
import scipy
import scipy.cluster.hierarchy as sch
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import time


def parse_tab_delim(list_of_lines):
    """Read in a tab-delimited file with vertical and horizontal headers

    :param list_of_lines: generator of lines, tab-delimited with floats as data
    :return: list of str, sample names; dictionary of
    {gene names (str): list of floats (values)}
    """
    sample_list, gene_dic = [], {}
    count = 0
    for line in list_of_lines:
        if not count:
            sample_list = line.strip().split("\t")[1:]
        else:
            reformat = line.strip().replace(',', '.').split("\t")
            gene_dic[reformat[0]] = reformat[1:]
        count += 1
    print("Read {} gene(s).".format(count))
    return sample_list, gene_dic


def parse_arguments():
    """Function to create an argument parser and parse arguments from cmd-line

    :return: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file',
                        help="Input file, tab delimited file with vertical and"
                             " horizontal headers.")
    parser.add_argument('-d', '--distance', default='correlation',
                        choices=['correlation', 'euclidean'],
                        help="Specify the hierarchical clustering distance "
                             "type.")
    parser.add_argument('-l', '--linkage', default='complete',
                        choices=['single', 'complete', 'average'],
                        help="Specify the hierarchical clustering "
                             "linkage type.")

    if len(argv) <= 1:
        parser.print_help()
        exit(1)
    else:
        return vars(parser.parse_args())


def run_pandas_clustering(pd_dataframe, linkage, metric):
    """Run pandas commands to be used in the exercise

    :param pd_dataframe: pandas dataframe with genes as columns and samples as
     rows
    :param linkage: str, indicating linkage type
    :param metric: str, indicating clustering distance metric
    :return: None, plotting function.
    """
    distances = sch.distance.pdist(pd_dataframe, metric='euclidean')
    print(distances.shape)
    clustering = sch.linkage(distances, method='complete')
    tree = sch.dendrogram(clustering)
    # plt.show()
    # print(pd_dataframe)
    df_tr = pd_dataframe.transpose()
    distances = sch.distance.pdist(df_tr, metric=metric)
    clustering = sch.linkage(distances, method=linkage)
    tree = sch.dendrogram(clustering, leaf_font_size=2, color_threshold=4,
                          labels=pd_dataframe.columns.tolist())
    plt.show()


def run_seaborn_clustermap(pd_dataframe, linkage, metric):
    """Run seaborn commands to be used in the exercise

    :param pd_dataframe: pandas dataframe with genes as columns and samples as
     rows
    :param linkage: str, indicating linkage type
    :param metric: str, indicating clustering distance metric
    :return: None, plotting function.
    """
    clustermap = sns.clustermap(pd_dataframe, method=linkage, metric=metric)
    plt.show()


if __name__ == '__main__':
    # if not len(argv) in [2, 4]:
    #     print(__doc__.format(argv[0]))

    args = parse_arguments()
    print(args['distance'])
    with open(args['input_file']) as ip_f:
        samples, genes = parse_tab_delim(ip_f)
    df = pd.DataFrame(genes, index=samples)
    run_pandas_clustering(df, linkage=args['linkage'], metric=args['distance'])
    df = df[df.columns].astype(float)
    # print(df.dropna())
    run_seaborn_clustermap(df.transpose(), linkage=args['linkage'],
                           metric=args['distance'])
