#!/usr/bin/env python3
"""
Exercise on clustering, making use of modules

"""

# import statements
import numpy
import scipy
import scipy.cluster.hierarchy as sch
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
from sys import argv, exit


def parse_expression_tsv(lines):
    start = True
    gene_names = []
    samples = {}
    for line in lines:
        if start:
            sample_ids = [a.strip() for a in line.split("\t")[1:]]
            for name in sample_ids:
                samples[name] = []
            start = False
        else:
            gene_names.append(line.split("\t")[0].strip())
            for i, name in enumerate(sample_ids):
                samples[name] += [line.split("\t")[i+1].strip()]
    return samples, gene_names


def pandas_plot_by_sample(samples, sample_ids, cluster_method='complete', distance_metric='euclidean'):
    panda_data = pd.DataFrame(samples, index=sample_ids).transpose()
    distances = sch.distance.pdist(panda_data, metric=distance_metric)
    # print(distances.shape)
    clustering = sch.linkage(distances, method=cluster_method)
    tree = sch.dendrogram(clustering)
    plt.savefig('plot_dendogram_by_sample.pdf', format='PDF')


def pandas_plot_by_gene(samples, sample_ids, cluster_method='complete', distance_metric='correlation'):
    panda_data = pd.DataFrame(samples, index=sample_ids).transpose().transpose()
    distances = sch.distance.pdist(panda_data, metric=distance_metric)
    clustering = sch.linkage(distances, method=cluster_method)
    tree2 = sch.dendrogram(clustering, leaf_font_size=1,
                           color_threshold=4, labels=sample_ids)
    plt.savefig('plot_dendogram_by_gene.pdf', format='PDF')


def sns_plot_map(samples, sample_ids, cluster_method='complete', distance_metric='correlation'):
    panda_data = pd.DataFrame(samples, index=sample_ids).transpose().transpose()
    distances = sch.distance.pdist(panda_data, metric=distance_metric)
    clustering = sch.linkage(distances, method=cluster_method)
    g = sns.clustermap(clustering, cmap="mako", robust=True)
    g.savefig('plot_map_gene_vs_sample.pdf', format='PDF')


if __name__ == '__main__':
    if len(argv) > 1:
        distance, method = 'correlation', 'complete'
        to_run = 'dend_sample'
        if '-i' not in argv:
            raise ValueError("No inputfile supplied. Exiting...")
        for i, arg in enumerate(argv):
            if arg == '--plot-type':
                try:
                    to_run = argv[i + 1]
                except IndexError:
                    print("Not specified which kind of plot to make (dend_sample, dend_gene or sns_map, \
                    defaulting to dend_sample")
            if arg == '-i':
                try:
                    input_file = argv[i + 1]
                except IndexError:
                    print("No input file given after -i option, exiting")
                    exit()
            if arg == '-d':
                try:
                    distance = argv[i + 1]
                except IndexError:
                    print("No distance parameter given after -d option, defaulting to correlation")
            if arg == '-m':
                try:
                    method = argv[i + 1]
                except IndexError:
                    print("No method parameter given after -m option, defaulting to complete")
        if '--plot-type' not in argv:
            print("Not specified what type of plot to make. Defaulting to dend_sample.\n"
                  "(to specify plot: cmd --plot-type [dend_sample | dend_gene | sns_map]")

        with open(input_file) as file:
            sample_dic, sample_ids = parse_expression_tsv(file)
            functions = {'dend_sample': pandas_plot_by_gene,
                         'dend_gene': pandas_plot_by_sample,
                         'sns_map': sns_plot_map}
            functions[to_run](sample_dic, sample_ids, distance_metric=distance, cluster_method=method)
