#!/usr/bin/env python3
"""Script to count the number of k-mers in a FASTA file

Output is compared to jellyfish output
Author: D.E. Bug; 
fixed by: Hidde Bleeker
931202071020
"""

import sys
import subprocess
import os.path


def parse_fasta(lines):
    """Return dictionary of {label:dna_seq}
    
    lines: list of lines in FASTA format
    """
    seqs = {}
    for line in lines:
        if not line.strip():
            continue
        if line.startswith('>'):
            label = line.strip()[1:]
            seqs[label] = ""
        else:
            seqs[label] += line.strip()
    return seqs


def extract_kmers(seqs, k=15, skip_unknown=True):
    """Return dict of {kmer_seq: kmer_count}

    seqs: dict of {label:dna_seq}
    k: int, specifying k-mer size
    skip_unknown: bool, specifying whether k-mers containing
        non-TGAC characters should be skipped
    """
    kmer_size = k
    ch = {}  # dict to store characters
    res = {}  # dict to store k-mers and counts
    for label, seq in seqs.items():
        for c in seq:
            if c not in ch:
                ch[c] = 0
            ch[c] += 1
        for i in range(len(seq) - kmer_size + 1):
            kmer = seq[i:i + kmer_size]
            count = True
            for c in kmer:
                if c not in "TGAC":
                    count = False
            if skip_unknown and not count:
                continue
            if kmer not in res:
                res[kmer] = 0
            res[kmer] += 1
    return res


def print_stats(kmer_table):
    """Print the kmer statistics to stdout
    
    kmer_table: dict of {kmer_seq: kmer_count}
    """
    print("MY OUTPUT")
    res = kmer_table
    unique = [i for i, j in res.items() if j == 1]
    print("Unique: {}".format(len(unique)))
    print("Distinct: {}".format(len(res)))
    total = sum(res.values())
    print("Total: {}".format(total))
    max_count = max(res.values())
    print("Max count: {}".format(max_count))
    for k, v in res.items():
        if v == max_count:
            print(k, v)
    print('----')
    return None


def run_jellyfish(input_fn, kmer_size=15):
    """Run jellyfish program on fasta file

    input_fn: string, filename of input FASTA file
    kmer_size: int, size of k-mers used by jellyfish
    """
    out_fn = 'tomato{}'.format(kmer_size)
    cmd = 'jellyfish count -m {} -s 1000000 -o {} {}'\
        .format(kmer_size, out_fn, input_fn)
    subprocess.check_call(cmd, shell=True)
    if os.path.exists(out_fn):
        cmd = 'jellyfish stats {}'.format(out_fn)
    else:
        cmd = 'jellyfish stats {}_0'.format(out_fn)
    res = subprocess.check_output(cmd, shell=True)
    return res


if __name__ == "__main__":

    # parse input data
    if len(sys.argv) != 2:
        print("Must include input file as command line argument. \n"
              "Quitting.")
        sys.exit(1)
    inp_fn = sys.argv[1]
    with open(inp_fn) as inp_file:
        dna_seqs = parse_fasta(inp_file)
    # extract k-mers of length 15 and print the results
    kmer_len = 15
    kmers = extract_kmers(dna_seqs, skip_unknown=True, k=kmer_len)
    print_stats(kmers)
    # run the tool jellyfish and print the results 
    jelly_out = run_jellyfish(inp_fn, kmer_size=kmer_len)
    print("JELLYFISH OUTPUT")
    print(str(jelly_out, encoding='utf-8'))

