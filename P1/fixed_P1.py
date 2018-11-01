#!/usr/bin/env python3
"""Script to count the number of k-mers in a FASTA file

Output is compared to jellyfish output
Author: D.E. Bug; 
fixed by: Hidde Bleeker
931202071020

Difference:
1a2
> '''Script to count the number of k-mers in a FASTA file
3,5d3
< '''
< Author: D.E. Bug
< Script to count the number of k-mers in a FASTA file
6a5,7
> Author: D.E. Bug;
> fixed by: Hidde Bleeker
> 931202071020
9c10
< from sys import argv
---
> import sys
12a14
>
20c22
<         if not line.stripit():
---
>         if not line.strip():
26c28
<             seqs[label] += line.strip()[1:]
---
>             seqs[label] += line.strip()
28a31
>
34c37
<     skip_unknown: bool, specifying wheter k-mers containing
---
>     skip_unknown: bool, specifying whether k-mers containing
38,40c41,43
<     ch = {} #dict to store characters
<     res = {} #dict to store k-mers and counts
<     for label, seq in seqs.items()
---
>     ch = {}  # dict to store characters
>     res = {}  # dict to store k-mers and counts
>     for label, seq in seqs.items():
45,46c48,49
<         for i in range(len(seq)-kmer_size):
<             kmer = seq[i:i+kmer_size]
---
>         for i in range(len(seq) - kmer_size + 1):
>             kmer = seq[i:i + kmer_size]
51c54
<             if skip_unknown is True and count is False:
---
>             if skip_unknown and not count:
57a61
>
65c69
<     unique = [i for i, j in res.items() if j==1]
---
>     unique = [i for i, j in res.items() if j == 1]
74c78
<         print(k, v)
---
>             print(k, v)
77a82
>
85c90
<     cmd = 'jellyphish count -m {} -s 1000000 -o {} {}'\
---
>     cmd = 'jellyfish count -m {} -s 1000000 -o {} {}'\
87c92
<     e = subprocess.check_output(cmd, shell=True)
---
>     subprocess.check_call(cmd, shell=True)
92c97
<     res = subprocess.check_call(cmd, shell=True)
---
>     res = subprocess.check_output(cmd, shell=True)
94a100
>
98,99c104,110
<     inp_fn = argv[1]
<     dna_seqs = parse_fasta(open(inp_fn))
---
>     if len(sys.argv) != 2:
>         print("Must include input file as command line argument. \n"
>               "Quitting.")
>         sys.exit(1)
>     inp_fn = sys.argv[1]
>     with open(inp_fn) as inp_file:
>         dna_seqs = parse_fasta(inp_file)
102c113
<     kmers = extract_kmers(dna_seqs, skip_unknown=True, k=14)
---
>     kmers = extract_kmers(dna_seqs, skip_unknown=True, k=kmer_len)
107c118
<     print(str(jelly_out,'utf-8'))
---
>     print(str(jelly_out, encoding='utf-8'))

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

