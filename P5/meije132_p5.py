#!/usr/bin/env python3
"""Name: David Meijer
Student number: 930108565100

Description: this script takes a FASTA file with a reference protein
sequence and a FASTA file with multiple related protein sequences. The
script will align the reference sequence with every related sequence
using Needle and return a table with original sequence length, the
Hamming distance between the two aligned sequences and the percentage
identity between the two aligned sequences.
"""

# Space to import packages:
from sys import argv
import sys
import subprocess
import os.path

# Space for definitions:
def parsing(infile):
    """Parsing infile in fasta format as {header:seq}.
    
    infile -- file object, lines describing fasta records.
    seq_dict -- dictionary, contains parsed infile records as 
    {header:seq}.
    """
    seq_dict = {}
    for line in infile:
        if line.startswith(">"):
            seq_header = line.strip()
            if seq_header not in seq_dict:
                seq_dict[seq_header] = []
                continue
        else:
            seq = line.strip()
            seq_dict[seq_header].append(seq)
    for seq_name, seq in seq_dict.items():
        seq_dict[seq_name] = ["".join(seq)]
    return(seq_dict)
    
def add_seqlength(infile_dict):
    """Takes parsed fasta records as dict and adds seq length.
    
    infile_dict -- dictionary, contains parsed infile reocords as 
    {header:[info]}. 
    dict_out -- dictionary, contains parsed infile fasta records 
    and seq length of the fasta record as {header:[seq,seq_length]}.
    """
    for header, info in infile_dict.items():
        seq_length = len(info[0])
        info.append(seq_length)
    dict_out = infile_dict
    return(dict_out)
    
def run_needle(ref_file, related_file, gapopen=8, gapextend=0.5):
    """Align seq in ref_file to seqs in related_file.
    
    ref_file -- file object, lines desribing fasta record.
    related_file -- file object, lines describing multiple fasta
    records.
    gapopen -- integer, gap open pentaly, default is 8.
    gapextend -- integer, gap extension pentaly, default is 0.5.
    """
    name_outfile = "out.needle"
    if os.path.isfile("./out.needle"):
        print("out.needle file already exists. Needle did not run!"\
        , file=sys.stderr)
        pass
    else:
        cmd = "needle -asequence {0} -bsequence {1} -gapopen {2} \
        -gapextend {3} -outfile {4} -aformat multiple".format(\
        ref_file, related_file, gapopen, gapextend, name_outfile)
        e = subprocess.check_call(cmd, shell=True)
        print("EXIT STATUS NEEDLE IS ",e, " AND EXIT TYPE IS ", \
        type(e), file=sys.stderr)

def calc_hamming_dist(seq1, seq2):
    """Counts mismatches between two aligned sequences of equal length.
    
    seq1 -- string, DNA/RNA/prot sequence.
    seq2 -- string, DNA/RNA/prot sequence.
    """
    count = 0
    hamming_list = list(zip(seq1, seq2))
    for elem in hamming_list:
        if elem[0] != elem[1]:
            count += 1
    return(count)

def perc_identity(aligned_seq, hamming_dist):
    """Calculates percentage identity of aligned sequences.
    
    aligned_seq -- string, one of the aligned seqs have the same length
    after alignment. 
    hamming_dist -- integer, mismatch count between two sequences. 
    """
    aligned_count = len(aligned_seq)
    matches = aligned_count - hamming_dist
    perc_identity = (matches/aligned_count)*100
    return(perc_identity)

def get_needle_records(infile):
    """Reads infile and passes on one aligned record at a time.
    
    infile -- file object, multiple needle alignments.
    """
    alignment_record = []
    for line in infile:
        if line.startswith("#=") or line.startswith("#-"):
            yield(alignment_record)
            alignment_record = []
        elif len(line) < 2:
            pass
        elif not line.startswith("#"):
            alignment_record.append(line.strip())
    if alignment_record:
        yield(alignment_record)

def parse_needle_records(record):
    """Takes alignment record and returns aligned sequences.
    
    record -- list of string, lines belonging to a single record.
    ref_seq -- string, aligned reference sequence.
    ali_seq -- string, aligned to reference sequence.
    seq_label_ref -- string, sequence label of reference sequence.
    seq_label_ali -- string, sequence label of aligned sequence to
    reference sequence.
    """
    count = 1
    seq_labels = []
    ref_seq = ""
    ali_seq = ""
    for line in record:
        line = line.strip()
        if count == 1:
            seq_label_ref = line[0:10]
            ref_seq += line[21:71]
            count += 1
            continue
        if count == 2:
            count += 1
            continue
        if count == 3:
            seq_label_ali = line[0:10]
            ali_seq += line[21:71]
            count = 1
    ref_seq_list = ref_seq.partition(" ")
    ali_seq_list = ali_seq.partition(" ")
    ref_seq = ref_seq_list[0]
    ali_seq = ali_seq_list[0]
    return(ref_seq, ali_seq, seq_label_ref, seq_label_ali)
    
# Space for the main code:
if __name__ == "__main__":
    
    # Open input ref file:
    with open(argv[1]) as ref_file:
        # Parse input ref file:
        ref_dict = parsing(ref_file)

    # Open input related file:
    with open(argv[2]) as related_file:
        # Parse input related file:
        related_dict = parsing(related_file)
    
    # Add sequence length to sequences in dicts and overwrite:
    related_dict = add_seqlength(related_dict)
    ref_dict = add_seqlength(ref_dict)
    
    # Run EMBOSS needle on ref.fasta and related.fasta:
    run_needle(argv[1],argv[2])

    # Print header for results table:
    print("Sequence1\tLength\tSequence2\tLength\tHamm\tIdent")
    
    # Open Needle file:
    with open("out.needle") as needle_file:
        # Get aligned sequences records from Needle file:
        for alignment_record in get_needle_records(needle_file):
            # No empty records!
            if len(alignment_record) > 1:
                # Get relevant parts from the record:
                ref_seq, ali_seq, seq_label_ref, seq_label_ali = \
                parse_needle_records(alignment_record)
                
                # Calculate Hamming distance between two seqs:
                hamm_dist = calc_hamming_dist(ref_seq, ali_seq)
                
                # Calculate percentage identity between two seqs:
                perc_ident = round(perc_identity(ref_seq, hamm_dist),2)
                
                # Find unaligned length for ref seq with header:
                for header, info in ref_dict.items():
                    if seq_label_ref[0:8] in header:
                        ref_length = info[1]
                
                # Find unaligned length for aligned seq with header:
                for header, info in related_dict.items():
                    if seq_label_ali[0:8] in header:
                        ali_length = info[1]

                # Printing the results into a table:
                print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(\
                seq_label_ref, ref_length, seq_label_ali, ali_length, \
                hamm_dist, perc_ident))
                
                
                
                
                

    
