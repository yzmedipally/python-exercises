#!/usr/bin/env python3
"""Exercise P3: Parse FASTQ file, get statistics, run trimmer and compare stats

Author : Sotiria Milia, Hidde Bleeker (931202071020)

usage: python {} input.fq
"""

from sys import argv
from os.path import exists
from subprocess import check_call, CalledProcessError


def parse_fastq(fastq_file):
    """Parses FASTQ file, returning dictionary with sequences and qualities

    :param fastq_file: str, fastq format filename
    :returns: dict of dicts, with {seq_id: {Sequence: sequence (str),
    Quality: quality (list of ints)}}
    """
    fastq_dict = {}
    counter = 0
    seq_id = ""
    with open(fastq_file, "r") as fastq:
        for line in fastq:
            counter += 1
            if counter == 1:
                seq_id = line.strip()
                fastq_dict[seq_id] = {}
            elif counter == 2:
                fastq_dict[seq_id]["Sequence"] = line.strip()
            elif counter == 4:
                qual = [ord(ch) - 64 for ch in line.strip()]
                fastq_dict[seq_id]["Quality"] = qual
                counter = 0
    return fastq_dict


def min_max_avg_length(fastq_dict):
    """Returns min, max and average length of sequences in FASTQ dictionary

    :param fastq_dict: dict of dicts, with {seq_id: {Sequence: sequence (str),
    Quality: quality (list of ints)}}
    :returns: tuple of floats, as (minimum, maximum, average length)
    """
    lengths = [len(e["Sequence"]) for e in fastq_dict.values()]
    return min(lengths), max(lengths), sum(lengths) / len(lengths)


def avg_pos_score(fastq_dict):
    """Returns the average quality score in each position in a FASTQ dictionary

    :param fastq_dict: dict of dicts, with {seq_id: {Sequence: sequence (str),
    Quality: quality (list of ints)}}
    :returns: list, of floats of average quality per position over all seqs
    """
    qual_list = []
    for v in fastq_dict.values():
        for i in range(len(v["Quality"])):
            try:
                qual_list[i].append(v["Quality"][i])
            except IndexError:
                qual_list.append([v["Quality"][i], ])
    avg_list = []
    for pos in qual_list:
        avg_list.append(sum(pos) / len(pos))

    return avg_list


def run_fastq_trimmer(fastq_file, output_fn, scale=64, qual_threshold=30):
    """Runs fastq_quality_trimmer with parameters and output to given file

    :param fastq_file: str, input FASTQ filename
    :param output_fn: str, output filename of trimmed FASTQ file
    :param scale: int, the scale of the quality data
    :param qual_threshold: int, the threshold to trim qualities below
    """
    if not exists(output_fn):
        cmd = "fastq_quality_trimmer -t {} -Q {} -i {} -o {}"\
            .format(scale, qual_threshold, fastq_file, output_fn)
        try:
            check_call(cmd, shell=True)
        except CalledProcessError as e:
            print("Failed to call fastq_quality_trimmer. Error:\n{}"
                  .format(str(e)))
            exit(1)


if __name__ == "__main__":
    if not len(argv) == 2:
        print(__doc__.format(argv[0]))
        exit(1)

    output_file_name = "trimmed.fq"

    # Step 1: Parse FASTQ file
    fastq_seqs = parse_fastq(argv[1])

    # Step 2: Calculate min, max, and avg sequence length
    len_stats = min_max_avg_length(fastq_seqs)

    # Step 3: Calculate average quality score in each sequence position
    avg_qual = avg_pos_score(fastq_seqs)

    # Step 4: Run fastq_quality_trimmer program with options
    run_fastq_trimmer(argv[1], output_file_name, scale=64, qual_threshold=30)

    # Step 5: Parse trimmed FASTQ file and calculate average scores of each
    # position. Calculate improvement to untrimmed file.
    fastq_trimmed_seq = parse_fastq(output_file_name)
    trim_len_stats = min_max_avg_length(fastq_trimmed_seq)
    avg_qual_trim = avg_pos_score(fastq_trimmed_seq)

    # Step 6: Report minimum, maximum and average sequence length for both
    print("ORIGINAL: min={:d}, max={:d}, avg={:.2f}".format(*len_stats))
    print("TRIMMED: min={:d}, max={:d}, avg={:.2f}".format(*trim_len_stats))

    # Step 7: For each read position, report the average quality score in the
    # original file, the average quality score in the trimmed file, and the
    # improvement in average quality, in tab-delimited columns.
    for i, (og, trim) in enumerate(zip(avg_qual, avg_qual_trim)):
        print("{:d}\t{:.2f}\t{:.2f}\t{:.2f}\t".format(i + 1, og, trim,
                                                      trim - og))
