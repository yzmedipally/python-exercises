#!/usr/bin/env python3

from sys import argv
import subprocess as sbp
import os.path


def decode_qual(seq_quality):
    qual_list = []
    for char in seq_quality:
        qual_list.append(ord(char) - 66)
    return qual_list



def parse_FASTQ(lines):
    sequence_dic={}
    plus=False
    for line in lines:
        if line.startswith('@'):
            plus=False
            current_label=line[1:].strip()
            sequence_dic[current_label]=['','',0]
        elif line.startswith('+'):
            plus=True
        elif plus:
            sequence_dic[current_label][1]=decode_qual(line.strip())
        elif not plus:
            sequence_dic[current_label][0]=line.strip()
            sequence_dic[current_label][2]=len(line.strip())
    return sequence_dic
    

def calc_seq_info(seq_dic):
    total_len = 0
    min_len, max_len = None, 0
    for val in seq_dic.values():
        total_len += val[2]
        if max_len < val[2]:
            max_len = val[2]
        elif min_len is None or min_len > val[2]:
            min_len = val[2]
    avg_len = total_len / len(seq_dic.keys())
    return (min_len, max_len, avg_len)


def calc_avg_qual(seq_dic):
    quality_list = [0] * calc_seq_info(seq_dic)[1]
    amount_list=[0]* calc_seq_info(seq_dic)[1]
    for val in seq_dic.values():
        for i in range(len(val[1])):
            quality_list[i] += val[1][i]
            amount_list[i]+=1
    avg_quality=[]
    for i in range(len(amount_list)):
        avg_quality.append(quality_list[i]/amount_list[i])
    return avg_quality

def comparing_dic(seq_dic,trim_dic):
    dif_qual=[]
    trim_qual=calc_avg_qual(trim_dic)
    seq_qual=calc_avg_qual(seq_dic)
    for i in range(len(seq_qual)):
        dif_qual.append(trim_qual[i]-seq_qual[i])
    return dif_qual
    
if __name__=="__main__":
    if len(argv)>1:
        if not os.path.exists('trimmed.fq'):
            e=sbp.check_call('fastq_quality_trimmer -t 30 -Q 64 -i '+argv[1]+' -o trimmed.fq', shell=True)
            if e:
                raise Exception("Can't trimm the sequences")
        else:
            user_input = input("Use exisiting trimmed.fq file? (y/n)")
            if not user_input in "Yy":
                e=sbp.check_call('fastq_quality_trimmer -t 30 -Q 64 -i '+argv[1]+' -o trimmed.fq', shell=True)
                if e:
                    raise Exception("Can't trimm the sequences")
                
        with open(argv[1],'r') as filex:
            seq_dic=parse_FASTQ(filex)
        if os.path.exists('trimmed.fq'):
            with open('trimmed.fq') as filey:
                trim_dic=parse_FASTQ(filey)
    else:
        raise Exception('Not enought arguments')
    
    print('ORIGINAL: min=%s, max=%s, avg=%s' % calc_seq_info(seq_dic))
    print('TRIMMED: min=%s, max=%s, avg=%s' % calc_seq_info(trim_dic))
    #print("{0:d}\t{1:d}\t{2:f}".format(calc_seq_info(trim_dic)))
    
    avg_qual_og = calc_avg_qual(seq_dic)
    avg_qual_trim = calc_avg_qual(trim_dic)
    avg_qual_diff = comparing_dic(seq_dic, trim_dic)
    for i in range(len(avg_qual_og)):
        print("%d\t%.2f\t%.2f\t%.2f" % (i+1, avg_qual_trim[i], avg_qual_og[i], avg_qual_diff[i]))
    #print("\t".join(list(map(str,comparing_dic(seq_dic,trim_dic)))))
    
    
    

