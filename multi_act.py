#!/usr/bin/env python3

import sys
import argparse
import pyfastaq
import subprocess
import os

parser = argparse.ArgumentParser(
    description = '',
    usage = '%(prog)s [options] <outprefix> <file1.fa> <file2.fa>  [<file3.fa ...]')
parser.add_argument('--blast_ops', help='blastall options [%(default)s]', default='-p blastn -m 8 -F F -e 0.01')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('fa_list', help='List of fasta files', nargs=argparse.REMAINDER)
options = parser.parse_args()
assert len(options.fa_list) > 1


def index_to_union(ops, i):
    return ops.outprefix + '.' + str(i) + '.union.fa'


def index_to_blast_file(ops, i):
    return ops.outprefix + '.blast.' + str(i) + '.vs.' + str(i+1)

# make union files
for i in range(len(options.fa_list)):
    seq = pyfastaq.sequences.Fasta('union', '')
    reader = pyfastaq.sequences.file_reader(options.fa_list[i])
    new_seq = []
    for s in reader:
        new_seq.append(s.seq)
    f = pyfastaq.utils.open_file_write(index_to_union(options, i))
    seq.seq = ''.join(new_seq)
    print(seq, file=f)
    pyfastaq.utils.close(f)


act_command = 'act ' + options.fa_list[0]

# run blast
for i in range(len(options.fa_list)-1):
    qry = index_to_union(options, i+1)
    ref = index_to_union(options, i)
    subprocess.check_output('formatdb -p F -i ' + ref, shell=True)
    cmd = ' '.join([
        'blastall', options.blast_ops,
        '-d', ref,
        '-i', qry,
        '-o', index_to_blast_file(options, i)
    ])
    subprocess.check_output(cmd, shell=True)
    act_command += ' ' + index_to_blast_file(options, i) + ' ' + options.fa_list[i+1]

os.unlink('formatdb.log')
print(act_command)
subprocess.check_output(act_command, shell=True)

