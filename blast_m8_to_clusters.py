#!/usr/bin/env python3.3

import sys
import utils
import copy
import os
import argparse
import fastn
import blast

parser = argparse.ArgumentParser(description = 'Clusters sequences from a blast hit file, where one input file was blasted against itself', usage='%(prog)s [options] <in.fasfa> <in.blast.m8> <outfile>')
parser.add_argument('--min_length', type=int, help='Minimum length of hits [%(default)s]', default=1, metavar='INT')
parser.add_argument('--min_pc_id', type=float, help='Minimum percent identity of hits [%(default)s]', default=0, metavar='FLOAT in [0,100]')
parser.add_argument('--min_length_proportion', type=float, help='This option is used if the blast intput file has the lengths of the query and ref added as extra columns (see blast_m8_add_seq_lengths.py). If this is set to X, hits will only be used if hit length / (length of shortest sequence) >= X [%(default)s]', default = 0.95, metavar='FLOAT in [0,1]')
parser.add_argument('seq_fastn', help='Name of fasta/q file of sequences that were blasted', metavar='in.fasta')
parser.add_argument('blastfile', help='Name of input file made by blastall with the -m 8 option', metavar='in.blast.m8')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

cluster_count = 0
clusters = [] # list of sets
strands = {}  # key: sequence name. value: '+' or '-' for the strands
all_seqs = {}
fastn.file_to_dict(options.seq_fastn, all_seqs)
unused_seqs = set(all_seqs.keys())

blast_reader = blast.file_reader(options.blastfile)
for hit in blast_reader:
    # check %identity
    if hit.percent_identity < options.min_pc_id or hit.alignment_length < options.min_length or hit.qry_name == hit.ref_name:
        continue

    # check length of hit ocmpared with the length of the sequences
    if hit.qry_length is None or hit.ref_length is None or (1.0 * hit.alignment_length / min(hit.qry_length, hit.ref_length)) < options.min_length_proportion:
        continue

    qry_index = None
    ref_index = None

    for i in range(len(clusters)):
        if hit.qry_name in clusters[i]:
            qry_index = i
        if hit.ref_name in clusters[i]:
            ref_index = i

        if qry_index is not None and ref_index is not None:
            break

    # work out the strands (qstart is always < qend)
    if hit.ref_start > hit.ref_end:
        same_strand = False
    else:
        same_strand = True

    if hit.qry_name in strands and hit.ref_name in strands:
        if (strands[hit.qry_name] == strands[hit.ref_name]) != same_strand:
            #print('Mismatching strand info', hit, strands[hit.qry_name], strands[hit.ref_name], sep='\n', file=sys.stderr)
            #sys.exit(1)
            if qry_index is not None  and ref_index is not None and qry_index == ref_index:
                print('Mismatching strand info', options.blastfile, hit, strands[hit.qry_name], strands[hit.ref_name], clusters[qry_index], sep='\n', file=sys.stderr)
                sys.exit(1)
            else:
                for x in clusters[qry_index]:
                    strands[x] = {'+':'-', '-':'+'}[strands[x]]

    elif hit.qry_name in strands:
        if same_strand:
            strands[hit.ref_name] = strands[hit.qry_name]
        else:
            strands[hit.ref_name] = {'+':'-', '-':'+'}[strands[hit.qry_name]]
    elif hit.ref_name in strands:
        if same_strand:
            strands[hit.qry_name] = strands[hit.ref_name]
        else:
            strands[hit.qry_name] = {'+':'-', '-':'+'}[strands[hit.ref_name]]
    else:
        strands[hit.ref_name] = '+'
        if same_strand:
            strands[hit.qry_name] = '+'
        else:
            strands[hit.qry_name] = '-'

    #print(strands[hit.qry_name], strands[hit.ref_name], hit)

    if qry_index is not None  and ref_index is not None:
        clusters[qry_index].add(hit.ref_name)
        if qry_index != ref_index:
            clusters[qry_index].update(clusters[ref_index])
            clusters.pop(ref_index)
    elif qry_index is not None:
        clusters[qry_index].add(hit.ref_name)
    elif ref_index is not None:
        clusters[ref_index].add(hit.qry_name)
    else:
        clusters.append(set([hit.qry_name, hit.ref_name]))


# print clusters ordered by size of cluster
clusters.sort(key=len, reverse=True)
f = utils.open_file_write(options.outfile)

for i in range(len(clusters)):
    for x in clusters[i]:
        unused_seqs.discard(x)
        print(i+1, x, strands[x], sep='\t', file=f)


counter = -1
for x in sorted(unused_seqs):
    print(counter, x, sep='\t', file=f)
    counter -= 1

utils.close(f)


assembled_seqs = []

for i in range(len(clusters)):
    reads_file = options.outfile + '.cluster.' + str(i+1)
    f = utils.open_file_write(reads_file)
    for id in clusters[i]:
        seq = all_seqs[id]
        if strands[id] == '-':
            seq = copy.copy(all_seqs[id])
            seq.revcomp()
        else:
            seq = all_seqs[id]

        print(seq, file=f)
    utils.close(f)

    utils.syscall('cap3 ' + reads_file)
    singlet_count = fastn.count_sequences(reads_file + '.cap.singlets')
    contig_count = fastn.count_sequences(reads_file + '.cap.contigs')
    if singlet_count == 0 and contig_count == 1:
        seq_reader = fastn.file_reader(reads_file + '.cap.contigs')
        for seq in seq_reader:
            seq.id = 'cluster.' + str(i+1) + '.contig'
            assembled_seqs.append(copy.copy(seq))

        for e in ['ace', 'contigs.links', 'contigs.qual', 'info', 'singlets', 'contigs']:
            os.unlink(reads_file + '.cap.' + e)
        os.unlink(reads_file)
    else:
        print('Got', singlet_count, 'singlets and', contig_count, 'contigs', file=sys.stderr)

f = utils.open_file_write(options.outfile + '.clustered.fa')
for id in sorted(unused_seqs):
    print(all_seqs[id], file=f)
for seq in assembled_seqs:
    print(seq, file=f)
utils.close(f)
