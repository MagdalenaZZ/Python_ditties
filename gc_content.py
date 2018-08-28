#!/usr/bin/env python


import sys
import argparse
import re

# Describe what the script does
parser = argparse.ArgumentParser(description='This script calculates GC content of fasta sequences', epilog= 'Please enjoy!')
parser.add_argument('-i', '--input', default=None, dest='fasta', action='store', help="FASTA file")

args  = parser.parse_args()
#print(args.fasta)


a_list = []
#a_list.append((1, 2))   


# Read in fasta file
with open(args.fasta) as f:
    for line in f:
        line = line.rstrip('\n')
        if re.match(r'^>', line):
            seq = f.readline()
            seq=seq.rstrip('\n')
            t=(line, seq)
            a_list.append((line,seq))
            #print ("head %s seq %s" % t)
        else:
            print ("Should not happen")

# open output 

of = args.fasta + ".gc"
out = open(of, 'w')

################################
# Calculate gc content

def compute_gc(sequence):
    gc = 0.0
    sequence_length = len(sequence)
    sequence = sequence.upper()
    for nt in sequence:
        if nt == 'C' or nt == 'G':
            gc += 1
    return round(gc / sequence_length,3)

#############################

for tup in a_list:
    gc = compute_gc(tup[1])
    out.write(''.join((tup[0],'\t',str(gc),'\n')))


    




    
