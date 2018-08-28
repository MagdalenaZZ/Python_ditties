#!/usr/bin/env python3.3


import argparse
import re
import copy
import os
import sys
from operator import attrgetter

import utils
import fastn
import sam
import external_progs

parser = argparse.ArgumentParser(
    description = 'Works out the layout of the contigs within scaffolds, using the file *.tags_and_sam.gz file made by the script scaffold_test_check_using_tags.py',
    usage = '%(prog)s [options] <infile> <outfile>')
parser.add_argument('infile', help='Name of *.tags_and_sam.gz file')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

WRONG_ORIENTATION = 1
WRONG_CHR = 2
WRONG_DIST = 4
WRONG_ORDER = 8
NA = 16
DUP = 32


f = utils.open_file_read(options.infile)
lines = f.readlines()
utils.close(f)

nodes = {}

i = 0
for i in range(0, len(lines), 3):
    flag = int(lines[i])
    sam1 = sam.SamRecord(lines[i+1].strip())
    sam2 = sam.SamRecord(lines[i+2].strip())

    if sam1.id not in nodes:
        nodes[sam1.id] = set()


    if sam2.id not in nodes:
        nodes[sam2.id] = set()

    if sam1.rname == sam2.rname:

        if sam1.pos < sam2.pos:
            if sam1.query_strand() == sam2.query_strand() == '+':
                nodes[sam1.id].add(sam2.id)
            elif sam1.query_strand() == sam2.query_strand() == '-':
                nodes[sam2.id].add(sam1.id)
            else:
                nodes[sam1.id].add(sam2.id)
                nodes[sam2.id].add(sam1.id)

        else:
            if sam1.query_strand() == sam2.query_strand() == '+':
                nodes[sam2.id].add(sam1.id)
            elif sam1.query_strand() == sam2.query_strand() == '-':
                nodes[sam1.id].add(sam2.id)
            else:
                nodes[sam1.id].add(sam2.id)
                nodes[sam2.id].add(sam1.id)

print(nodes)


cmd = 'echo "digraph G {'
first  = True

for node, l in sorted(nodes.items()):
    if first:
        first = False
    else:
        cmd += ';'

    if len(l):
        cmd += ';'.join([node + '->' + x for x in l])
    else:
        cmd += node


cmd += '}"  | dot -Tpdf > ' + options.outprefix + '.pdf'

make_graph = options.outprefix + '.make_graph.sh'
f = utils.open_file_write(make_graph)
print(cmd, file=f)
utils.close(f)
utils.syscall('bash ' + make_graph)

