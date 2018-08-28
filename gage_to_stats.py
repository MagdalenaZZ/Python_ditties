#!/usr/bin/env python3


import argparse
import re
import copy
import os
import sys
from operator import attrgetter

import utils

parser = argparse.ArgumentParser(
    description = 'Gets GAGE stats from bsub stdout file',
    usage = '%(prog)s <gage.o>')
parser.add_argument('infile', help='Name of input gage.o bsub stout file')
options = parser.parse_args()


contigs = -1
scaffs = -1
contig_N50 = -1
scaff_N50 = -1
contig_corr_N50 = -1
scaff_corr_N50 = -1
contig_errs = -1
scaff_errs = -1


f = utils.open_file_read(options.infile)
lines = f.readlines()
utils.close(f)


i = 0
while i < len(lines):
    line = lines[i].rstrip()

    if line == 'Contig Stats':
        contigs = int(lines[i+1].split()[-1])
        if lines[i+8].startswith('N50'):
            contig_N50 = int(lines[i+8].split()[1])
    elif line == 'Scaffold Stats':
        scaffs = int(lines[i+1].split()[-1])
        if lines[i+8].startswith('N50'):
            scaff_N50 = int(lines[i+8].split()[1])
    elif line == 'Corrected Contig Stats':
        if lines[i+8].startswith('N50'):
            contig_corr_N50 = int(lines[i+8].split()[1])
    elif line == 'Corrected Scaffold Stats':
        data = lines[i+1].split(',')
        scaff_errs = int(data[1])
        scaff_corr_N50 = int(data[-1])

    i += 1

#contig_N50 = round(0.001 * contig_N50, 1)
#scaff_N50 = round(0.001 * scaff_N50)
#contig_corr_N50 = round(0.001 * contig_corr_N50, 1)
#scaff_corr_N50 = round(0.001 * scaff_corr_N50)



print(contigs, contig_N50, contig_errs, contig_corr_N50, scaffs, scaff_N50, scaff_errs, scaff_corr_N50, sep='\t')
