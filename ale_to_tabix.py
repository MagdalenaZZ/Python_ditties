#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Convert ALE output file to tabix file, so that plots for artemis can be made easliy',
    usage = '%(prog)s <infile> <outfile>')
parser.add_argument('infile', help='Name of ALE output file')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

fin = utils.open_file_read(options.infile)
fout = utils.open_file_write(options.outfile)

for line in fin:
    if line.startswith('#'):
        if line.startswith('# Reference:'):
            ref_name = line.split()[2]
    else:
        a = line.rstrip().split()
        a[0] = ref_name
        a[1] = str(int(a[1])+1)
        print('\t'.join(a), file=fout)

utils.close(fin)

