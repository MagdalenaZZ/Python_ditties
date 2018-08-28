#!/usr/bin/env python3.3

import sys
import argparse
import utils

parser = argparse.ArgumentParser(
    description = 'Combines failes line-by-line, can add, take min or max of columns (or do nothing)',\
    usage = '%(prog)s [options] <comma separated list of actions in min,max,same,sum> <outfile> <list of files>')
parser.add_argument('actions', help='Comma separated list of what to do to each column.  Can be same,min,max or sum')
parser.add_argument('outfile', help='Name of output file')
parser.add_argument('files', nargs='*', help='List of files to combine')
options = parser.parse_args()

actions = options.actions.split(',')
allowed_actions = set(['min', 'max', 'sum', 'same'])
for a in actions:
    if a not in allowed_actions:
        print('Action not recognised:', a, file=sys.stderr)
        sys.exit(1)

utils.file_column_combine(actions, options.files, options.outfile)
