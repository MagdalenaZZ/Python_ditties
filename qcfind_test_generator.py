#!/usr/bin/env python3.3

import fastaq
import utils
import argparse
import itertools

parser = argparse.ArgumentParser(
    description = 'Generates test code for qcfind',
    usage = '%(prog)s [options] <outprefix>')
parser.add_argument('outprefix', help='prefix of output files')
options = parser.parse_args()


options_keys = [
    'types', 
    'symlinks',
    'archives',
    'summarys',
    'levels',
    'counts',
    'assigned_directly',
]


options = {
    'types': [
        'study',
        'lane',
        'file',
        'library',
        'sample',
        'species'
    ],
    'symlinks': [True, False],
    'archives' = [True, False],
    'summarys' = [True, False],
    'levels' = [True, False],
    'counts' = [True, False],
    'assigned_directly' = [True, False],
}



# combinations of not giving an id:




