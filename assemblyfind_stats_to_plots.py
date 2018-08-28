#!/usr/bin/env python3

import argparse
import os
from fastaq import utils

parser = argparse.ArgumentParser(
    description = 'Makes plots etc from spreadsheet made by assemblyfind',
    usage = '%(prog)s [options] <in.csv> <outprefix')
parser.add_argument('infile')
parser.add_argument('outprefix')
parser.add_argument('--ref', help='Reference fasta file')
options = parser.parse_args()



def get_ref_stats(options):
    if not options.ref:
        return None

    stats = {}
    tmp = utils.syscall_get_stdout('stats -S ' + options.ref)
    for line in tmp:
        key, val = line.rstrip().split()[1:]
        if key == 'mean_length':
            val = float(val)
        else:
            val = int(val)

        stats[key] = val

    return stats

def get_assembly_stats(options):
    f = utils.open_file_read(options.infile)
    csv_headers = []
    stats = {}
    float_headers = set([
        'Avg Contig Length',
        'Average Quality',
        'Insert Size Average',
        'Insert Size Std Dev'
    ])

    for line in f:
        if len(csv_headers) == 0:
            csv_headers = line.rstrip().split('\t')[2:]
            stats = {k:[] for k in csv_headers}
        else:
            data = line.rstrip().split('\t')[2:]
            assert len(data) == len(csv_headers) == len(stats)
            for i in range(len(data)):
                if csv_headers[i] in float_headers:
                    stats[csv_headers[i]].append(float(data[i]))
                else:
                    stats[csv_headers[i]].append(int(data[i]))

    utils.close(f)
    return stats


def plot_assembly_stat(options, ref_stats, assembly_stats, stat):
    x_in_mb = stat in ['Total Length']

    if ref_stats is not None:
        ref_length = ref_stats['total_length'] / 1000000 if x_in_mb else ref_stats['total_length']
    if x_in_mb:
        r_vector = 'c(' + ','.join([str(x/1000000) for x in assembly_stats[stat]]) + ')'
    else:
        r_vector = 'c(' + ','.join([str(x) for x in assembly_stats[stat]]) + ')'

    files_prefix = options.outprefix + '.' + stat.translate(str.maketrans(' ()', '___'))
    r_script = files_prefix + '.R'
    f = utils.open_file_write(r_script)
    print('x=' + r_vector, file=f)
    print('', file=f)
    print('pdf("' + files_prefix + '.pdf")', file=f)
    print('hist(x, xlab="', stat, '", main="", breaks=100)', sep='', file=f)
    print('dev.off()', file=f)
    utils.close(f)
    utils.syscall('R CMD BATCH ' + r_script)
    os.unlink(r_script)
    os.unlink(r_script + 'out')



ref_stats = get_ref_stats(options)
assembly_stats = get_assembly_stats(options)


if ref_stats is not None:
    assembly_stats['Total Length Percent'] = [round(100 * float(x / ref_stats['total_length']), 2) for x in assembly_stats['Total Length']]


assembly_stats['Reads Mapped Percent'] = [round(100 * float(assembly_stats['Reads Mapped'][i] / assembly_stats['Total Raw Reads'][i]), 2) for i in range(len(assembly_stats['Reads Mapped']))]

assembly_stats['Reads Paired Percent'] = [round(100 * float(assembly_stats['Reads Paired'][i] / assembly_stats['Total Raw Reads'][i]), 2) for i in range(len(assembly_stats['Reads Paired']))]


for stat in assembly_stats:
    plot_assembly_stat(options, ref_stats, assembly_stats, stat)

