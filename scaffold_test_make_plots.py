#!/usr/bin/env python3.3

import sys
import argparse
import os
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Used to make plots for scaffolder testing',
    usage = '%(prog)s [options] <directory of scaff results> <output_prefix>')
parser.add_argument('--use_sspace_iter', action='store_true', help='Use this to include sspace_iter results')
parser.add_argument('dir_in', help='Input directory where scaffolding was run')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

scaffolders = ['MIP.bowtie-v0',
    'MIP.bowtie-v3',
    'MIP.bowtie2',
    'MIP.bwa',
    'OPERA.bowtie',
    'OPERA.bwa',
    'SCARPA.bowtie-v0',
    'SCARPA.bowtie-v3',
    'SCARPA.bowtie2',
    'SCARPA.bwa',
    'SGA.bowtie2',
    'SGA.bwa',
    'SOAP',
    'SOAP2',
    'SOPRA.bowtie-v0',
    'SOPRA.bowtie-v3',
    'SOPRA.bowtie2',
    'SOPRA.bwa',
    'SSPACE.bowtie-v0',
    'SSPACE.bowtie-v3']

r_symbols = {'MIP.bowtie-v0': '0',
    'MIP.bowtie-v3': '12',
    'MIP.bowtie2': '5',
    'MIP.bwa': '1',
    'OPERA.bowtie': '12',
    'OPERA.bwa': '1',
    'SCARPA.bowtie-v0': '0',
    'SCARPA.bowtie-v3': '12',
    'SCARPA.bowtie2': '5',
    'SCARPA.bwa': '1',
    'SGA.bowtie2': '5',
    'SGA.bwa': '1',
    'SOAP': '4',
    'SOAP2': '3',
    'SOPRA.bowtie-v0': '0',
    'SOPRA.bowtie-v3': '12',
    'SOPRA.bowtie2': '5',
    'SOPRA.bwa': '1',
    'SSPACE.bowtie-v0': '0',
    'SSPACE.bowtie-v3': '12'}


r_colours = {'MIP.bowtie-v0': 'black',
    'MIP.bowtie-v3': 'black',
    'MIP.bowtie2': 'black',
    'MIP.bwa': 'black',
    'OPERA.bowtie': 'red',
    'OPERA.bwa': 'red',
    'SCARPA.bowtie-v0': 'blue',
    'SCARPA.bowtie-v3': 'blue',
    'SCARPA.bowtie2': 'blue',
    'SCARPA.bwa': 'blue',
    'SGA.bowtie2': 'green',
    'SGA.bwa': 'green',
    'SOAP': 'purple',
    'SOAP2': 'purple',
    'SOPRA.bowtie-v0': 'brown',
    'SOPRA.bowtie-v3': 'brown',
    'SOPRA.bowtie2': 'brown',
    'SOPRA.bwa': 'brown',
    'SSPACE.bowtie-v0': 'orange',
    'SSPACE.bowtie-v3': 'orange'}

if options.use_sspace_iter:
    scaffolders.append('SSPACE.iter')
    r_symbols['SSPACE.iter'] = 4
    r_colours['SSPACE.iter'] = 'orange'


#r_colours = {}
#
#for k in r_symbols:
#    if k.endswith('bwa'):
#        r_colours[k] = 'red'
#    elif k.endswith('bowtie'):
#        r_colours[k] = 'blue'
#    elif k.endswith('bowtie2'):
#        r_colours[k] = 'green'
#    else:
#        r_colours[k] = 'black'

r_symbol_list = [r_symbols[x] for x in scaffolders]
r_colour_list = [r_colours[x] for x in scaffolders]
r_symbol_vector = 'c(' + ', '.join(r_symbol_list) + ')'
r_colour_vector = 'c(' + ', '.join(['"' + x + '"' for x in r_colour_list]) + ')'

r_legend = 'legend("topleft", ' \
                + 'c(' + ', '.join(['"' + x + '"' for x in scaffolders]) + '), ' \
                + 'col=' + r_colour_vector + ', ' \
                + 'pt.bg=' + r_colour_vector + ', ' \
                + 'pch=c(' + ', '.join(r_symbol_list) + ')' \
                + ')'


possible_flags = [0,1,2,4,5,8,12,16]


def get_scaff_results(dir):
    flag_counts = {k:0 for k in possible_flags}
    flag_counts['skipped'] = 0
    flag_counts['lost'] = 0

    log_file = dir + '/check_scaffolds.log'

    if os.path.exists(log_file):
        f = utils.open_file_read(dir + '/check_scaffolds.log')
        for line in f:
            a = line.split()

            if a[0].isdigit():
                flag_counts[int(a[0])] = int(a[1])
            elif a[0] in ['lost', 'skipped']:
                flag_counts[a[0]] = int(a[1])

        utils.close(f)
    else:
        print('Warning: no log file', log_file, file=sys.stderr)
        flag_counts['bad_joins'] = 0

    flag_counts['bad_joins'] = sum([flag_counts[x] for x in flag_counts.keys() if x not in [0,16, 'skipped']])

    return flag_counts





# gather all the counts for each scaffolding run
results = {k:{} for k in scaffolders}


for scaff in scaffolders:
    print(scaff)
    scaff_dir = options.dir_in + '/' + scaff
    results[scaff]['flag_counts'] = get_scaff_results(scaff_dir)

    bsub_outfile = scaff_dir + '/' + scaff.split('.')[0].lower() + '.o'
    bsub_out = utils.syscall_get_stdout('bsub-out2stats.py -s ' + bsub_outfile)
    assert len(bsub_out) == 1
    (attempt_no, exit_code, wall_hrs, cpu_secs, cpu_hrs, mem, swap, filename) = bsub_out[0].split('\t')
    assert exit_code == '0'

    results[scaff]['CPU'] = int(round(float(cpu_secs), 0))
    results[scaff]['mem'] = mem


# make a tsv file of all the stats
f = utils.open_file_write(options.outprefix + '.stats.tsv')
print('Scaffolder', 'Good joins',
      '\t'.join([str(x) for x in possible_flags if x not in [0,16]]),
      'Bad joins', 'Total joins', '% correct joins', 'Lost tags', 'Skipped tags', 'CPU', 'Mem', 'Extra CPU', 'Extra Mem', sep='\t', file=f)

for scaff in scaffolders:
    r = results[scaff]
    fc = r['flag_counts']
    total_joins = sum([fc[x] for x in possible_flags if x != 16])
    try:
        percent_good = round(100 * (1 - float(fc['bad_joins'] / total_joins)), 2)
    except ZeroDivisionError:
        percent_good = 0

    print(scaff,
          '\t'.join([str(fc[x]) for x in possible_flags if x != 16]),
          fc['bad_joins'],
          total_joins,
          percent_good,
          fc['lost'],
          fc['skipped'],
          r['CPU'], r['mem'], sep='\t', file=f)


utils.close(f)


# make a scatter plot
r_script = options.outprefix + '.plot.R'
f = utils.open_file_write(r_script)


x_coords = [results[scaff]['flag_counts'][0] for scaff in scaffolders]
y_coords = [results[scaff]['flag_counts']['bad_joins'] for scaff in scaffolders]
x_max = max(x_coords)
y_max = max(y_coords)


for type in ['pdf', 'png']:
    print(type + '("' + r_script + '.scatter.' + type + '")', file=f)

    print('plot(c(' + ','.join(str(x) for x in x_coords), '), ',
          'c(', ','.join(str(x) for x in y_coords), '), ',
          'xlab="Correct joins", ',
          'ylab="Incorrect joins", ',
          'xlim=c(0,', x_max, '), ',
          'ylim=c(0,', y_max, '), ',
          'col=', r_colour_vector, ', ',
          'pch=', r_symbol_vector, ', ',
          'bg=', r_colour_vector,
          ')', sep='', file=f)

    print(r_legend, file=f)

    print('dev.off()', file=f)

utils.close(f)


utils.syscall('R CMD BATCH ' + r_script)


