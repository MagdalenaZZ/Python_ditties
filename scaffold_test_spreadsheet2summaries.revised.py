#!/usr/bin/env python3

import sys
import argparse
import os
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Used to make plots for scaffolder testing',
    usage = '%(prog)s [options] <results.tsv> <output_prefix>')
parser.add_argument('--noplots', action='store_true', help='Don\'t make any plots')
parser.add_argument('--use_sspace_iter', action='store_true', help='Use this to include sspace_iter results')
parser.add_argument('infile', help='Input directory where scaffolding was run')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()


def run_r_script(script):
    if not options.noplots:
        utils.syscall('R CMD BATCH ' + script)
        os.unlink(script  + 'out')
        os.unlink(script)



scaffolders = [
    'ABySS',
    'Bambus2.bowtie2',
    'Bambus2.bwa',
    'MIP.bowtie-v0',
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
    'SOAP2',
    'SOPRA.bowtie-v0',
    'SOPRA.bowtie-v3',
    'SOPRA.bowtie2',
    'SOPRA.bwa',
    'SSPACE.bowtie-v0',
    'SSPACE.bowtie-v3']

r_symbols = {
    'ABySS': '3',
    'Bambus2.bowtie2': '6',
    'Bambus2.bwa': '2',
    'MIP.bowtie-v0': '0',
    'MIP.bowtie-v3': '15',
    'MIP.bowtie2': '6',
    'MIP.bwa': '2',
    'OPERA.bowtie': '15',
    'OPERA.bwa': '2',
    'SCARPA.bowtie-v0': '0',
    'SCARPA.bowtie-v3': '15',
    'SCARPA.bowtie2': '6',
    'SCARPA.bwa': '2',
    'SGA.bowtie2': '6',
    'SGA.bwa': '2',
    'SOAP2': '4',
    'SOPRA.bowtie-v0': '0',
    'SOPRA.bowtie-v3': '15',
    'SOPRA.bowtie2': '6',
    'SOPRA.bwa': '2',
    'SSPACE.bowtie-v0': '0',
    'SSPACE.bowtie-v3': '15'}


r_colours = {
    'ABySS': 'cyan2',
    'Bambus2.bowtie2': 'aquamarine4',
    'Bambus2.bwa': 'aquamarine4',
    'MIP.bowtie-v0': 'black',
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
    'SOAP2': 'purple',
    'SOPRA.bowtie-v0': 'saddlebrown',
    'SOPRA.bowtie-v3': 'saddlebrown',
    'SOPRA.bowtie2': 'saddlebrown',
    'SOPRA.bwa': 'saddlebrown',
    'SSPACE.bowtie-v0': 'orange',
    'SSPACE.bowtie-v3': 'orange'}


if options.use_sspace_iter:
    scaffolders.append('SSPACE.iter')
    r_symbols['SSPACE.iter'] = '11'
    r_colours['SSPACE.iter'] = 'orange'



r_symbol_list = [r_symbols[x] for x in scaffolders]
r_colour_list = [r_colours[x] for x in scaffolders]
r_symbol_vector = 'c(' + ', '.join(r_symbol_list) + ')'
r_colour_vector = 'c(' + ', '.join(['"' + x + '"' for x in r_colour_list]) + ')'

def r_legend(position):
    return 'legend("' + position + '", ' \
                + 'c(' + ', '.join(['"' + x + '"' for x in scaffolders]) + '), ' \
                + 'col=' + r_colour_vector + ', ' \
                + 'pt.bg=' + r_colour_vector + ', ' \
                + 'pch=c(' + ', '.join(r_symbol_list) + ')' \
                + ')'
possible_flags = [0,1,2,4,5,8,12,16]
bad_runs = set()
#bad_runs.add(('3D7.combi', 'ABySS'))
bad_runs.add(('3D7.combi', 'Bambus2.bowtie2'))
bad_runs.add(('3D7.combi', 'Bambus2.bwa'))
#bad_runs.add(('human.long', 'ABySS'))
bad_runs.add(('human.long', 'SOPRA.bowtie2'))
#bad_runs.add(('human.combi', 'ABySS'))
bad_runs.add(('human.combi', 'MIP.bowtie2'))
bad_runs.add(('human.combi', 'MIP.bwa'))
bad_runs.add(('human.combi', 'SOPRA.bowtie2'))
bad_runs.add(('human.combi', 'SOPRA.bwa'))

class Scaffold_result:


    def __init__(self, headers, results_list):
        self.headers = headers.rstrip().split('\t')
        self.results = {}
        self.data_type = None

        for r in results_list:
            res = dict(zip(self.headers, r.rstrip().split('\t')))

            if self.data_type is None:
                self.data_type = res['Dataset']
            else:
                assert self.data_type == res['Dataset']

            del res['Dataset']
            scaffolder = res['Scaffolder']
            del res['Scaffolder']

            try:
                res['% correct joins'] = str(int(res['Correct joins']) / (int(res['Bad joins']) + int(res['Correct joins'])))
            except:
                res['% correct joins'] = '0'
            self.results[scaffolder] = res


    def __repr__(self):
        s = '\t'.join(self.headers) + '\n'  + self.data_type + '\n'
        for scaff in self.results:
            s += scaff + '\t' + repr(self.results[scaff]) + '\n'

        return s.rstrip()


    def get_r_vectors(self):
        symbol_list = [r_symbols[x] for x in scaffolders if (self.data_type, x) not in bad_runs]
        colour_list = [r_colours[x] for x in scaffolders if (self.data_type, x) not in bad_runs]
        symbol_vector = 'c(' + ', '.join(symbol_list) + ')'
        colour_vector = 'c(' + ', '.join(['"' + x + '"' for x in colour_list]) + ')'
        return(symbol_vector, symbol_list, colour_vector, colour_list)


    def plot_scatter(self, stat1, stat2, outprefix, legend=False, main=''):
        r_script = outprefix + '.R'
        f = utils.open_file_write(r_script)

        x_coords = [int(self.results[scaff][stat1]) for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        y_coords = [int(self.results[scaff][stat2]) for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        x_max = max(x_coords)
        y_max = max(y_coords)
        r_syms_v, r_syms_l, r_cols_v, r_cols_l = self.get_r_vectors()

        for type in ['pdf', 'png', 'svg']:
            print(type + '("' + outprefix + '.' + type + '")', file=f)

            print('plot(c(' + ','.join(str(x) for x in x_coords), '), ',
                  'c(', ','.join(str(x) for x in y_coords), '), ',
                  'xlab="', stat1, '", ',
                  'ylab="', stat2, '", ',
                  #'xlim=c(0,', x_max, '), ',
                  #'ylim=c(0,', y_max, '), ',
                  'main="', main, '",',
                  'col=', r_cols_v, ', ',
                  'pch=', r_syms_v, ', ',
                  'bg=', r_cols_v,
                  ')', sep='', file=f)

            if legend:
                print(r_legend('topleft'), file=f)

            print('dev.off()', file=f)

        utils.close(f)
        run_r_script(r_script)


    def plot_scatter_with_chulls(self, stat1, stat2, outprefix, legend=False, main=''):
        r_script = outprefix + '.R'
        f = utils.open_file_write(r_script)


        #for type in ['pdf', 'png']:
        for type in ['pdf', 'svg']:
            x_coords = [float(self.results[scaff][stat1]) for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
            y_coords = [float(self.results[scaff][stat2]) for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
            x_max = max(x_coords)
            y_max = max(y_coords)
            x_min = min(x_coords)
            y_min = min(y_coords)
            x_extra =  (x_max - x_min) * 0.03
            y_extra =  (y_max - y_min) * 0.03
            r_syms_v, r_syms_l, r_cols_v, r_cols_l = self.get_r_vectors()

            if stat1 in ['TPR', 'FDR']:
                x_max = 1
                y_max = 1

            if type == 'pdf':
                print(type + '("' + outprefix + '.' + type + '", useDingbats=F)', file=f)
            else:
                print(type + '("' + outprefix + '.' + type + '")', file=f)


            print('plot(c(' + ','.join(str(x) for x in x_coords), '), ',
                  'c(', ','.join(str(x) for x in y_coords), '), ',
                  'xlab="', stat1, '", ',
                  'ylab="', stat2, '", ',
                  #'xlim=c(0,', x_max, '), ',
                  #'ylim=c(0,', y_max, '), ',
                  'main="', main, '",',
                  'col=', r_cols_v, ', ',
                  'pch=', r_syms_v, ', ',
                  'bg=', r_cols_v,
                  ')', sep='', file=f)

            # get x and y coords grouped by colour
            x_coords = {x:[] for x in set(r_cols_l)}
            y_coords = {x:[] for x in set(r_cols_l)}

            for scaff in scaffolders:
                if (self.data_type, scaff) in bad_runs:
                    continue
                x_coords[r_colours[scaff]].append(float(self.results[scaff][stat1]))
                y_coords[r_colours[scaff]].append(float(self.results[scaff][stat2]))


            print(r'''
col2trans = function(colour) {
    x=as.vector(col2rgb(colour))
    return(rgb(x[1],x[2],x[3],20,maxColorValue=255))
}

''', file=f)


            for col in set(r_cols_l):
                xpts = []
                ypts = []
                if x_extra == 0:
                    x_extra = 5
                if y_extra == 0:
                    y_extra = 0.1

                for p in x_coords[col]:
                    xpts.extend([p-x_extra,p,p,p+x_extra])
                for p in y_coords[col]:
                    ypts.extend([p,p-y_extra,p+y_extra,p])

                xpts = 'c(' + ','.join([str(x) for x in xpts]) + ')'
                ypts = 'c(' + ','.join([str(x) for x in ypts]) + ')'

                print('pts=structure(list(x=', xpts, ',y=', ypts, '))', sep='', file=f)
                #print('points(pts, col="', r_colours[scaff], '")', file=f)
                print('chuld = lapply(pts, "[", chull(pts))', file=f)
                #print('polygon(spline.poly(as.matrix(as.data.frame(chuld)),100), border="', col, '",lwd=2)', file=f)
                print('polygon(chuld, border=col2trans("', col, '"), col=col2trans("', col, '"), lwd=2)', sep='', file=f)

            if legend:
                print(r_legend('topright'), file=f)

            print('dev.off()', file=f)

        utils.close(f)
        run_r_script(r_script)


    def barplot_two_stacked(self, stat1, stat2, y_label, outprefix, main=''):
        r_script = outprefix + '.R'
        f = utils.open_file_write(r_script)
        r_syms_v, r_syms_l, r_cols_v, r_cols_l = self.get_r_vectors()

        bar_heights1 = [self.results[scaff][stat1] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        bar_heights2 = [self.results[scaff][stat2] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        #names = 'c(' + ','.join(['"' + x + '"' for x in scaffolders]) + ')'

        all_data = [(int(bar_heights1[i]) + int(bar_heights2[i]), bar_heights1[i], bar_heights2[i], scaffolders[i], r_colours[scaffolders[i]]) for i in range(len(bar_heights1)) if (self.data_type, scaffolders[i]) not in bad_runs]
        all_data.sort()
        bar_heights1 = [all_data[i][1] for i in range(len(all_data))]
        bar_heights2 = [all_data[i][2] for i in range(len(all_data))]
        names = 'c(' + ','.join(['"' + all_data[i][3] + '"' for i in range(len(all_data))]) + ')'
        #cols = [r_colours[x] for x in names]

        print(r''' col2trans = function(colour) {
    x=as.vector(col2rgb(colour))
    return(rgb(x[1],x[2],x[3],20,maxColorValue=255))
}
''', file=f)
        print('cols=', r_cols_v, file=f)
        print('trans_cols=sapply(cols, col2trans)', file=f)
        for type in ['pdf', 'png', 'svg']:
            print(type + '("' + outprefix + '.' + type + '")', file=f)
            print('par(mar=c(10,4,4,2) + 0.1)', file=f)
            print('bar_heights1 = c(' + ','.join(bar_heights1) + ')',
                  'bar_heights2 = c(' + ','.join(bar_heights2) + ')',
                  'bar_heights = t(as.matrix(data.frame(bar_heights1, bar_heights2)))', sep='\n', file=f)
            print('barplot(bar_heights, names.arg=', names, ', '
                  ' ylab="', y_label, '", ',
                  #'col=c(cols, col2trans) ',
                  'col=c("black", "gray"), ',
                  'main="', main, '",',
                  'las=2)', sep='', file=f)

            print('dev.off()', file=f)

        utils.close(f)
        run_r_script(r_script)

    def barplot_beside(self, stat1, stat2, y_label, outprefix, main=''):
        r_script = outprefix + '.R'
        f = utils.open_file_write(r_script)
        r_syms_v, r_syms_l, r_cols_v, r_cols_l = self.get_r_vectors()

        bar_heights1 = [self.results[scaff][stat1] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        bar_heights2 = [self.results[scaff][stat2] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        #names = 'c(' + ','.join(['"' + x + '"' for x in scaffolders]) + ')'

        #all_data = [(max(int(bar_heights1[i]), int(bar_heights2[i])), bar_heights1[i], bar_heights2[i], scaffolders[i]) for i in range(len(bar_heights1))]
        all_data = [(int(bar_heights1[i]), bar_heights1[i], bar_heights2[i], scaffolders[i]) for i in range(len(bar_heights1)) if (self.data_type, scaffolders[i]) not in bad_runs]
        all_data.sort()
        bar_heights1 = [all_data[i][1] for i in range(len(all_data))]
        bar_heights2 = [all_data[i][2] for i in range(len(all_data))]
        names = 'c(' + ','.join(['"' + all_data[i][3] + '"' for i in range(len(all_data))]) + ')'
        #cols = [r_colours[x] for x in names]

        print(r''' col2trans = function(colour) {
    x=as.vector(col2rgb(colour))
    return(rgb(x[1],x[2],x[3],20,maxColorValue=255))
}
''', file=f)
        print('cols=', r_cols_v, file=f)
        print('trans_cols=sapply(cols, col2trans)', file=f)

        for type in ['pdf', 'png', 'svg']:
            print(type + '("' + outprefix + '.' + type + '")', file=f)
            print('par(mar=c(10,4,4,2) + 0.1)', file=f)
            print('bar_heights1 = c(' + ','.join(bar_heights1) + ')',
                  'bar_heights2 = c(' + ','.join(bar_heights2) + ')',
                  'bar_heights = t(as.matrix(data.frame(bar_heights1, bar_heights2)))', sep='\n', file=f)
            print('barplot(bar_heights, names.arg=', names, ', '
                  ' ylab="', y_label, '", ',
                  #'col=c(cols, col2trans) ',
                  'col=c("black", "gray"), ',
                  'beside=T,',
                  'main="', main,  '",',
                  'las=2)', sep='', file=f)

            print('dev.off()', file=f)

        utils.close(f)
        run_r_script(r_script)




    def barplot_of_one_stat_coloured(self, stat, outprefix, main=''):
        r_script = outprefix + '.R'
        f = utils.open_file_write(r_script)
        r_syms_v, r_syms_l, r_cols_v, r_cols_l = self.get_r_vectors()

        bar_heights = [self.results[scaff][stat] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        names = 'c(' + ','.join(['"' + x + '"' for x in scaffolders if (self.data_type, x) not in bad_runs]) + ')'

        for type in ['pdf', 'png', 'svg']:
            print(type + '("' + outprefix + '.' + type + '")', file=f)

            print('barplot(c(' + ','.join(bar_heights), '), ',
                  'names.arg=', names, ', '
                  'main="', main, '",',
                  ' ylab="', stat, '", ',
                  'col=', r_cols_v, ', ',
                  ')', sep='', file=f)

            print('dev.off()', file=f)

        utils.close(f)
        run_r_script(r_script)

    def barplot_of_one_stat_sorted(self, stat, outprefix, main='', stat2=None):
        r_script = outprefix + '.R'
        f = utils.open_file_write(r_script)

        if stat2 is None:
            bar_heights = [self.results[scaff][stat] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]
        else:
            bar_heights = [self.results[scaff][stat] + self.results[scaff][stat2] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]

        all_data = list(zip([int(x) for x in bar_heights], [scaff for scaff in scaffolders if (self.data_type, scaff) not in bad_runs], [r_colours[scaff] for scaff in scaffolders if (self.data_type, scaff) not in bad_runs]))


        all_data.sort()
        bar_heights = [str(x[0]) for x in all_data]
        names = 'c(' + ','.join(['"' + x[1] + '"' for x in all_data]) + ')'
        cols = ['"' + t[2] + '"' for t in all_data]

        for type in ['pdf', 'png', 'svg']:
            print(type + '("' + outprefix + '.' + type + '")', file=f)
            print('par(mar=c(10,4,4,2) + 0.1)', file=f)

            print('barplot(c(' + ','.join(bar_heights), '), ',
                  'names.arg=', names, ', '
                  'main="', main, '",',
                  ' ylab="', stat, '", ',
                  'col=c(' + ','.join(cols) + '), ',
                  'las=2',
                  ')', sep='', file=f)

            print('dev.off()', file=f)

        utils.close(f)
        run_r_script(r_script)

# dump input file into list
f = utils.open_file_read(options.infile)
infile_contents = f.readlines()
utils.close(f)

infile_header = infile_contents[0]
infile_by_test = {}
results = {} # test data -> results for that test
test_data_types = []

# make a scaffold_result object for each of the datasets
for x in infile_contents[1:]:
    test_data = x.split('\t')[0]
    if test_data not in infile_by_test:
        infile_by_test[test_data] = []
        test_data_types.append(test_data)

    infile_by_test[test_data].append(x)


for test_data in test_data_types:
    results[test_data] = Scaffold_result(infile_header, infile_by_test[test_data])

texfile = options.outprefix + '.tex'
f_tex = utils.open_file_write(texfile)

print(r'''\documentclass[11pt, a4paper]{article}
\usepackage[dvipdfm, left=1.5cm, right=1.5cm, top=1.5cm, bottom=1.5cm]{geometry}
\usepackage{graphicx}
\usepackage{grffile}
\begin{document}
''', file=f_tex)

datasets = [
    'Staph.perfect.3kb.short',
    'Staph.perfect.3kb.long',
    'Staph.perfect.10kb.short',
    'Staph.perfect.10kb.long',
    'Staph.velvet',
    'Rhodobacter',
    '3D7.short',
    '3D7.long',
    '3D7.combi',
    'human.short',
    'human.long',
    'human.combi'
]
# make some plots
for d in datasets:
    res = results[d]
    outprefix = options.outprefix + '.' + res.data_type.replace(' ', '_')
    res.plot_scatter('Correct joins', 'Bad joins', outprefix + '.good_v_bad_scatter', main='Good vs bad joins')
    #res.plot_scatter_with_chulls('Correct joins', 'Bad joins', outprefix + '.good_v_bad_scatter_with_chulls', main='Good vs bad joins')
    res.plot_scatter_with_chulls('Correct joins', 'Bad joins', outprefix + '.good_v_bad_scatter_with_chulls')
    res.plot_scatter('Correct joins', 'Bad joins', outprefix + '.good_v_bad_scatter.with_legend', legend=True, main='Good vs bad joins')
    #res.plot_scatter_with_chulls('FDR', 'TPR', outprefix + '.roc_with_chulls', main='ROC plot', legend=True)

    res.barplot_of_one_stat_coloured('% correct joins', outprefix + '.percent_good_joins_barchart', main='Proportion of correct joins')
    #res.barplot_of_one_stat_coloured('Skipped tags', outprefix + '.skipped_tags_barchart', main='Skipped tags')
    #res.barplot_of_one_stat_sorted('Skipped tags', outprefix + '.skipped_tags_barchart', main='Skipped tags')
    res.barplot_of_one_stat_sorted('Skipped tags', outprefix + '.skipped_tags_barchart')
    #res.barplot_of_one_stat_sorted('Lost tags', outprefix + '.lost_tags_barchart', main='Lost tags')
    res.barplot_of_one_stat_sorted('Lost tags', outprefix + '.lost_tags_barchart')
    #res.barplot_of_one_stat_coloured('Lost tags', outprefix + '.lost_tags_barchart', main='Lost tags')
    #res.barplot_of_one_stat_sorted('Total CPU', outprefix + '.cpu', main='CPU usage')
    res.barplot_of_one_stat_sorted('Total CPU', outprefix + '.cpu')
    #res.barplot_two_stacked('CPU', 'Extra CPU', 'CPU (s)', outprefix + '.cpu', main='CPU usage')
    #res.barplot_beside('Memory', 'Extra Memory', 'Memory (MB)', outprefix + '.mem', main='Memory usage')
    res.barplot_beside('Memory', 'Extra Memory', 'Peak Memory (MB)', outprefix + '.mem')

    # update tex file
    print(r'''\newpage''', file=f_tex)
    print(r'''\section*{''', res.data_type, '}', sep='', file=f_tex)

    def insert_pdf(fname):
        return r'''\includegraphics[width=8cm]{''' + fname + '}'

    print(insert_pdf('gather_all.boxplot.' + res.data_type + '.R.pdf'), file=f_tex)
    print(insert_pdf(outprefix + '.good_v_bad_scatter_with_chulls.pdf'), file=f_tex)
    #print(insert_pdf(outprefix + '.roc_with_chulls.pdf'), file=f_tex)
    print('', file=f_tex)
    print(r'''\noindent''', file=f_tex)
    print(insert_pdf(outprefix + '.skipped_tags_barchart.pdf'), file=f_tex)
    print(insert_pdf(outprefix + '.lost_tags_barchart.pdf'), file=f_tex)
    print('', file=f_tex)
    print(r'''\noindent''', file=f_tex)
    print(insert_pdf(outprefix + '.cpu.pdf'), file=f_tex)
    print(insert_pdf(outprefix + '.mem.pdf'), file=f_tex)
    #print(insert_pdf(outprefix + '.percent_good_joins_barchart.pdf'), file=f_tex)


print(r'''\end{document}''', file=f_tex)
utils.close(f_tex)
if not options.noplots:
    utils.syscall('pdflatex ' + texfile)
    utils.syscall('pdflatex ' + texfile)


#def get_data_by_scaffolder(data_type):
#    d = {s: [] for s in scaffolders}
#
#    for test_type in test_data_types:
#        for scaff in scaffolders:
#            d[scaff].append(results[test_type].results[scaff][data_type])
#
#    return d
#
#percent_good_by_scaff = get_data_by_scaffolder('% correct joins')


def barplot_by_scaffolder(stat_to_plot, prefix, plot_width=21, plot_height=7):
    bar_heights = []
    colours = []

    for scaff in scaffolders:
        for t in test_data_types:
            bar_heights.append(results[t].results[scaff][stat_to_plot])
            colours.append(r_colours[scaff])

    r_script = prefix + '.R'
    f = utils.open_file_write(r_script)

    for type in ['png', 'pdf']:
        print(type + '("' + outprefix + '.' + type + '", width=', plot_width, ', height=', plot_height, ')', file=f)

        print('barplot(c(' + ','.join(bar_heights), '), ',
                      #'names.arg=', names, ', '
                      ' ylab="', stat_to_plot, '", ',
                      'col=c(', ','.join(['"' + x + '"' for x in colours]), ') ',
                      ')', sep='', file=f)

        print('dev.off()', file=f)

    utils.close(f)
    run_r_script(r_script)


def barplot_by_input_data(stat_to_plot, prefix, plot_width=21, plot_height=7):
    bar_heights = []
    colours = []

    for t in test_data_types:
        for scaff in scaffolders:
            bar_heights.append(results[t].results[scaff][stat_to_plot])
            colours.append(r_colours[scaff])

    r_script = prefix + '.R'
    f = utils.open_file_write(r_script)

    for type in ['png', 'pdf']:
        print(type + '("' + outprefix + '.' + type + '", width=', plot_width, ', height=', plot_height, ')', file=f)

        print('barplot(c(' + ','.join(bar_heights), '), ',
                      #'names.arg=', names, ', '
                      ' ylab="', stat_to_plot, '", ',
                      'col=c(', ','.join(['"' + x + '"' for x in colours]), ') ',
                      ')', sep='', file=f)

        print('dev.off()', file=f)

    utils.close(f)
    run_r_script(r_script)


for x in [('% correct joins', 'percent_correct_joins'),
          ('Lost tags', 'lost_tags'),
          ('Skipped tags', 'skipped_tags')]:
    pass
    barplot_by_scaffolder(x[0], options.outprefix + '.' + x[1] + '_by_scaff')
    barplot_by_input_data(x[0], options.outprefix + '.' + x[1] + '_by_data')



