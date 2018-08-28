#!/usr/bin/env python3

import sys
import os
import fastaq




def rename_seqs(infile, outfile, prefix):
    s = fastaq.sequences.file_reader(infile)
    f = fastaq.utils.open_file_write(outfile)
    counter = 1

    for seq in s:
        seq.id = prefix + '.' + str(counter)
        counter += 1
        print(seq, file=f)

    fastaq.utils.close(f)

scaffolds_files = []


dirs = os.walk(os.getcwd()).__next__()[1]

for d in dirs:
    old_fname = os.path.join(d, 'scaffolds.final.gapfilled.fa')
    new_fname = os.path.join(d, 'scaffolds.final.gapfilled.renamed.fa')
    if os.path.exists(old_fname):
        rename_seqs(old_fname, new_fname, d) 
        scaffolds_files.append(new_fname)
    



    
