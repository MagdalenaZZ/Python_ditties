#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Converts a fastq file to fasta + qual file',
    usage = '%(prog)s [options] <fastq_in> <fasta_out>')
parser.add_argument('fastq_in', help='Name of input fastq file')
parser.add_argument('fasta_out', help='Name of output fasta (fasta_out.qual will also be created)')
options = parser.parse_args()


seq_reader = fastn.file_reader(options.fastq_in)
fasta_out = utils.open_file_write(options.fasta_out)
qual_out = utils.open_file_write(options.fasta_out + '.qual')
fastn.Fasta.line_length = 0

for seq in seq_reader:
    fa, qual = seq.to_Fasta_and_qual()
    print(fa, file=fasta_out)
    print('>' + fa.id, ' '.join([str(x) for x in qual]), sep='\n', file=qual_out)


utils.close(fasta_out)
utils.close(qual_out)

