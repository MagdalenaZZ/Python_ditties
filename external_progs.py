import os
import utils

class Error (Exception): pass

bin = os.path.join(os.path.expanduser('~mh12'), 'bin')

bowtie2_bin = os.path.join(bin, 'bowtie2-2.0.5')
bowtie2_align = os.path.join(bowtie2_bin, 'bowtie2-align')
bowtie2_build = os.path.join(bowtie2_bin, 'bowtie2-build')
bowtie2_index_file_extensions = ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']

bwa = os.path.join(bin, 'bwa-0.6.2', 'bwa')

mummer_dir = os.path.join(bin, 'MUMmer3.23')
nucmer = os.path.join(mummer_dir, 'nucmer')
delta_filter = os.path.join(mummer_dir, 'delta-filter')
show_coords = os.path.join(mummer_dir, 'show-coords')

samtools = os.path.join(bin, 'samtools-0.1.18', 'samtools')

smalt = os.path.join(bin, 'smalt-0.7.0.1', 'smalt_x86_64')


def bwa_index_clean(prefix):
    for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']:
        os.unlink(prefix + '.' + ext)


def is_bowtie2_indexed(file):
    for f in [file + '.' + x for x in bowtie2_index_file_extensions]:
        if not os.path.exists(f):
            return False

    return True


def index_with_bowtie2(file):
    if not is_bowtie2_indexed(file):
        utils.syscall(bowtie2_build + ' ' + file + ' ' + file)


def clean_bowtie2_index(file):
    for f in [file + '.' + x for x in bowtie2_index_file_extensions]:
        os.unlink(f)

