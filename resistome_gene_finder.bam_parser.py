#!/usr/bin/env python3

import subprocess
import pysam
import copy
import argparse
import fastaq

parser = argparse.ArgumentParser(
    description = 'Helper script for resistome_gene_finder.2.py',
    usage = '%(prog)s [options] <in.bam> <outprefix>')
parser.add_argument('bam_in', help='Name of input BAM file', metavar='FILENAME')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()


def sam_record_to_soft_clipped(s):
    if s.cigar is None or len(s.cigar) == 0:
        return False, False
    return (s.cigar[0][0] == 4, s.cigar[-1][0] == 4)



def write_tabix_file(d, fname):
    f = fastaq.utils.open_file_write(fname)
    for seq in d:
        lines = [(x,d[seq][x][0],d[seq][x][1]) for x in d[seq]]
        lines.sort()
        for l in lines:
            print(seq, l[0], l[1], l[2], sep='\t', file=f)

    fastaq.utils.close(f)
    fastaq.utils.syscall('bgzip ' + fname)
    fastaq.utils.syscall('tabix -s 1 -b 2 -e 2 ' + fname + '.gz')

def decode(s):
    try:
        s = s.decode()
    except:
        return s

    return s


def mpileup_to_coverage(mpileup_cmd):
    depth = {}
    mpileup_cmd += ' | cut -f 1,2,4'
    print(mpileup_cmd)
    mpileup_out = decode(subprocess.Popen(mpileup_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).communicate()[0]).split('\n')[:-1]
    for line in mpileup_out:
        seq, pos, cov = line.rstrip().split()
        if seq not in depth:
            depth[seq] = {}
        depth[seq][int(pos)] = int(cov)

    return depth
       
    

def combine_coverage(d1, d2):
    d = {}
    refs = set(list(d1.keys()) + list(d2.keys()))

    for ref in refs:
        keys = []

        if ref in d1:
            keys.extend(d1[ref].keys())

        if ref in d2:
            keys.extend(d2[ref].keys())

        keys = set(keys)
        d[ref] = {}

        for p in keys:
            d[ref][p] = [0, 0]
            if ref in d1 and p in d1[ref]:
                d[ref][p][0] = d1[ref][p]

            if ref in d2 and p in d2[ref]:
                d[ref][p][1] = d2[ref][p]

    return d
               
        

soft_clipped = {}

sam_reader = pysam.Samfile(options.bam_in, "rb")



for s in sam_reader.fetch(until_eof=True):
    if s.is_unmapped:
        continue

    ref_name = sam_reader.getrname(s.tid)
    left_clip, right_clip = sam_record_to_soft_clipped(s)

    if left_clip or right_clip:
        if ref_name not in soft_clipped:
            soft_clipped[ref_name] = {}

        if left_clip:
            p = s.pos
            if p not in soft_clipped[ref_name]:
                soft_clipped[ref_name][p] = [0, 0]
            soft_clipped[ref_name][p][0] += 1

        if right_clip:
            p = s.aend
            if p not in soft_clipped[ref_name]:
                soft_clipped[ref_name][p] = [0, 0]
            soft_clipped[ref_name][p][1] += 1




write_tabix_file(soft_clipped, options.outprefix + '.artemis_plot.soft_clipped')

# 0x0008 = mate unmapped.
# 0x0010 = strand of query
# 0x18 - mate unmapped and strand of query
unmapped_mate_depth_fwd = mpileup_to_coverage('samtools view -b -f 8 -F 20 ' + options.bam_in + ' | samtools mpileup -A - ')
unmapped_mate_depth_rev = mpileup_to_coverage('samtools view -b -f 24 -F 4 ' + options.bam_in + ' | samtools mpileup -A - ')
unmapped_mate_depth = combine_coverage(unmapped_mate_depth_fwd, unmapped_mate_depth_rev)
write_tabix_file(unmapped_mate_depth, options.outprefix + '.artemis_plot.unmapped_mate')
