#!/usr/bin/python

import os.path
import subprocess
from collections import defaultdict

'''
Functions of general use for NGS pipelines
'''

def read_samples(sample_file):
	'''
	Read a tsv file containing the sample name and either one or two
	fastq files to align to a reference genome
	'''
	samples = defaultdict(list)
	with open(sample_file, "rU") as f:
		for line in f:
			line = line.rstrip("\n\r")
			a = line.split("\t")
			if len(a) == 2:
				samples[a[0]].append(a[1])
			else:
				samples[a[0]].append(a[1])
				samples[a[0]].append(a[2])
	return samples

def ensure_dir(f):
	'''
	Check if a directory exists and if not, create it
	'''
	d = os.path.dirname(f)
	if not os.path.exists(d):
		os.makedirs(d)







