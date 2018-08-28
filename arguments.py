#!/usr/bin/env python


"""

Nice descriptive header

"""


import sys
import argparse


# Make sure you are working on the right version of Python
if sys.version_info[0] == 3:
    print ('\nCAVA does not run on Python 3.\n')
    quit()



# Command line argument parsing


parser = argparse.ArgumentParser(description='Process some integers.')

descr = 'OpEx (Optimised Exome) pipeline ' + ver + '.'
parser = OptionParser(usage='python opex.py <options>', version=ver, description=descr)
parser.add_option('-i', "--input", default=None, dest='fastq', action='store', help="fastq.gz files")
parser.add_option('-o', "--output", default=None, dest='name', action='store', help="Sample name (output prefix)")
parser.add_option('-b', "--bed", default=None, dest='bed', action='store', help="Bed file")
parser.add_option('-r', "--reference", default=None, dest='reference', action='store', help="Reference genome file")
parser.add_option('-t', "--threads", default=1, dest='threads', action='store', help="Number of processes to use")
parser.add_option('-f', "--full", default=False, dest='full', action='store_true',help="Output full CoverView output [default value: %default]")
parser.add_option('-c', "--config", default=None, dest='config', action='store', help="Configuration file")
parser.add_option('-k', "--keep", default=False, dest='keep', action='store_true', help="Keep temporary files")
parser.add_option('-h', "--help", default=None, help="This is a help message")

(options, args) = parser.parse_args()
checkInputs(options)

args = parser.parse_args()

# complain if something is missing
if not 'REFERENCE' in params.keys():
    if options.reference is None:
        print ('Error: no reference genom provided.')
        
        quit()






