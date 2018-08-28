#!/usr/bin/env python
from __future__ import print_function
import hashlib
import sys

""" 

This script 

"""

if (len(sys.argv)<2):
    print ("Usage: md5_do.py   <multirunID>      ")
    print ("Example: md5_do.py LP3000944-DNA_D09 LP3000947-DNA_F06  ")
    print("Returns the name of the multisample folder these samples are in ")
    quit(1)

i = sys.argv[1:]

#print(i)
#ids=("hi","there")
#ids=("LP3000944-DNA_D09","LP3000947-DNA_F06")

md5_hash = hashlib.md5('-'.join([str(x) for x in sorted(i)])).hexdigest()

print ("/genomes/analysis/multisample/",md5_hash, sep='')







