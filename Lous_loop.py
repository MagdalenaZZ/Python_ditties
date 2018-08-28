#!/usr/bin/env python

"""

Script to help Lou's group do custom annotation of neuroblastoma and medulloblastoma datasets
magdalena.z@icr.ac.uk

Written: Feb 2018



"""



import sys
import argparse
import re
import os.path
from subprocess import call

epi = ('\
    \n\
    This script is to help do annotation of the factors most important in neuroblastoma and medulloblastoma\n\
    \n\
    The input is a folder with prepared gene lists of the type:\n\
    ENSEMBL_ID <tab> annotation\n\
    \n\
    And a dataset, expression dataset or log2-fold changes\n\
    \n\
    For log2fold changes, and GSEA gene enrichment calculation is performed on the custom gene lists\n\
    For expression values, a heatmap is created of those particular gene lists\n\
    For a gene list, overlaps are calculated\n\
    \n\
')


#print(epi)

# Describe what the script does
parser = argparse.ArgumentParser(description='This script annotates output', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-f', '--folder', default='/data/rds/shared/PEDBIO/REFERENCES/NB_annotation', dest='fol', action='store', required=False, help="Folder containing gene lists")
parser.add_argument('-i', '--input', default=None, dest='inf', action='store', required=True, help="Input data file to be annotated")
parser.add_argument('-t', '--type', default=None, dest='t', action='store', required=True, help="Type of input data: EXP/L2F/GENES")
parser.add_argument('-a', '--alias', default='/data/rds/shared/PEDBIO/REFERENCES/alias.txt', dest='al', action='store', required=False, help="Translation of unique identifiers to original identifiers")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()
#print(args.rank)

# Check if input files exist
if not os.path.isdir(args.fol)==True:
    print("Cannot find input directory ",args.fol)
    sys.exit(1)

if not os.path.isfile(args.inf)==True:
    print("Cannot find input file ",args.inf)
    sys.exit(1)

chip ='/home/breakthr/mzarowiecki/REF/GSEA/GENE_SYMBOL.chip'
if not os.path.isfile(chip)==True:
    print("Cannot find gene symbol chip file ",chip)
    sys.exit(1)


# Create dictionary of aliases for unique identifiers



# Predict which species it is

import subprocess

species = subprocess.check_output( ' '.join(("head"," -2",args.inf,"| tail -1 | cut -f1")))

print("Species ", species)

if args.t=="EXP" or args.t=="E":
    print("You have input an expression matrix")

# Now, call a script that can handle that expression matrix and create a heatmap from it

# For each file in this folder

# that matches right species


# Input is log2fold changes
elif args.t=="L2F" or args.t=="L":
    print("You have input a log2fold list with: gene <tab> l2f")

# Input is a list of genes
elif args.t=="GENES" or args.t=="G":
    print("You have input a gene list")

# Input is what?
else:
    print("I don't understand your input ",args.t)
    print("Please choose EXP;expression, L2F;log2fold changes or GENES; list of genes")
    sys.exit(1)









"""
print (''.join(("\njava -Xmx2048m  -cp /Users/mzarowiecki/bin/gsea-3.0_beta_3.jar xtools.gsea.GseaPreranked -gmx ",args.gs," -collapse false -mode Max_probe -norm meandiv -nperm ",args.iter, " -rnk ",args.rank," -scoring_scheme weighted -rpt_label my_analysis -chip ",chip," -include_only_symbols true -make_sets true -plot_top_x 4726 -rnd_seed timestamp -set_max 5000 -set_min 4 -zip_report false -out ",args.rank,"_OUT -gui false > ",args.rank,".err \n\n")))

print ("Starting gsea\n")
call (''.join(("java -Xmx2048m  -cp /Users/mzarowiecki/bin/gsea-3.0_beta_3.jar xtools.gsea.GseaPreranked -gmx ",args.gs," -collapse false -mode Max_probe -norm meandiv -nperm ",args.iter, " -rnk ",args.rank," -scoring_scheme weighted -rpt_label my_analysis -chip ",chip," -include_only_symbols true -make_sets true -plot_top_x 4726 -rnd_seed timestamp -set_max 5000 -set_min 4 -zip_report false -out ",args.rank,"_OUT -gui false > ",args.rank,".err")))
print ("Finished gsea\n")








of = args.rank + ".R"



orig_stdout = sys.stdout
f = open('of', 'w')
sys.stdout = f

for i in range(2):
    print( 'i = ', i)

sys.stdout = orig_stdout
f.close()


"""




    
