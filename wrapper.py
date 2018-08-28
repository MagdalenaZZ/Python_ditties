#!/usr/bin/env python

"""

Nice descriptive help text

"""



import sys
import argparse
import re
import os.path
from subprocess import call

epi = ('\
    \n\
    ranked-list has this structure:\n\
    \n\
    Gene_ID   <TAB>  3.21\n\
    Gene_ID   <TAB>  -3.1\n\
    Gene_ID   <TAB>  -1.6\n\
    Gene_ID   <TAB>  0\n\
    \n\
    GeneID will be translated to symbols, i.e. RSF1	ARID4A	NR6A1	TIMM50\n\
    using $rosetta, which points to popularly used conversion_files.chip\n\
    More conversion_files.chip are found in the GSEA distribution.\n\
    \n\
    Geneset is downloaded from MSigDB .gmt format\n\
    default is msigdb.v6.0.symbols.gmt\n\
    \n\
')


#print(epi)

# Describe what the script does
parser = argparse.ArgumentParser(description='This script wraps something', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-n', '--iterations', default=10, dest='iter', action='store', required=True, help="Number of GSEA iterations")
parser.add_argument('-i', '--input', default=None, dest='rank', action='store', required=True, help="Ranked gene list, tab delimited, gene name and rank")
parser.add_argument('-g', '--geneset', default=None, dest='gs', action='store', required=True, help="Gene set database .gmt")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()
#print(args.rank)

# Check if input files exist
if not os.path.isfile(args.rank)==True:
    print("Cannot find input file ",args.rank)
    sys.exit(1)

if not os.path.isfile(args.gs)==True:
    print("Cannot find input file ",args.gs)
    sys.exit(1)

chip ='/Users/mzarowiecki/bin/GSEA/GENE_SYMBOL.chip'
if not os.path.isfile(chip)==True:
    print("Cannot find gene symbol chip file ",chip)
    sys.exit(1)



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







    
