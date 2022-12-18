#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
import subprocess
import requests
import csv
import re
import json
#sys.path.append('/Users/mz3/bin/Python_ditties/Variants_curation')
from wbpreader.get_fulltext import fulltext_wbp
from wbpreader.id_conversions import wbp2pmid
#from xml.dom import minidom
#from html.parser import HTMLParser
from bs4 import BeautifulSoup

#print ("Finished importing subs")

'''

Script to parse the NN file and highlight sentences


'''

epi = ('\
    \n\
	Script to parse the NN file and highlight sentences\n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes the commands for running GATK CNV pipeline from tumour and normal data', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='inp', action='store', required=True, help="http://svm.textpresso.org/celegans/svm_results/20201206/120620_091120_structcorr")
parser.add_argument('-t', '--tempdir', default='Tmpdir', dest='tmp', action='store', required=True, help="Directory for temporary results")
parser.add_argument('-o', '--output', default=None, dest='out', action='store', required=True, help="output file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

if not os.path.exists(args.tmp):
    os.makedirs(args.tmp)
    print ("Created folder", args.tmp)

# Check if input files exist
if not os.path.isfile(args.out)==True:
    print("Cannot find output file ",args.out)
    sys.exit(1)

# Retrieve input file from URL
url = args.inp
r = requests.get(url)
y=BeautifulSoup(r.content, features="html.parser")
#print (y.body.a.text)
#print (y.find("a", href=True)["href"])


# Read in NN file and parse each paper
for tag in y.find_all(re.compile("^a")):
    #print(tag.text, tag["href"])

        # Is a WBPapaer
        if re.match( 'WBPaper', tag.text):
            #print ("Match",tag.text )
            # Try to retrieve the full-text document from WB, and look for valid sentences
            ft = fulltext_wbp(tag.text, args.tmp)

        else:
            print ("ERROR: The input url page is misformed, should not happen", tag.text )



#print ("Trying to open", args.out)

#exit;


# Read in NN file
#mydoc = minidom.parse(args.out)
#url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=' + pmid 
#r = requests.get(url)
#y=BeautifulSoup(r.content)


quit()

'''
with open(args.out, 'rt') as csvfile:
    reader = csv.reader(csvfile, delimiter ='\t')
    for row in reader:
        #print ("Row",row[0])
        wbp = row[0].split('')
        #print (row[0], wbp[0])

        # Is a WBPapaer
        if re.match( 'WBPaper', wbp[0]):
            print ("Match",wbp[0] )

            # Try to retrieve the full-text document from WB, and look for valid sentences
            ft = fulltext_wbp(wbp[0], args.tmp)
            for sen in ft:
                sen =str(sen)
                #print ("SEN", str(sen))
                #pattern = re.compile(r"exon", )
                #pattern = re.compile(r'RNA', flags=re.IGNORECASE)
                #genes=rf'\b{'gene'}\b' 
                #genes = /\b(?:one|two|three)\b/gi
                #re.compile(r"\gene\b",flags=re.IGNORECASE )
                if re.search(r'exon', sen,flags=re.IGNORECASE):
                    print (wbp[0], "Exon", sen)
                elif re.search( 'intron', sen, flags=re.IGNORECASE):
                    print (wbp[0], "Intron", sen)
                elif re.search( 'promoter', sen,flags=re.IGNORECASE):
                    print (wbp[0], "Promoter", sen)
                elif re.search( gene, sen):
                    print (wbp[0], "Gene", sen)
                elif re.search( 'splice', sen,flags=re.IGNORECASE):
                    print (wbp[0], "Splice", sen)
                elif re.search( 'replicate', sen,flags=re.IGNORECASE):
                    print (wbp[0], "Replicate", sen)

                else:
                    pass

            ft = fulltext_pmid(row[2], args.tmp)
            print ("FT", ft)

            
            # Try to retrieve the full-text document from PMC
            if re.match( '\d+', row[2]):

                fn = args.tmp + '/' + wbp[0] + '.json'
                #command = 'curl -o ' + fn + ' https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:' + row[2] + '&metadataPrefix=pmc' 
                command = 'curl -o ' + fn + ' https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=' + row[2] 
                comlist = command.split()
                #print (row[2])
                #print (comlist)
                #proc = subprocess.Popen(comlist, stdout=subprocess.PIPE)
                #(out, err) = proc.communicate()
                                #open(fn, 'wb').write(out)
            else:
                print ("Missing pubmed ref", wbp[0], row[2])
            

        else:
            print ("NO ",wbp[0] , args.tmp)

'''



quit()



# This will be our resulting structure mapping compound ChEMBL IDs into target uniprot IDs
compounds2targets = dict()

# First, let's just parse the csv file to extract compounds ChEMBL IDs:
with open('compounds_list.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        compounds2targets[row[0]] = set()

# OK, we have our source IDs, let's process them in chunks:
chunk_size = 50
keys = compounds2targets.keys()

for i in range(0, len(keys), chunk_size):
    # we jump from compounds to targets through activities:
    activities = new_client.activity.filter(molecule_chembl_id__in=keys[i:i + chunk_size])
    # extracting target ChEMBL IDs from activities:
    for act in activities:
        compounds2targets[act['molecule_chembl_id']].add(act['target_chembl_id'])

# OK, now our dictionary maps from compound ChEMBL IDs into target ChEMBL IDs
# We would like to replace target ChEMBL IDs with uniprot IDs

for key, val in compounds2targets.items():
    # We don't know how many targets are assigned to a given compound so again it's
    # better to process targets in chunks:
    lval = list(val)
    uniprots = set()
    for i in range(0, len(val), chunk_size):
        targets = new_client.target.filter(target_chembl_id__in=lval[i:i + chunk_size])
        uniprots |= set(sum([[comp['accession'] for comp in t['target_components']] for t in targets],[]))
    compounds2targets[key] = uniprots

# Finally write it to the output csv file
with open('compounds_2_targets.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for key, val in compounds2targets.items():
        writer.writerow([key] + list(val))

