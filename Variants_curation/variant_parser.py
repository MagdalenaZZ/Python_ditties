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
from wbpreader.readerdef import fulltext_wbp
from wbpreader.readerdef import fulltext_pmid
#print ("Finished importing subs")

'''

Script to parse the variants file and highlight sentences


'''

epi = ('\
    \n\
	Give an input file with variants IDs\n\
    Your genome.intervals_list is output from GATK pipeline, and a panel-of-normals (from PON.py) \n\
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
'''

# Retrieve input file from URL
url = args.inp
r = requests.get(url, allow_redirects=True)
#fo = open("foo.txt", "wb")
open(args.out, 'wb').write(r.content)
#print (r.content, url)
#res.close()

# Check if output files exist
if not os.path.isfile(args.out)==True:
    print("Cannot find output file ",args.out)
    sys.exit(1)
'''

# Open output file
outf = open(args.out,'w')
        #print (com1, file=curlf)
        #outf.close()
#print ("Trying to open", args.out)


# Read in variants file
with open(args.inp, 'rt') as csvfile:
    reader = csv.reader(csvfile, delimiter ='\t')
    for row in reader:
        #print ("Row", row)
        wbp = row[0].split('\t')
        #print (row[0], wbp[0])
        seen=dict()

        # Is a WBPapaer
        if re.match( 'WBPaper', wbp[0]):
           #print ("Match",wbp[0] )

            # Try to retrieve the full-text document from WB, and look for valid sentences
            ft = fulltext_wbp(wbp[0], args.tmp)

            # Also check pmid data, if textpresso didn't work
            if len(ft) < 10:
                if bool(re.search(r'\d', row[4])):
                    ft = fulltext_pmid(row[4], args.tmp)
                    print ("GET FT from pmid", row[4])
                else:
                    print ("No number in pmid", row[4])


            for sen in ft:
                sen=str(sen)
                sen = sen.strip('\n')
                # Check the sentence has not already been processed
                if sen in seen:
                    pass
                else:
                    seen[sen] = 1
                    #print ("SEN", str(sen))
                    #pattern = re.compile(r"exon", )
                    #pattern = re.compile(r'RNA', flags=re.IGNORECASE)
                    #genes=rf'\b{'gene'}\b' 
                    #genes = /\b(?:one|two|three)\b/gi
                    #re.compile(r"\gene\b",flags=re.IGNORECASE )
                    var= "Variation-" + row[3]
                    if re.search(r'SNP', sen,flags=re.IGNORECASE):
                        print (wbp[0], "SNP", sen, sep='\t', file=outf)
                    elif re.search( 'substitution', sen,flags=re.IGNORECASE):
                        print (wbp[0], "Substitution", sen, sep='\t', file=outf)
                    elif re.search( 'deletion', sen,flags=re.IGNORECASE):
                        print (wbp[0], "Deletion", sen, sep='\t', file=outf)
                    elif re.search( 'insertion', sen,flags=re.IGNORECASE):
                        print (wbp[0], "Insertion", sen, sep='\t', file=outf)
                    elif re.search( 'mutation', sen,flags=re.IGNORECASE):
                        print (wbp[0], "Mutation", sen, sep='\t', file=outf)
                    elif re.search( 'point mutation', sen,flags=re.IGNORECASE):
                        print (wbp[0], "PointMutation", sen, sep='\t', file=outf)
                    elif re.search( 'amino acid', sen,flags=re.IGNORECASE):
                        print (wbp[0], "Amino acid", sen, sep='\t', file=outf)
                    elif re.search( 'stop codon', sen,flags=re.IGNORECASE):
                        print (wbp[0], "Stop codon", sen, sep='\t', file=outf)
                    elif row[3]:
                        if re.search( row[3], sen, flags=re.IGNORECASE):
                            print (wbp[0], var, sen, sep='\t', file=outf)
                    else:
                        pass

            print ( sep='\t', file=outf)




            '''
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
            '''

        else:
            print ("NO ",wbp[0] , args.tmp)


outf.close()


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

