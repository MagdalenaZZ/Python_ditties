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
import nltk.data
import os
from xml.dom import minidom
from bs4 import BeautifulSoup
import requests


############ Convert WBPaper to pmid ########
''' Subroutine to convert a WBPaper to pmid '''

def wbp2pmid (pmid):
    # First try to get it from existing reference
    with open('WpaXref', 'rt') as csvfile:
        reader = csv.reader(csvfile, delimiter ='\t')
        for row in reader:
            print ("Row",row[0])
            wbp = row[0].split('')

    # Then try to refresh current if it is a few days old
    url = 'http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref' 
    r = requests.get(url)
    y=BeautifulSoup(r.content)
    return y.record["pmcid"]




############ Convert pmid to pmcid ########
''' Subroutine to convert a pmid to pmcid '''

def pmid2pmcid (pmid):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=' + pmid 
    r = requests.get(url)
    y=BeautifulSoup(r.content)
    return y.record["pmcid"]



############ Convert pmid to pmcid ########
''' Subroutine to convert a pmid to doi '''

def pmid2doi (pmid):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=' + pmid 
    r = requests.get(url)
    y=BeautifulSoup(r.content)
    return y.record["doi"]




'''
https://cgc.umn.edu/static/cgc-strains.txt '''


'''



https://doi.org/10.17912/kedf-yn42


curl 'http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/paper_display.cgi?action=Search+!&data_number=00054673'   # this is how to get a doi
curl 'http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/paper_display.cgi?action=Search+!&data_number=00050739'   # this is how to get a doi


# Get dois for ALL WBPapers
curl 'http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref' 



https://www.ncbi.nlm.nih.gov/pmc/tools/get-full-text/
https://www.ncbi.nlm.nih.gov/pmc/tools/oai/

This E-Summary call retrieves article metadata for a PMCID:
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pmc&id=3539452&retmode=json

The returned data includes the PubMed ID (PMID) for this article, so you can do a further call to get PubMed metadata for the same article:
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=22368089&retmode=json.

You can use E-Fetch with PubMed to get even more:
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=22368089

Use the OAI-PMH service to get metadata in Dublin Core format:
    https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:3539452&metadataPrefix=oai_dc

Or in PMC front-matter format:
    https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:3539452&metadataPrefix=pmc_fm



# Id converter between pubmed, pmcid and doi
https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=10848627,10848628 # webpage 2 pmids to the others
https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=PMC3283037.2&format=json # json format link from PMC to doi and PMID
https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=3283037&idtype=pmcid&tool=my_tool&format=json # forcing PMC id on the ID format
# Root tool
https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/
# API docu
https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/

# Entrez/pubmed REST documentation
https://www.ncbi.nlm.nih.gov/books/NBK25500/

# Pubmed abstract and json format search terms and PMCID
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=31014992
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=33092268 # under ids; other; tag str is PMC7589727

# OAI XML of fulltext paper im PMC example
https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:4304705&metadataPrefix=pmc
https://www.ncbi.nlm.nih.gov/pmc/tools/oai/


# PMC REST service documentation
https://europepmc.org/RestfulWebService
https://europepmc.org/developers







'''
