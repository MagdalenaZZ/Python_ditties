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





######### Search for a list of words in a dictionary ##############

def searchwords(words, sen):
    """This sub takes a pmid eg 30342085 and returns the fulltext paper in sentences"""
    
    for reg in words:

        #print (reg, words[reg])

        if re.search( words[reg], sen):
            print (reg, sen)




'''
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
