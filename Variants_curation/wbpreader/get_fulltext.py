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


######### Try to retrieve the full-text document from WB ##############

def fulltext_wbp(wbpid,repo):
    """This sub takes a wbpid eg WBPaper00056731 and returns the fulltext paper in sentences"""
    
    ft=[0];

    # Check that wbpid is a valid WBPaper
    if not re.match( 'WBPaper', wbpid):
        print (wbpid, "is not a valid WBPaper ID")
        return ft


    # Download paper if it doesn't exist
    fn = repo + '/' + wbpid + '.json'

    if os.path.exists(fn) and os.path.getsize(fn) > 16:
        pass
        #print (fn,"exists already, so skipping")
    else:
        #print (fn,"does not exist, so attempting to download")
        #com1 = '-k '+ '\n'+'-d "{"token":"DYpds7m8El78T9n8qBKW", "query": {"accession": "' + wbpid +'", "type": "document", "corpora": ["C. elegans"]}, "include_fulltext": true}"'
        com1 = '-o '+fn +'\n-k '+ '\n'+'-d "{\\"token\\":\\"DYpds7m8El78T9n8qBKW\\", \\"query\\": {\\"accession\\": \\"' + wbpid +'\\", \\"type\\": \\"document\\", \\"corpora\\": [\\"C. elegans\\"]}, \\"include_fulltext\\": true}"'

        configf= wbpid + '.tmp.config'
        curlf = open(configf,'w')
        print (com1, file=curlf)
        curlf.close()
        command = 'curl -o '+ fn +' -K '+ configf+' https://textpressocentral.org:18080/v1/textpresso/api/search_documents' 
        comlist = command.split()
        #print (command)
        os.system(command)
        #proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #(out, err) = proc.communicate()
        #print ("OUT", out, "ERR", err)
        #print ()
        #OF= open(fn, "wb")
        #OF.write(out)
        #OF.close()

    # Remove the % sign at the end of the file

    '''with open(fn, 'rb+') as fh:
        fh.seek(-1, os.SEEK_END)
        fh.truncate()
        print ("Truncated", fh)
    '''

    # Read the paper, and split into sentences

    if os.path.exists(fn) and os.path.getsize(fn) > 20:

        # Open our JSON file and load it into python
        input_file = open (fn)
        json_array = json.load(input_file)


        for item in json_array:
            abs = item["abstract"]
            fullt =  item["fulltext"]
            tokenizer = nltk.data.load('tokenizers/punkt/english.pickle')
            #fp = open("test.txt")
            #data = fp.read()
            ft = tokenizer.tokenize(abs)
            ftt=tokenizer.tokenize(fullt)
            ft = ft +ftt

        #print ("FT", ft)

    else:
        print ("File", fn, "does not contain data")

    return ft



######### Try to retrieve the full-text document from PMID ##############

def fulltext_pmid(pmid,repo):
    """This sub takes a pmid eg 30342085 and returns the fulltext paper in sentences"""
    
    ft=[0];

    # Check that wbpid is a valid WBPaper
    if not re.match( '\d+', pmid):
        print (pmid, "is not a valid pubmedID")
        return ft


    # Download paper if it doesn't exist
    fn = repo + '/' + pmid + '.pmid' + '.json'
    
    # Paper already exists, so load
    if os.path.exists(fn) and os.path.getsize(fn) > 0:
        print (fn,"exists already, so skipping")
        with open(fn) as json_file:
            #y = json.load(json_file)
            y = json.load(json_file)
    # Paper does not exist, so attempt to download 
    else:
        #print (fn,"does not exist, so attempting to download")

        # Try to retrieve the full-text document from PMC, and save it as a json file
        #url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=' + pmid 
        url= 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=' + pmid
        url = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/'+ pmid + '/unicode'
        r = requests.get(url)

        # Save json
        curlf=open(fn, 'wb')
        curlf.write(r.content)
        curlf.close()

        y = json.loads(r.content)
    
    tokenizer = nltk.data.load('tokenizers/punkt/english.pickle')

    for docs in y["documents"]:
        for infons in docs["passages"]:
            ftt = tokenizer.tokenize(infons["text"])
            for item in ftt:
                ft.append(str(item))

    return ft


######### Try to retrieve the full-text document from PMCID ##############

def fulltext_pmcid(pmcid,repo):
    """This sub takes a pmcid eg PMC7252405 and returns the fulltext paper in sentences"""
    
    ft=[0];

    # Check that pmcid is a valid PMCID
    if not re.match( 'PMC\d+', pmcid):
        print (pmcid, "is not a valid PMCID")
        return ft
    else:
        print (pmcid, "is a valid PMCID")






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
