#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("http://intermine.wormbase.org/tools/wormmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Allele")

query.add_constraint("primaryIdentifier", "=", "WBVar00604270", code="A")

query.add_view(
    "amber_UAG", "detectionMethod", "missenseChange", "linkedTo", "ochre_UAA",
    "opal_UGA", "typeOfMutation", "gene.symbol", "gene.CDSs.primaryIdentifier",
    "gene.CDSs.symbol", "gene.CDSs.protein.CDSs.primaryIdentifier",
    "gene.CDSs.protein.CDSs.symbol"
)

# Just shows the query
#print (query)

for row in query.rows():
    print (row)



quit()

query.add_view(
    "amber_UAG", "detectionMethod", "missenseChange", "linkedTo", "ochre_UAA",
    "opal_UGA", "typeOfMutation", "gene.symbol", "gene.CDSs.primaryIdentifier",
    "gene.CDSs.symbol", "gene.CDSs.protein.CDSs.primaryIdentifier",
    "gene.CDSs.protein.CDSs.symbol"
)
query.add_constraint("Allele", "LOOKUP", "WBVar00604270", "C. elegans", code="A")

for row in query.rows():
    print(row["amber_UAG"], row["detectionMethod"], row["missenseChange"], row["linkedTo"], \
        row["ochre_UAA"], row["opal_UGA"], row["typeOfMutation"], row["gene.symbol"], \
        row["gene.CDSs.primaryIdentifier"], row["gene.CDSs.symbol"], \
        row["gene.CDSs.protein.CDSs.primaryIdentifier"], row["gene.CDSs.protein.CDSs.symbol"])



Y74C9A.3 (I216T)"

# The view specifies the output columns
query.add_view("gene.symbol", "gene.allele.symbol", "typeOfMutation", "linkedTo")

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Allele.gene.symbol", "ASC")

for row in query.rows():
    print(row["gene.symbol"], row["gene.allele.symbol"], row["typeOfMutation"], row["linkedTo"])

