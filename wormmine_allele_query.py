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

# The view specifies the output columns
query.add_view(
    "linkedTo", "otherName", "species", "productionMethod", "sequenceStatus",
    "status", "type", "KOConsortiumAllele", "primaryIdentifier", "symbol",
    "typeOfMutation"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Allele.linkedTo", "ASC")


f = open('Alleles.txt', 'w')

for row in query.rows():
    print(row["linkedTo"], row["otherName"], row["species"], row["productionMethod"], \
        row["sequenceStatus"], row["status"], row["type"], row["KOConsortiumAllele"], \
        row["primaryIdentifier"], row["symbol"], row["typeOfMutation"], file=f)



f.close()



from __future__ import print_function
from intermine.webservice import Service
service = Service("http://intermine.wormbase.org/tools/wormmine/service")
query = service.new_query("Publication")
query.add_view("pubMedId", "firstAuthor", "title", "year", "journal", "volume", "pages")

for row in query.rows():
    print(row["pubMedId"], row["firstAuthor"], row["title"], row["year"], row["journal"], row["volume"], row["pages"], file=f)




