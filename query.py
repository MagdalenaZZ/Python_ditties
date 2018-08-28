#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("http://www.humanmine.org/humanmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "synonyms.value", "symbol", "name", "primaryIdentifier",
    "secondaryIdentifier", "organism.name", "briefDescription"
)

for row in query.rows():
    print row["synonyms.value"], row["symbol"], row["name"], row["primaryIdentifier"], \
        row["secondaryIdentifier"], row["organism.name"], row["briefDescription"]

