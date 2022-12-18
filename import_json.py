
# Load up JSON Function
import json

# Open our JSON file and load it into python
input_file = open ('tmpdir/WBPaper00051486.json')
json_array = json.load(input_file)

# Create a variable that will take JSON and put it into a python dictionary
store_details = [
    ["name"],
    ["city"]
]

# Learn how to loop better =/
#for stores in [item["store_details"] for item in json_array]

for item in json_array:
    print (item["abstract"])

#for thing in cow["data"]:
#    print(thing["id"])

# Print my results
#print(store_details)


