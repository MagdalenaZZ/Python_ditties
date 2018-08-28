# print random data for a number of years

import random

# do a random seed
random.seed(123)

# do a starting year
year=1970

# print header
print "year number float"

# print random numbers
while year < 2015:
    n = random.randint(1,100000)
# the dot afterwards gives the digits after the dot
    x = random.randint(1,100000)/1000.
    print "%d %5d %6.3f" % (year, n ,x)
    year = year+1

print "Done"


