#!/usr/bin/env python


import sys
import argparse

# Describe what the script does
parser = argparse.ArgumentParser(description='Process some integers.')


# Add the numbers passed
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                   help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                   const=sum, default=max,
                   help='sum the integers (default: find the max)')

# Save the input integer in args
args = parser.parse_args()


print(args.integers)

# use function accumulate to add the input integers
print(args.accumulate(args.integers))





