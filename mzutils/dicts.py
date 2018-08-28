#!/usr/bin/env python

"""

This script reads in 


"""



def printDict(d):
    for k, v in d.items():
        if type(v) is dict:
            printDict(v)
        else:
            print("{0} : {1}".format(k, v))

