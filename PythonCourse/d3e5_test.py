

import d3e5

f = open('random_words.1000.txt','r')

L = f.readlines()

f.close()

d3e5.selSort(L)


print L[0],L[1],L[999]

