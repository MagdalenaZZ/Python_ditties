#!/usr/bin/python3


import sys
import re

start = sys.argv[1]
length= sys.argv[2]
length= int(length)-1

# Increment char (a -> b, az -> ba)
def inc_char(text, chlist = 'abcdefghijklmnopqrstuvwxyz'):
    # Unique and sort
    chlist = ''.join(sorted(set(str(chlist))))
    chlen = len(chlist)
    if not chlen:
        return ''
    text = str(text)
    # Replace all chars but chlist
    text = re.sub('[^' + chlist + ']', '', text)
    if not len(text):
        return chlist[0]
    # Increment
    inc = ''
    over = False
    for i in range(1, len(text)+1):
        lchar = text[-i]
        pos = chlist.find(lchar) + 1
        if pos < chlen:
            inc = chlist[pos] + inc
            over = False
            break
        else:
            inc = chlist[0] + inc
            over = True
    if over:
        inc += chlist[0]
    result = text[0:-len(inc)] + inc
    return result

##################


print(start)

for i in range(length):
    print(inc_char(start))
    start = inc_char(start)





