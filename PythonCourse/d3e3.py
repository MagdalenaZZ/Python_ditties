
"""Word search module"""


fileName = 'random_words.txt'
file  = open(fileName, 'r')
wordList = file.readlines()
file.close()
refDict={}
chars = 'abcdefghij'
for c1 in chars:
        for c2 in chars:
            for c3 in chars:
                for c4 in chars:
                    refDict[c1+c2+c3+c4] = []

for i in xrange(len(wordList)):
    refDict[wordList[i][0:4]].append(i)



def wordSearch(word):

    """Word search module"""
    print "Initiating wordSearch"
    matches=[]
    n = len(word)

    if n<4:
        print "if"
        return matches

    for i in refDict[word[:4]]:
        print i
        if word == wordList[i][:n]:
            matches.append(wordList[i][:-1])
    return matches



