
def wordSearch(word, fileName = 'random_words.txt'):
    file  = open(fileName, 'r')
    wordList = file.readlines()
    #print type(wordList)
    file.close()
    n= len(word)
    row =1
    for x in wordList:
        if len(x) > n and x[:n] == word:
            return "Match found: %s at line %d" % (x[:-1],row)
        row +=1

res = wordSearch('beegea', 'random_words.10.txt')

print res

