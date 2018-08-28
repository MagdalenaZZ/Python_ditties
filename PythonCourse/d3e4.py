


class wordFile:

    """An object aassociated with a text file"""
    # The initiation takes a file and turns it into a dictionary with first 4 letters of word as key, and row in file as value
    def __init__ (self, fileName = 'random_words.txt'):
        file  = open(fileName, 'r')
        self.wordList = file.readlines()
        file.close()
        self.refDict={}
        chars = 'abcdefghij'
        for c1 in chars:
            for c2 in chars:
                for c3 in chars:
                    for c4 in chars:
                        self.refDict[c1+c2+c3+c4] = []

        for i in xrange(len(self.wordList)):
            self.refDict[self.wordList[i][0:4]].append(i)

    def wordSearch(self, word):
        
        """Word search module"""
        print "Initiating wordSearch"
        n = len(word)

        if n<4:
            print "Too short word", word
            return []

        matches=[]

        for i in self.refDict[word[:4]]:
            #print i
            if word == self.wordList[i][:n]:
                matches.append(self.wordList[i][:-1])
        return matches


