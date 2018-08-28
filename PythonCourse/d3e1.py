import random


def wordGenerator(maxLength=12):
    """This function returns a random set of letters length as you define, but between 4 and 12"""
    s=''
    wordLength=random.randint(4,maxLength)
    for i in range(wordLength):
        # return random integer
        s += chr(random.randint(ord('a'), ord('j')))
    s += "\n"
    return s

#print (wordGenerator(14))



def writeWords(fileName = 'random_words.txt', wordNumber=10):
    file = open(fileName, 'w')
    for i in range(wordNumber):
        file.write(wordGenerator())
    file.close()




