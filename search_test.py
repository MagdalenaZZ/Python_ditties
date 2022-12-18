import re
from wbpreader.readerdef import searchwords

sen = 'The RNA binding proteins PUF-5, PUF-6, and PUF-7 reveal multiple systems for maternal mRNA regulation during C. elegans oogenesis of one two three genes gene.'


# Ignore case
pattern = re.compile(r'RNA', flags=re.IGNORECASE)

# Any of one or two
pattern2 = re.compile(r'/\b(?:one)\b/gi')

# Gene only full word
gene='gene'
pattern3=re.compile(rf'\b{gene}\b')


# Gene and genes 
pattern4 = re.compile(r'(\b(gene|genes)\b.*){2,}')

# Gene and genes as part of word
pattern5 = re.compile(r'/protein/')

'''
if re.search(pattern, sen):
    print ( "RNA", sen)

if re.search(pattern2, sen):
    print ( "123", sen)

if re.search(pattern3, sen):
    print ( "gene word", sen)

if re.search(pattern4, sen):
    print ( "Multipe words" , sen)

# Just returns true/false
if 'prot' in sen:
    print ( "prot", sen)

if re.search('prot', sen):
    print ( "protS", sen)
'''

words =	{
    "exon":  re.compile(r'exon', flags=re.IGNORECASE),
    "intron": re.compile(r'intron', flags=re.IGNORECASE),
    "gene": re.compile(rf'\b{gene}\b'),
    "protein": re.compile(r'protein', flags=re.IGNORECASE),
    "prot": re.compile(r'prot', flags=re.IGNORECASE),
    "geneS": re.compile(r'(\b(gene|genes)\b.*){2,}')
}


res = searchwords(words, sen)

print (res)

