#! python3

import sys
import gzip
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
import string

# Get the  file name from the commandline.
if len(sys.argv) == 3:
    patient_list = sys.argv[1]
    gene_list = sys.argv[2]

    output_name =  patient_list + "." + gene_list + ".generes"
    output_name2 =  patient_list + "." + gene_list + ".res"    


else:
    print("\n\n"
          "Usage: python3 patient_list gene_list\n\nResult in patient_list.gene_list.res\n")
    exit()


##########################


PatDict={}
PatSet=set()
fh = gzip.open('CosmicSample.tsv.gz', 'r')
for line in fh:
    line=str(line).strip()
    val=line.split('\\t')
    can=0

    try:
        i=int(val[28])
        can=1
        #print("try",i,val[3],val[0])
    except:
        #print("ex",val[28],val[3],val[0])
        pass
    #print(can)

    if can > 0:
        #print("if",val[28],val[3],val[0])
        if (i < int(patient_list)):
            print("if",int(patient_list),i,val[3],val[0])
            #PatDict[val[0]]=val[3]
            #PatSet.add(val[0])
        else:
                        pass
            #print("else",int(patient_list),i,val[3],val[0])
    else:
        pass
        #print("else",val[28],val[3],val[0])


"""
        if can > 0:
            if (i < patient_list):
                        print("if",i,val[3],val[0])
                else:
                        print("else",i,val[3],val[0])
"""
"""
    if val[28]==int:
        print("if",val[28],val[3],val[0])
    else:
        print("else",val[28],val[3],val[0])
"""
"""
    try:
        i = int(val[28])
        print("if",val[28],val[3],val[0])
        if (i < patient_list):
            print("if2",val[28],val[3],val[0])
        #    PatDict[val[0]]=val[3]
        #    PatSet.add(val[0])
        else:
            print("else",val[28],val[3],val[0])
    except:
        print("pass",val[28],val[3],val[0])
        pass

    #if (int(val[28]) < patient_list):
    #    print(val[28],val[3],val[0])
                #        PatDict[val[0]]=val[3]
                #        PatSet.add(val[0])
"""

"""

#try:
    #fh = open( patient_list, 'r' )
#    PatSet = set(line.strip() for line in open(patient_list,'r'))
#except IOError:
#   print ("Error: Can\'t find file or read data from patient file", patient_list)
"""

try:
    #fh2 = open( gene_list, 'r' )
    GeneSet = set(line.strip() for line in open(gene_list,'r'))
except IOError:
   print ("Error: Can\'t find file or read data from gene file", gene_list)


try:
    fh3 = gzip.open('CosmicCompleteCNA.tsv.gz', 'r')
    #with ZipFile('CosmicCompleteCNA.tsv.gz, 'r') as fh3:
    #    fh3.open('CosmicCompleteCNA.tsv.gz, 'r') as fh4

except IOError:
   print ("Error: Can\'t find file or read data from gene file CosmicCompleteCNA.tsv.gz")


try:
    out = open( output_name, 'w' )
except IOError as e:
    print ("Error: Can\'t write data to output file", output_name, "\n", e)
    sys.exit()
except IsADirectoryError as e:
    print ("Error: Can\'t write data to output file", output_name, "\n", e)
    sys.exit()

try:
    out2 = open( output_name2, 'w' )
except IOError as e:
    print ("Error: Can\'t write data to output file", output_name2, "\n", e)
    sys.exit()
except IsADirectoryError as e:
    print ("Error: Can\'t write data to output file", output_name2, "\n", e)
    sys.exit()


##########################

#print(GeneSet)

# For each CNA, keep line only if it is the right patient and gene

for line in fh3:
    line = str(line.strip())
    #print("no1" + line)
    fi = line.split('\\t')
    #print("no2" + fields[0])
    if fi[3] in PatSet:
        if fi[1] in GeneSet:
            print(fi[1],fi[3],PatDict[fi[3]],fi[10],fi[11],fi[12],sep='\t',file=out) 
        else:
            print(fi[3],PatDict[fi[3]],fi[10],fi[11],fi[12],sep='\t',file=out2 )



