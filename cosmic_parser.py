#! python3

import sys
import gzip
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
import string

# Get the  file name from the commandline.
if len(sys.argv) == 3:
    patient_list = sys.argv[2]
    gene_list = sys.argv[1]

    output_name =  gene_list + "." + patient_list + ".generes"
    output_name2 =  gene_list  + "." + patient_list + ".res"    
    output_name3 =  gene_list  + "." + patient_list + ".pats"

else:
    print("\n\n"
          "Usage: python3 gene_list patient_age\n\nResult in gene_list.patient_age.res\n")
    exit()


##########################

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

try:
    out3 = open( output_name3, 'w' )
except IOError as e:
    print ("Error: Can\'t write data to output file", output_name3, "\n", e)
    sys.exit()
except IsADirectoryError as e:
    print ("Error: Can\'t write data to output file", output_name3, "\n", e)
    sys.exit()

######################

#print("Start parse")
PatDict={}
PatSet=set()
fh = gzip.open('CosmicSample.tsv.gz', 'r')
print("Age","SAMPLE_ID","Patient_ID",file=out3)

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
			val[0]=val[0].replace('b\'','')
			val[3]=val[3].replace('b\'','')
			print(i,val[3],val[0],file=out3)
			PatDict[val[0]]=val[3]
			PatSet.add(val[0])
	
		else:
                        pass
			#print("else",int(patient_list),i,val[3],val[0])
	else:
		pass
		#print("else",val[28],val[3],val[0])

try:
    #fh2 = open( gene_list, 'r' )
    GeneSet = set(line.strip() for line in open(gene_list,'r'))
except IOError:
   print ("Error: Can\'t find file or read data from gene file", gene_list)


try:
	fh3 = gzip.open('CosmicCompleteCNA.tsv.gz', 'r')
	#with ZipFile('CosmicCompleteCNA.tsv.gz, 'r') as fh3:
	#	fh3.open('CosmicCompleteCNA.tsv.gz, 'r') as fh4

except IOError:
   print ("Error: Can\'t find file or read data from gene file CosmicCompleteCNA.tsv.gz")

#print("End parse")

##########################

#print(GeneSet)

# For each CNA, keep line only if it is the right patient and gene
#for ent in GeneSet:
	#ent=ent.replace('b\'','')
	#print(type(ent),ent)
	#pass

print("ID_GENE","ID_SAMPLE","ID_Patient","Hist","Hist1","Hist2","Hist3","TOTAL_CN","MUT_TYPE",sep='\t' ,file=out)
print("ID_SAMPLE","ID_Patient","Hist","Hist1","Hist2","Hist3",sep='\t' ,file=out2)

for line in fh3:
	line = str(line.strip())
	#print("no1" + line)
	fi = line.split('\\t')
	#print("no2" + fields[0])
	if fi[3] in PatDict.keys():
		#print("Match",fi[3])
		if fi[1] in GeneSet:
			print(fi[1],fi[3],PatDict[fi[3]],fi[9],fi[10],fi[11],fi[12],fi[14],fi[16],sep='\t',file=out)
		else:
			print(fi[3],PatDict[fi[3]],fi[9],fi[10],fi[11],fi[12],sep='\t',file=out2 )
	else:
		pass
		#print("No match",type(fi[3]),fi[3])

exit()







