#! python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

# Get the xlsm file name from the commandline.
if len(sys.argv) == 2:
    input_name = sys.argv[1]
    output_name =  input_name+".s.fa" 


else:
    print("\n\n"
          "Usage: python3 fasta2singleLine.py\n\nResult in infile.s.fa\n")
    exit()


##########################

try:
    fh = open( input_name, 'r' )
except IOError:
   print ("Error: Can\'t find file or read data from input file", input_name)

try:
    out = open( output_name, 'w' )
except IOError as e:
    print ("Error: Can\'t write data to output file", output_name, "\n", e)
    sys.exit()
except IsADirectoryError as e:
    print ("Error: Can\'t write data to output file", output_name, "\n", e)
    sys.exit()

##########################


#print("Input:",input_name)

OUT = open(output_name, "w")

for seq_record in SeqIO.parse(input_name, "fasta"):

    my_head=seq_record.id
    my_seq = str(seq_record.seq)
    my_out = ">"+my_head+"\n"+my_seq  
    out.write(str(my_out) + "\n")

exit()



