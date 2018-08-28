#! python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# Get the file name from the commandline.
if len(sys.argv) == 2:
    input_name = sys.argv[1]
    output_name =  input_name+".rc.fa" 


else:
    print("\n\n"
          "Usage: python3 revcomp.py\n\nResult in infile.rc.fa\n")
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

for seq_record in SeqIO.parse(input_name, "fasta"):
    my_head=seq_record.id
    my_seq = seq_record.seq

    #my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
    #my_seq
    #my_seq.complement()
    #out.write(">",my_head,"\n",my_seq.reverse_complement(), sep='')
    #my_seq = Seq(str(my_seq), IUPAC.unambiguous_dna)
    #print (my_seq)
    #print (str(my_seq))

    my_out = ">"+my_head+"\n"+my_seq.reverse_complement()    
    out.write(str(my_out) + "\n")


exit()



