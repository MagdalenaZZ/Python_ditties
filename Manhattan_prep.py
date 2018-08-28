#! python3

import sys



# Get the xlsm file name from the commandline.
if len(sys.argv) == 2:
    input_name = sys.argv[1]
    bed = sys.argv[2]

    output_name =  input_name+".prot.fa" 
    output_name2 =  input_name+".untranslated.fa" 

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

OUT = open(output_name, "w")
OUT2 = open(output_name2, "w")

for seq_record in SeqIO.parse(input_name, "fasta"):
    my_head=seq_record.id
    try:
        my_seq = seq_record.seq.translate()
        my_out = ">"+my_head+"\n"+my_seq
        OUT.write(str(my_out) + "\n")
    except:
        print ("The following record have a character that is not ACGT", my_head)
        my_seq = seq_record.seq
        my_out = ">"+my_head+"\n"+my_seq
        OUT2.write(str(my_out) + "\n")

OUT.close()
OUT2.close()


