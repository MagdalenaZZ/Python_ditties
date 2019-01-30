#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import os.path
import argparse
import csv
import pysam
import re
from subprocess import call
from pathlib import Path


"""

Script for taking strelka2 output and creating ISO-validation of small variants

"""


epi = ('\
    \n\
	File parser,  VCF files\n\
    \n\
')

# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a Strelka2 VCF file and creates recall and precision stats', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, nargs = '+' , help="VCF.gz file(s)")
#parser.add_argument('-b', '--bed', default=None, dest='bed', action='store', required=False, help="BED boundary file")
#parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, help="VCF.gz file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)


args  = parser.parse_args()

prefix= args.vcf[0].split('/')[-1]
px=prefix.split('.')[0]



# First take the files, and do filtering, formatting
python ~/bin/Python_ditties/VCF_filter_format.py -i LP3000396-DNA_D01.somatic.indels.vcf.gz LP3000396-DNA_D01.somatic.snvs.vcf.gz

# Run the resulting script
bash LP3000396-DNA_B01.doTx.sh

# Merge with Tx data
bcftools merge  --force-samples  -O z -o LP3000396-DNA_A01.merge.vcf.gz /genomes/analysis/cancer/new_pipeline/ISO_validation/TracerXtruthSet/LP3000396-DNA_A01.tx.vcf.gz LP3000396-DNA_A01.compTx.vcf.gz


source activate py27.12
# Run recall stats and move files
ls *merge.vcf.gz | awk '{print "python ~/git/TRACERx_validation/VCF_summary_stats.py -i "$1" > "$1".bed" }' > recall.sh
bash recall.sh
mkdir Recall
mv *merge.vcf.gz.* Recall/
cd Recall
paste *merge.vcf.gz.sstats  | cut -f1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48 > merge.stats
cd ..

# Run precision stats and move files
ls *merge.vcf.gz  | awk '{print "python ~/git/TRACERx_validation/VCF_summary_precision.py -i "$1" > "$1".bed" }' > prec.sh
bash prec.sh
mkdir Precision
mv *merge.vcf.gz.* Precision/
cd Precision
paste *merge.vcf.gz.sstats  | cut -f1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48 > merge.stats
cd ..



bed='/genomes/scratch/magz/REF/SureSelectV5_TRACERx_Edition.padded.reduced.hg38.bed'


# create an object of new bed file and open in to write data.
output=px+".doTx.sh"
out = open(output, 'w')

#print (prefix,px)

vcfs= ' '.join(args.vcf)

# Collate input files, if there are several
print ("bcftools concat -a -D -O z -o %s.somatic.vcf.gz  %s" % (px, vcfs),file=out)
print ("bcftools index  %s.somatic.vcf.gz" % (px),file=out)

# Filter by BED, if BED exists
#bcftools view -R /genomes/scratch/magz/REF/SureSelectV5_TRACERx_Edition.padded.reduced.hg38.bed -O z -o ../3_VCFs_2/LN3999999-DNA_A01_LP3000396-DNA_A01.TxPASS.vcf.gz ../3_VCFs_2/LN3999999-DNA_A01_LP3000396-DNA_A01.somatic.vcf.gz
print ("bcftools view -R %s -O z -o %s.filter.vcf.gz %s.somatic.vcf.gz" % (bed,px,px),file=out)

# Keep only the tumour
# bcftools view -s LP3000396-DNA_E01 -O z -o LN3999999-DNA_E01_LP3000396-DNA_E01.Txfil.vcf.gz LN3999999-DNA_E01_LP3000396-DNA_E01.TxPASS.vcf.gz
print ("bcftools view -s TUMOR -O z -o %s.filt.to.vcf.gz %s.filter.vcf.gz" % (px,px ),file=out)

# Reformat to general VCF format - to be compatile with Tx
# python ~/git/TRACERx_validation/VCF_our_reform_to_general.py -i LN3999999-DNA_A01_LP3000396-DNA_A01.Txfil.vcf.gz
print ("python ~/git/TRACERx_validation/VCF_our_reform_to_general.py -i %s.filt.to.vcf.gz" % (px),file=out)

# Remove redundant fields
# bcftools  annotate -x "INFO/ALTMAP,INFO/ALTPOS,INFO/DP,INFO/IC,INFO/IHP,INFO/MQ,INFO/MQ0,INFO/NT,INFO/OVERLAP,INFO/PNOISE,INFO/PNOISE2,INFO/QSI,INFO/QSI_NT,INFO/QSS_NT,INFO/RC,INFO/ReadPosRankSum,INFO/RU,INFO/SGT,INFO/SNVSB,INFO/TQSI,INFO/TQSI_NT,INFO/TQSS,INFO/TQSS_NT,INFO/AF1000G,INFO/AA,INFO/GMAF,INFO/cosmic,INFO/clinvar,INFO/EVS,INFO/RefMinor,INFO/phyloP,INFO/CSQT,INFO/CSQR,INFO/OLD_MULTIALLELIC,INFO/OLD_VARIANT,FORMAT/AU,FORMAT/CU,FORMAT/GU,FORMAT/TU,FORMAT/FDP,FORMAT/SDP,FORMAT/SUBDP,FORMAT/DP2,FORMAT/TAR,FORMAT/TIR,FORMAT/TOR,FORMAT/DP50,FORMAT/FDP50,FORMAT/SUBDP50" LN3999999-DNA_A01_LP3000396-DNA_A01.Txfil.vcf.gz.pysam.vcf.gz | bgzip > LN3999999-DNA_A01_LP3000396-DNA_A01.TIN.vcf.gz
print ("bcftools  annotate -x \"INFO/DP,INFO/IC,INFO/IHP,INFO/MQ,INFO/MQ0,INFO/NT,INFO/OVERLAP,INFO/PNOISE,INFO/PNOISE2,INFO/QSI,INFO/QSI_NT,INFO/QSS_NT,INFO/RC,INFO/ReadPosRankSum,INFO/RU,INFO/SGT,INFO/SNVSB,INFO/TQSI,INFO/TQSI_NT,INFO/TQSS,INFO/TQSS_NT,FORMAT/AU,FORMAT/CU,FORMAT/GU,FORMAT/TU,FORMAT/FDP,FORMAT/SDP,FORMAT/SUBDP,FORMAT/DP2,FORMAT/TAR,FORMAT/TIR,FORMAT/TOR,FORMAT/DP50,FORMAT/FDP50,FORMAT/SUBDP50\" %s.filt.to.vcf.gz.pysam.vcf.gz | bgzip > %s.pruned.vcf.gz"  % (px, px),file=out)
print ("vt sort -o %s.compTx.vcf.gz %s.pruned.vcf.gz" % (px,px),file=out)
print ("bcftools index  %s.compTx.vcf.gz" % (px),file=out)
print ("rm %s.filter.vcf.gz %s.filt.to.vcf.gz %s.filt.to.vcf.gz.pysam.vcf.gz %s.pruned.vcf.gz" % (px,px,px,px),file=out)


out.close()


print ("bash %s" % (output))


"""



# Merge the file with ALT variants
alts='/home/mzarowiecki/scratch/REF/allASDPs.SNV.50_10.valid.vcf.gz'

call(["bcftools", "merge","--force-samples","-O","z","-o",args.vcf+".asdp.vcf.gz",args.vcf,alts])









# read the input file
myvcf = pysam.VariantFile(args.vcf+".asdp.vcf.gz", "r")

# create an object of new bed file and open in to write data.
output=args.vcf+".asdp.res.vcf.gz"
out = open(output +'.review', 'w')
vaf = open(output +'.vaf', 'w')



myvcf.header.info.add("ALT", "1", "String", "Is variant on ALT or primary")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(output, 'w', header=myvcf.header)



# First parse through VCF file and pick out all SNVs in ALTs.

res={}
# res[alt][hom/het/amb/none]=number


for r in myvcf:

    #### FILTER OUT #####
    # Shared called total
    # Filter out sites which
    chr = r.chrom
    pos = r.pos
    id = str(r.id)
    #varID=':'.join([id.split(":")[0],id.split(":")[1]])
    #altb = r.ref
    #altb = r.alts
    score = r.qual
    filter = r.filter
    info = r.info
    format = r.format
    samples = r.samples
    end = r.stop # r.info["END"]
    strand='.'
    r.info['ALT']="NAN"
    
    # Test if asdp exists
    #if 'AL' in r.info.keys():
        #print (r.info['AL'])
    # Test if we have that variant
#    if 'DP' in r.format.keys():
#        print (r.format.keys())


    # Test if both info exists
    if 'AL' in r.info.keys() and 'AD' in r.format.keys():


        # Filter complex variants
        if len( r.samples[0]['AD']) > 2:
            continue
            #print(r.samples[0]['AD'])

        # Classify variant
        fq = r.samples[0]['AD'][1]/(r.samples[0]['AD'][1]+r.samples[0]['AD'][0])

        # If it is certainly a REF
        if fq<0.1:
            r.info['ALT']="REF"
        # If it is certainly a HOM ALT
        elif fq>=0.9:
            r.info['ALT']="HOM"
        # If it is certainly a HET PA/ALT
        elif fq>0.30 and fq<0.7:
            r.info['ALT']="HET"
        else:
            r.info['ALT']="AMB"

        print(r.info['AL'], r.info['ALT'],fq,r.samples[0]['AD'][0],r.samples[0]['AD'][1], file=vaf, sep='\t')
        vcf_out.write(r)
        

    elif 'AL' in r.info.keys():
        r.info['ALT']="REF"

    else:
        continue

# Populate summary
    if r.info['AL'] in res:
            
        if r.info['ALT'] in res[r.info['AL']]:
            res[r.info['AL']][r.info['ALT']] += 1
        else:
            res[r.info['AL']][r.info['ALT']] = 1

    else:
        res[r.info['AL']]={}
        #res[r.info['AL']]['NON']=0
        res[r.info['AL']]['HOM']=0
        res[r.info['AL']]['HET']=0
        res[r.info['AL']]['AMB']=0
        res[r.info['AL']]['REF']=0
        res[r.info['AL']][r.info['ALT']] = 1



# Then summarise the score

print ("Alt","Max","AMB","HET","HOM","REF","Sum","Verdict",sep='\t',end='',file=out)

for key in sorted(res):
    #print(key,end='\t',file=out) 
    max2 = max(res[key], key=res[key].get)
    print('\n',key,'\t',max2,sep='',end='',file=out)
    sum2=0
    for key2 in sorted(res[key]):
                print('\t',res[key][key2],sep='',end='',file=out) 
                sum2+=res[key][key2]

    print('\t',sum2,sep='',end='',file=out)


    # Try to classify
    fq=res[key][max2]/sum2

    # The max value is larger than all the others together
    if fq>0.8:
        print('\t',max2,sep='',end='',file=out)
    elif fq>0.6 and max2=='HET':
        print('\t',max2,sep='',end='',file=out)
    else:
        print('\t','AMBI',sep='',end='',file=out)


    

print ('\n',sep='\t',end='',file=out)







 
vaf.close()

vcf_out.close()


"""

exit(0)
