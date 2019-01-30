#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
import pysam
import re


"""

Script for parsing a VCF file

"""


epi = ('\
    \n\
	VCF file parser, turning data into Mutect1-like format for deTIN. Have to give input that has both normal and tumour sample in it\n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and creates a Mutect1 like format for deTIN', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, help="VCF file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)


# Read input
output=args.vcf+".mutect1.stats"
out = open(output, "w")
#print("Input: ",args.vcf, "Output: " , output)


# read the input file
myvcf = pysam.VariantFile(args.vcf, "r")


print ("contig","position","ref_allele","alt_allele","tumor_name","normal_name","t_ref_count","t_alt_count",
"n_ref_count","n_alt_count","failure_reasons","judgement",sep="\t",file=out)

filen = re.split("[/]+", args.vcf)
names = re.split("[._]+", filen[-1])
tumour_name=names[2]+ "_" + names[3]
normal_name=names[0]+ "_" + names[1]
# print (normal_name, tumour_name)

for r in myvcf:

    #### FILTER OUT #####
    # Shared called total
    # Filter out sites which
    chr=r.chrom
    pos=r.pos
    id =r.id
    ref = r.ref
    alt = r.alts
    qual=r.qual
    filter=r.filter.keys()
    info=r.info.keys()
    failure_reasons=str(filter[0])
    judgement="KEEP"


    # Calculate failure and judgement
    # If it is PASS, make it keep
    # If it is failed, make it a "fail" which allows for reassignment
    if re.match("PASS",failure_reasons):
        failure_reasons=''
        judgement = "KEEP"
    # Check if there is LowQScore fail
    elif re.match("LowQscore",failure_reasons):
        qss=0
        vqsr=0
        # If it is, check if INFO['QSS'] > 20 and if INFO['VSQR'] > 2.5
        #print (info['QSS'], info['VSQR'])
        if "QSS" in r.info.keys():
            qss=r.info["QSS"]
        else:
            qs=0

        if "VQSR" in r.info.keys():
            vqsr=r.info["VQSR"]
        else:
            vqsr=0

        if (vqsr>2.5) and (qss> 20) :
            failure_reasons="alt_allele_in_normal"
            judgement = "REJECT"
        else:
            failure_reasons="str_contraction"
            judgement = "REJECT"

    else:
        failure_reasons="t_lod_fstar"
        judgement = "REJECT"


    # Make more sane alt
    if alt:
        alt=str(alt[0])
    else:
        alt="."

    altb = alt

    trefd = 0
    taltd = 0
    nrefd = 0
    naltd = 0
    # Switch between D1 and D2 reads counts
    d0=1

    # If is SNP, which Mutect1 outputs, calculate read coverage
    if 'CU' in r.samples[0].keys():

        if (r.ref == 'A'):
            nrefd = r.samples[0]['AU'][d0]
            trefd = r.samples[1]['AU'][d0]
        elif (r.ref == 'C'):
            nrefd = r.samples[0]['CU'][d0]
            trefd = r.samples[1]['CU'][d0]
        elif (r.ref == 'G'):
            nrefd = r.samples[0]['GU'][d0]
            trefd = r.samples[1]['GU'][d0]
        elif (r.ref == 'T'):
            nrefd = r.samples[0]['TU'][d0]
            trefd = r.samples[1]['TU'][d0]
        else:
            print("WARN: ", r.ref, " is not A,C,G or T :", r.id)

        if (altb == 'A'):
            naltd = r.samples[0]['AU'][d0]
            taltd = r.samples[1]['AU'][d0]
        elif (altb == 'C'):
            naltd = r.samples[0]['CU'][d0]
            taltd = r.samples[1]['CU'][d0]
        elif (altb == 'G'):
            naltd = r.samples[0]['GU'][d0]
            taltd = r.samples[1]['GU'][d0]
        elif (altb == 'T'):
            naltd = r.samples[0]['TU'][d0]
            taltd = r.samples[1]['TU'][d0]
        # If ALT is unknown/N, add all tier1 counts together and subtract REF counts
        elif (altb == '.'):
            naltd = r.samples[0]['AU'][d0]+ r.samples[0]['CU'][d0]+ r.samples[0]['GU'][d0]+ r.samples[0]['TU'][d0]-nrefd
            taltd = r.samples[1]['AU'][d0]+ r.samples[1]['CU'][d0]+ r.samples[1]['GU'][d0]+ r.samples[1]['TU'][d0]-trefd
        else:
            print("WARN: ", r.ref, " is not A,C,G or T :", r.id)

    t_ref_count=trefd
    t_alt_count=taltd
    n_ref_count=nrefd
    n_alt_count=naltd




    #r.info["VT"] = "None"
    print (chr,pos,ref,alt,tumour_name,normal_name,t_ref_count,t_alt_count,n_ref_count,n_alt_count,failure_reasons,judgement,sep="\t",file = out)


out.close()
quit()



# Complete mutect1 output
#print ("contig","position","context","ref_allele","alt_allele","tumor_name","normal_name","score","dbsnp_site","covered","power",
#       "tumor_power","normal_power","normal_power_nsp","normal_power_wsp","total_reads","map_Q0_reads","init_t_lod","t_lod_fstar",
#       "t_lod_fstar_forward","t_lod_fstar_reverse","tumor_f","contaminant_fraction","contaminant_lod","t_q20_count","t_ref_count",
#       "t_alt_count","t_ref_sum","t_alt_sum","t_ref_max_mapq","t_alt_max_mapq","t_ins_count","t_del_count","normal_best_gt",
#       "init_n_lod","normal_f","n_q20_count","n_ref_count","n_alt_count","n_ref_sum","n_alt_sum",
#       "power_to_detect_positive_strand_artifact","power_to_detect_negative_strand_artifact","strand_bias_counts",
#       "tumor_alt_fpir_median","tumor_alt_fpir_mad","tumor_alt_rpir_median","tumor_alt_rpir_mad","observed_in_normals_count",
#       "failure_reasons","judgement","total_counts_coverage_less_8","total_counts_coverage_greater_8",
#       "alt_count_greater1_af_greater_01percent","alt_count_greater2_af_greater_03percent","alt_count_greater3_af_greater_1percent",
#       "alt_count_greater3_af_greater_3percent","alt_count_greater3_af_greater_20percent","alt_count_greater10_af_greater_20percent",
#       "PoN_Germline","PoN_Artifact","bad",sep="\t")



