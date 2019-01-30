#!/usr/bin/env python
from __future__ import print_function
import sys
import os.path
import argparse
import re
from subprocess import call

"""

Script for preparing and running deTIN

"""


epi = ('\
    \n\
	Make the preparation for deTIN, test tumour and normal sample\n\
     \n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes the commands for running deTIN tumour-in-normal subtraction', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-v', '--vcf', default=None, dest='vcf', action='store', required=True, help="Somatic VCF file")
parser.add_argument('-s', '--SVvcf', default=None, dest='SV', action='store', required=False, help="Somatic VCF SV file")
parser.add_argument('-g', '--genome', default=None, dest='gfa', action='store', required=True, help="genome.fa")
#parser.add_argument('-d', '--genome.dict', default=None, dest='gdic', action='store', required=True, help="genome.dict file")
parser.add_argument('-n', '--normal_bam', default=None, dest='norm', action='store', required=True, help="Normal.bam file full path")
parser.add_argument('-t', '--tumour_bam', default=None, dest='tum', action='store', required=True, help="Tumour.bam file full path")
#parser.add_argument('-p', '--pon', default=None, dest='pon', action='store', required=True, help="Panel of Normal CNV hdf5 file")
parser.add_argument('-o', '--out_folder', default=None, dest='out', action='store', required=True, help="Output folder full path")


gc_prof = '/home/mzarowiecki/bin/hartwigmedicalfoundation/GC_profile.hg38.1000bp.cnp'
jar_loc = '/home/mzarowiecki/bin/hartwigmedicalfoundation/hmftools_pipeline_v4_3'
circos_loc = '/home/mzarowiecki/bin/circos-0.69-6/bin/circos'
# Needs a BED file with human variants
human_bed = '/home/mzarowiecki/scratch/REF/af-only-gnomad.hg38.01.chr.bed'


# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()


# Check if input files exist
if not os.path.isfile(args.norm)==True:
    print("Cannot find input file ",args.norm)
    sys.exit(1)

# Check if input files exist
if not os.path.isfile(args.tum)==True:
    print("Cannot find input file ",args.tum)
    sys.exit(1)

# Create prefixes and files
norm=args.norm.rstrip()
normpx = norm.split("/")[-1]
normpx = ''.join(normpx.split(".bam")[0:-1])
tum=args.tum.rstrip()
tumpx = tum.split("/")[-1]
tumpx = ''.join(tumpx.split(".bam")[0:-1])


print ('mkdir %s' % (args.out))

# Tidy up VCF
print ('python /home/mzarowiecki/bin/Python_ditties/VCF_strelka_for_PURPLE.py -i %s' % (args.vcf) )
print ('cat %s.pysam.vcf  | tr -d \'$\' | bgzip > %s.pysam.vcf.gz ' % (args.vcf,args.vcf) )
#print ('python /home/mzarowiecki/bin/Python_ditties/VCF_strelka_for_PURPLE.py -i %s' % (args.SV) ')

# Run COBALT

print ('source activate r-3.5.0  # or use another method to get R-package biocLite("copy number") in your path')

# Fairly small and fast
print ('bsub -q bio -J COB%s  -o COB%s.o -e COB%s.e -n 1 -P Analysis java -jar %s/cobalt-1.5.jar -reference %s -reference_bam %s -tumor %s \
-tumor_bam %s -output_dir %s/cobalt -threads 24 -gc_profile %s' % (tumpx,tumpx,tumpx, jar_loc, normpx, args.norm, tumpx, args.tum, args.out,gc_prof ))


# Run AMBER


# Needs bsub
print ('module add Sambamba/0.6.6')
print ('bsub -P Analysis -q bio -J SB%s  -o SB%s.o -e SB%s.e -R "span[hosts=1]" -R "select[largedsk]" -R "hname!=\'hpc-prod-grid-lsfexec-001\' && hname!=\'hpc-prod-grid-lsfexec-002\'" -n 4 \
"sambamba mpileup -t 20 -L %s %s --samtools -q 1 -f %s \
> %s/%s.mpileup" ' % (normpx, normpx, normpx, human_bed, args.norm, args.gfa, args.out, normpx))

# Needs bsub
print ('bsub -q bio -P Analysis -J SB%s  -o SB%s.o -e SB%s.e -R "span[hosts=1]" -R "select[largedsk]" -R "hname!=\'hpc-prod-grid-lsfexec-001\' && hname!=\'hpc-prod-grid-lsfexec-002\'" -n 4 \
"sambamba mpileup -t 20 -L %s %s --samtools -q 1 -f %s \
> %s/%s.mpileup" ' % (tumpx, tumpx, tumpx, human_bed, args.tum, args.gfa, args.out,tumpx))


print ('bsub -q bio -P Analysis -J AMB%s  -o AMB%s.o -e AMB%s.e -n 1  java -jar %s/amber-1.7.jar -sample %s -output_dir %s/amber -reference %s/%s.mpileup \
-tumor %s/%s.mpileup' % (tumpx,tumpx,tumpx,jar_loc, tumpx,args.out,args.out,normpx,args.out,tumpx))



# Run PURPLE

if args.SV:
    print ('bsub -q bio -J PU%s  -o PU%s.o -e PU%s.e -n 1 java -jar %s/purple-2.16.jar  \
-run_dir %s -ref_sample %s -tumor_sample %s -threads 6 \
-gc_profile %s -amber %s/amber -cobalt %s/cobalt -baf %s \
-somatic_vcf %s.pysam.vcf.gz -structural_vcf %s -circos %s  \
-min_purity 0.1 -ref_genome hg38' % (tumpx,tumpx,tumpx,jar_loc,args.out,normpx,tumpx,gc_prof,args.out, args.out, args.out, args.vcf, args.SV,circos_loc))

else:
    print ('bsub -q bio -J PU%s  -o PU%s.o -e PU%s.e -n 1 java -jar %s/purple-2.16.jar  \
-run_dir %s -ref_sample %s -tumor_sample %s -threads 6 \
-gc_profile %s -amber %s/amber -cobalt %s/cobalt -baf %s \
-somatic_vcf %s.pysam.vcf.gz -circos %s  \
-min_purity 0.1 -ref_genome hg38' % (tumpx,tumpx,tumpx,jar_loc,args.out,normpx,tumpx,gc_prof,args.out, args.out, args.out, args.vcf,circos_loc))



quit()








