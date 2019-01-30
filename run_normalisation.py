#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import argparse
import shutil

"""

Script for doing Bertha normalisation outside Bertha pipeline


"""


epi = ('\
    \n\
	Allowing for Bertha advanced filtering and normalisation of VCF files\n\nDo NOT worry about the error messages, the files will be produced anyway\n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and applies Bertha filtering\nFirst do: source /tools/bertha-test/bertha-compute/1.15.0/bin/activate\n\n', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default='file.vcf.gz', dest='vcf', action='store', required=False, help="full path to VCF file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)


# Add paths if needed
sys.path.append('/home/mzarowiecki/git/bertha/code/bertha-commons')
os.environ["TIME_COMMAND"] = "/usr/bin/time"


# Create default params dict

print ("Making params...")
params = {'input_files': {'delivery_type': 'Cancer_Normal', 'files': {'input_vcf':'file.vcf.gz' , 'vcf_file': 'sort.vcf.gz' , 'non_ambiguous_ref_genome_fasta' : '/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa' , 'md5_manifest' : 'md5.manifest'} ,   'delivery_format_version' : 'DeliveryFormatVersion.V1_5' },  'output_files': {'files': {'filter_output': 'filtered.vcf.gz', 'decompose_output_uncompressed': 'decomposed.uncompressed.vcf', 'vt_log': 'vt.log', 'bgzip_warning': 'bgzip.log', 'decompose_output': 'decomposed.uncompressed.vcf.gz' , 'bcftools_warning': 'bcftools.log' , 'left_align_indels_output_uncompressed' : 'left.vcf', 'left_align_indels_output' : 'left.vcf.gz' , 'remove_duplicates_output_uncompressed' : "rmdup.vcf", 'remove_duplicates_output': 'rmdup.vcf.gz', 'filter_output_uncompressed' : 'fil2.vcf' , 'filter_output' :  'fil2.vcf.gz' , 'sort_vcf_output' : 'sort.vcf.gz' , 'vt_sort_log' : 'vt.sort.log' , 'decompose_output_index' : 'decomposed.uncompressed.vcf.gz.csi' , 'left_align_indels_output_index' : 'left.vcf.gz.csi' , 'remove_duplicates_output_index' : 'rmdup.vcf.gz.csi' , 'sort_vcf_output_index' : 'sort.vcf.gz.csi' , 'vcf_file': 'final.vcf.gz' , 'vcf_index' : 'final.vcf.gz.tbi' , 'bcftools_log' : 'bcftools.log2' , 'tabix_log': "logfol/tabix.log" }, 'dirs': {'base_dir': '/home/mzarowiecki' } },  'config': {'bcftools': 'bcftools', 'vt': '/tools/bertha/vt-0.ee9a751/bin/vt', 'bgzip': 'bgzip', 'normalisation_window_size' : 1000 , 'tabix' : 'tabix'  } , 'sample_id' : ['LP3000000-DNA_A01','LP3000000-DNA_A01'] , 'variant_caller' : 'strelka' , 'selection_method' : 'af', 'sample_name' : 'LP3000000-DNA_A01'  }


# Get the name of the input samples and replace details in params

from executionutils.helpers import (
    timed_command,
    check_file_exists_and_nonempty,
)

# Get the sample names from the file
profile, samples = timed_command( '{bcftools} view -h {vcf_file} | tail -1 | cut -f10- '.format(
    bcftools=params['config']['bcftools'],
    vcf_file=args.vcf,
    ),
    params['output_files']['dirs']['base_dir'],
    params['output_files']['files']['bcftools_log'],
)
samples = samples.strip()
params['sample_name'] = samples.split('\t')[1]
params['sample_id'][0] = samples.split('\t')[0]
params['sample_id'][1] = samples.split('\t')[1]


here  = os.getcwd()
params['output_files']['dirs']['base_dir']=here
#params['output_files']['dirs']['base_dir']='/'.join((here, params['sample_id'][1], 'normalise_vcf', '12345', '1' ))
params['output_files']['dirs']['log_dir']='/'.join((here, params['sample_id'][1], 'normalise_vcf', '12345', '1', 'logs' ))
params['output_files']['dirs']['output_dir']='/'.join((here, params['sample_id'][1], 'normalise_vcf', '12345', '1', 'output' ))

logfol =params['output_files']['dirs']['log_dir'] 
od =params['output_files']['dirs']['output_dir'] 

if not os.path.isdir(params['output_files']['dirs']['base_dir']):
    os.makedirs(params['output_files']['dirs']['base_dir'])


if not os.path.isdir(params['output_files']['dirs']['log_dir']):
    os.makedirs(params['output_files']['dirs']['log_dir'])


if not os.path.isdir(params['output_files']['dirs']['output_dir']):
    os.makedirs(params['output_files']['dirs']['output_dir'])



#shutil.copy(args.vcf, od)
#shutil.copy(''.join((args.vcf,'.tbi')), od)
#os.chdir(od)

#params['input_files']['delivery_type']='Cancer_Normal'
params['input_files']['files']['input_vcf'] = args.vcf 
params['input_files']['files']['vcf_file'] = ''.join((params['sample_id'][1],'.lexsort.duprem.left.split.vcf.gz')) #
params['input_files']['files']['non_ambiguous_ref_genome_fasta'] = '/genomes/bertha-prod/resources/bertha/data/GRCh38Decoy/reference/GRCh38Decoy_no_alt_no_ambiguous.fa' 
params['input_files']['files']['md5_manifest'] = 'md5.manifest' 
params['input_files']['files']['delivery_format_version']='DeliveryFormatVersion.V1_5' 

#params['output_files']['files']['filter_output'] = ''.join((params['sample_id'][1],'.filtered.vcf.gz')) #

# For filter_variant_records_with_missing_alternative_allele_values(params)
params['output_files']['files']['filter_output_uncompressed']=''.join((params['sample_id'][1],'.filtered.vcf')) #
params['output_files']['files']['filter_output']=''.join((params['sample_id'][1],'.filtered.vcf.gz')) 
# "/genomes/bertha-prod/analysis/singlesample/LP3000396-DNA_A01/CANCP41877/normalise_vcf/116332/1/output/LP3000396-DNA_A01.filtered.vcf.gz"

# For decomp_multi_allelic_variants(params)
params['output_files']['files']['decompose_output_uncompressed']=''.join((params['sample_id'][1],".split.vcf"))
params['output_files']['files']['decompose_output']=''.join((params['sample_id'][1],".split.vcf.gz"))
params['output_files']['files']['decompose_output_index']=''.join((params['sample_id'][1],".split.vcf.gz.csi")) #

params['output_files']['files']['vt_log']='/'.join((logfol, 'vt.log')) #
params['output_files']['files']['bgzip_warning']='/'.join((logfol,  'bgzip.warning')) #
params['output_files']['files']['bcftools_warning']='/'.join((logfol,  'bcftools.warning')) #
# For left_align_indels(params)
params['output_files']['files']['left_align_indels_output_uncompressed']= ''.join((params['sample_id'][1], '.left.split.vcf')) #
params['output_files']['files']['left_align_indels_output' ]=''.join((params['sample_id'][1], '.left.split.vcf.gz')) #
params['output_files']['files']['left_align_indels_output_index']=''.join((params['sample_id'][1], '.left.split.vcf.gz.csi')) #
# "left_align_indels_output_index": "/genomes/bertha-prod/analysis/singlesample/LP3000396-DNA_A01/CANCP41877/normalise_vcf/116332/1/output/LP3000396-DNA_A01.left.split.vcf.gz.csi"

# For remove_duplicates(params)
params['output_files']['files']['remove_duplicates_output_uncompressed']=''.join((params['sample_id'][1], '.lexsort.duprem.left.split.vcf')) #
params['output_files']['files']['remove_duplicates_output']=''.join((params['sample_id'][1], '.lexsort.duprem.left.split.vcf.gz')) #
params['output_files']['files']['remove_duplicates_output_index']=''.join((params['sample_id'][1], '.lexsort.duprem.left.split.vcf.gz.csi')) #
# "remove_duplicates_output_index": "/genomes/bertha-prod/analysis/singlesample/LP3000396-DNA_A01/CANCP41877/normalise_vcf/116332/1/output/LP3000396-DNA_A01.lexsort.duprem.left.split.vcf.gz.csi"

# For sort(params)
params['output_files']['files']['sort_vcf_output']=''.join((params['sample_id'][1],'.duprem.left.split.vcf.gz')) #
params['output_files']['files']['vt_sort_log']='/'.join((logfol, 'vt_sort.log')) #
params['output_files']['files']['sort_vcf_output_index']=''.join((params['sample_id'][1],'.duprem.left.split.vcf.gz.csi')) #
# "sort_vcf_output_index": "/genomes/bertha-prod/analysis/singlesample/LP3000396-DNA_A01/CANCP41877/normalise_vcf/116332/1/output/LP3000396-DNA_A01.duprem.left.split.vcf.gz.csi"

# Final
params['output_files']['files']['vcf_file']=''.join((params['sample_id'][1], '.duprem.left.split.reheadered.vcf.gz')) # final
params['output_files']['files']['vcf_index']=''.join((params['sample_id'][1], '.duprem.left.split.reheadered.vcf.gz.tbi')) # final
params['output_files']['files']['bcftools_log']='/'.join((logfol,'bcftools.log')) ##
params['output_files']['files']['tabix_log']='/'.join((logfol,"tabix.log")) ##

#params[ 'config': {'bcftools': 'bcftools', 'vt': '/tools/apps/vt/ee9a751/bin/vt', 'bgzip': 'bgzip', 'normalisation_window_size' : 1000 , 'tabix' : 'tabix'  } ,
#params[    'sample_id' : ['LP3000000-DNA_A01','LP3000000-DNA_A01'] , 
#params[    'variant_caller' : 'strelka' , 
#params[    'selection_method' : 'af', 
#params[    'sample_name' : 'LP3000000-DNA_A01'  }


# Import all the Bertha modules for normalise_VCF

print ("Normalise VCF...")
sys.path = ['/home/mzarowiecki/git/bertha/code/bertha-compute/analysis/pipeline/components/normalise_vcf'] + sys.path
print (sys.path)

#import component
from component import main
from component import filter_variant_records_with_missing_alternative_allele_values
from component import fix_vcf_header
from component import decomp_multi_allelic_variants
#from component import decompose_output_uncompressed
from component import left_align_indels
from component import remove_duplicates
from component import sort
from component import check_all_files_exist

# Run normalise_VCF
main(params)

# Delete the function main, so it can be replaced
del locals()['main']
del locals()['filter_variant_records_with_missing_alternative_allele_values']
del locals()['fix_vcf_header']
del locals()['decomp_multi_allelic_variants']
del locals()['left_align_indels']
del locals()['remove_duplicates']
del locals()['sort']
del locals()['check_all_files_exist']
#del component



# Move outputfiles to output
files = os.listdir(here)

for f in files:
    if (f.startswith( ''.join((params['sample_id'][1] , '.' )) )):
        shutil.move(f, od)






#del globals()['main']

# Import all the bertha modules for vcf_cleanup_somatic
print ("VCF cleanup...")
sys.path.remove('/home/mzarowiecki/git/bertha/code/bertha-compute/analysis/pipeline/components/normalise_vcf')
sys.path = ['/home/mzarowiecki/git/bertha/code/bertha-compute/analysis/pipeline/components/vcf_cleanup_somatic'] + sys.path
print (sys.path)
import imp
foo = imp.load_source('component', '/home/mzarowiecki/git/bertha/code/bertha-compute/analysis/pipeline/components/vcf_cleanup_somatic/component.py')
from component import main


# Get the right inputfile
#"files": {"vcf_file": "/genomes/bertha-prod/analysis/singlesample/LP3000396-DNA_A01/CANCP41877/normalise_vcf/116332/1/output/LP3000396-DNA_A01.duprem.left.split.vcf.gz"
params['input_files']['files']['vcf_file'] = '/'.join((od,  ''.join((params['sample_id'][1],'.duprem.left.split.vcf.gz')) ))


params['output_files']['dirs']['log_dir']='/'.join((here, params['sample_id'][1], 'vcf_cleanup_somatic', '12345', '2', 'logs' ))
params['output_files']['dirs']['output_dir']='/'.join((here, params['sample_id'][1], 'vcf_cleanup_somatic', '12345', '2', 'output' ))

logfol =params['output_files']['dirs']['log_dir'] 
od =params['output_files']['dirs']['output_dir'] 

if not os.path.isdir(params['output_files']['dirs']['base_dir']):
    os.makedirs(params['output_files']['dirs']['base_dir'])


if not os.path.isdir(params['output_files']['dirs']['log_dir']):
    os.makedirs(params['output_files']['dirs']['log_dir'])


if not os.path.isdir(params['output_files']['dirs']['output_dir']):
    os.makedirs(params['output_files']['dirs']['output_dir'])










# And then run it twice, once for normal variants, and once for SV

# Run vcf_cleanup_somatic - NOT IMPLEMENTED
# "files": {"vcf_file": "/genomes/bertha-prod/by_name/LP3000396-DNA_A01/CANCP41877/SomaticVariations/LP3000417-DNA_A01_LP3000396-DNA_A01.somatic.SV.vcf.gz"
# main(params)


# Run vcf_cleanup_somatic
main(params)

files = os.listdir(here)

for f in files:
    if (f.startswith( ''.join((params['sample_id'][1] , '.' )) )):
        shutil.move(f, od)

print ("Finished successfully!!!\n\nDont worry about the error messages\n\n")



