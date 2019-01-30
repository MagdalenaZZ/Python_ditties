"""
+---------------------------+----------------------------------------------------------------------+
| Display Name	            | Somatic VC                                                           |
+---------------------------+----------------------------------------------------------------------+
| Short Name                | somatic_vc                                                           |
+---------------------------+----------------------------------------------------------------------+
| Stage                     | Cancer variant calling                                               |
+---------------------------+----------------------------------------------------------------------+
| Since Bertha Version      | 2.x                                                                  |
+---------------------------+----------------------------------------------------------------------+
| Sample Types	            | Cancer                                                               |
+---------------------------+----------------------------------------------------------------------+
| Component Type	    | multi-sample                                                         |
+---------------------------+----------------------------------------------------------------------+
| Component Parallelisation | single-application                                                   |
+---------------------------+----------------------------------------------------------------------+


Description
===========

Run Strelka for somatic variant calls 

Github Link
https://github.com/Illumina/strelka

Version History
===============

+---------+-------------+----------------------------------------------------------------+---------+
| Version | Date        | Comments                                                       | Since   |
|         |             |                                                                | Bertha  |
|         |             |                                                                | Version |
+=========+=============+================================================================+=========+
| 1.0.0   | 17/12/2018  | - initial version                                              | NA      |
+---------+-------------+----------------------------------------------------------------+---------+


External Tools
==============

This component does not use any external tools.

+----------------+----------------------------------------------------------------------+----------+
| Name           | Description                                                          | Version  |
+================+======================================================================+==========+
|  Strelka       |   module add strelka                                                 | 2.9.9    |
+----------------+----------------------------------------------------------------------+----------+



"""

from __future__ import print_function
import sys
import shutil
import os.path
import re
from subprocess import call


# Bertha components import

VERSION = '0.0.1'

#from bertha.execution.utils.helpers import (
#    timed_command,
#    check_file_exists_and_nonempty
#)

### Independent runnable section #########################################



import argparse



# Describe what the script does
parser = argparse.ArgumentParser(description='This script writes the commands for running Strelka somatic VC', formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-n', '--normal_bam', default=None, dest='norm', action='store', required=True, help="Normal.bam file full path")
parser.add_argument('-t', '--tumour_bam', default=None, dest='tum', action='store', required=True, help="Tumour.bam file full path")
parser.add_argument('-g', '--genome', default=None, dest='gfa', action='store', required=True, help="genome.fa")
parser.add_argument('-o', '--out_folder', default=None, dest='out', action='store', required=True, help="Output folder full path")
parser.add_argument('-c', '--bed_boundaries', default=None, dest='bed', action='store', required=True, help="BED file with VC calling boundaries")



# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()


######### Job submission and environment ###########################################

#   -n 12 -R "span[hosts=1]" -R "select[largedsk]" -R "hname!='hpc-prod-grid-lsfexec-001' && hname!='hpc-prod-grid-lsfexec-002'"
# module add strelka   # currently strelka/2.9.9
execfile('/usr/local/Modules/3.2.10/init/python.py')
module('load','strelka')


######### Creating params ###########################################


sample_id=args.out.split("/")[-1]

params= {}
params['config']={}
params['config']['strelka']='configureStrelkaSomaticWorkflow.py'
params['config']['bcftools']='bcftools'
params['input_files']={}
params['input_files']['files']={}
params['input_files']['files']['normal_bam_file']=args.norm
params['input_files']['files']['tumour_bam_file']=args.tum
params['input_files']['files']['reference_genome_fasta']=args.gfa
params['input_files']['files']['regions_bed']=args.bed
params['output_files']={}
params['output_files']['folders']={} 
params['output_files']['folders']['out_folder']=args.out
params['config']['strelka_runscript']= args.out + '/runWorkflow.py'
params['home']='/home/mzarowiecki/scratch/New_deliveries/genomes/by_name'
params['output_files']['files']={} 

# Output files generated by strelka
params['output_files']['dirs']={}
params['output_files']['dirs']['base_dir']=args.out



params['output_files']['files']['callable'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/regions/somatic.callable.regions.bed.gz' )
params['output_files']['files']['callable_index'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/regions/somatic.callable.regions.bed.gz.tbi' )
params['output_files']['files']['somatic_vcf_snv'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/variants', 'somatic.snvs.vcf.gz' )
params['output_files']['files']['somatic_vcf_snv_index'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/variants', 'somatic.snvs.vcf.gz.tbi' )
params['output_files']['files']['somatic_vcf_indel'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/variants', 'somatic.indels.vcf.gz' )
params['output_files']['files']['somatic_vcf_indel_index'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/variants', 'somatic.indels.vcf.gz.tbi' )
# Output files generated by this script
params['output_files']['files']['somatic_vcf'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/variants', ''.join([sample_id,'.somatic.vcf.gz']) )
params['output_files']['files']['somatic_vcf_index'] = os.path.join(params['output_files']['folders']['out_folder'], 'results/variants', ''.join([sample_id,'.somatic.vcf.gz.tbi']) )
# Core "delivery" folder which final results will be moved to
params['output_files']['folders']['core_folder']=os.path.join(params['home'],sample_id)

# Fake error and log files
params['output_files']['logfiles']={} 
params['output_files']['logfiles']['strelka_log']=os.path.join(params['output_files']['folders']['core_folder'],'strelka.log')
params['output_files']['logfiles']['strelka_config_log']=os.path.join(params['output_files']['folders']['core_folder'],'strelka_config.log')
params['output_files']['logfiles']['bcftools_warning']=os.path.join(params['output_files']['folders']['core_folder'],'bcftools_merge.log')








try:
    os.makedirs(params['output_files']['folders']['core_folder'])
    open('w', params['output_files']['logfiles']['strelka_log']).close()
    open('w', params['output_files']['logfiles']['strelka_config_log']).close()
    open('w', params['output_files']['logfiles']['bcftools_warning']).close()

except:
    pass


#print (params)

#for file_name in params['output_files']['files']:
#    print (file_name)

#for filetype,file_name in params['output_files']['files'].iteritems():
#    print (file_name, filetype)


####################################



import imp
sample = imp.load_source('run_strelka', '/home/mzarowiecki/git/mertha/somatic_vc/component.py')
sample.main(params)





