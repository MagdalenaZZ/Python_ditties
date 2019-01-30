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

######### Creating params ###########################################


sample_id=args.out.split("/")[-1]

params= {}
params['config']={}
params['config']['strelka']='configureStrelkaSomaticWorkflow.py'

params['input_files']={}
params['input_files']['files']={}
params['input_files']['files']['normal_bam_file']=args.norm
params['input_files']['files']['tumour_bam_file']=args.tum
params['input_files']['files']['reference_genome_fasta']=args.gfa
params['input_files']['files']['regions_bed']=args.bed
params['output_files']={}
params['output_files']['folders']={} 
params['output_files']['folders']['out_folder']=args.out
params['config']['strelka_runscript']= args.out + 'runWorkflow.py'
params['home']='/home/mzarowiecki/scratch/New_deliveries/genomes/by_name'
params['output_files']['files']={} 

# Output files generated by strelka
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


try:
    os.makedirs(params['output_files']['folders']['core_folder'])
except:
    pass


# print (params)

#for file_name in params['output_files']['files']:
#    print (file_name)

#for filetype,file_name in params['output_files']['files'].iteritems():
#    print (file_name, filetype)


####################################
import imp
sample = imp.load_source('run_strelka', '/home/mzarowiecki/bin/Python_ditties/run_strelka.py')

sample.main(params)




def main(params):
    
    """
    Executes the component tasks.

    :param params: dictionary of component parameters
    :return: list of timed command profiles and  list of files to register in catalog in addition
             to the files in the params dictionary
    :rtype: tuple
    """
    
    profiles = []

    # Call strelka to initiate the run, and create a python run-script
    profile, _ = timed_command(
        '{strelka} --normalBam={normal_bam_file} --tumorBam={tumour_bam_file} --outputCallableRegions --referenceFasta={genome_ref} --callRegions={call_regions} --runDir={out_folder} ' 
        '> {strelka_all}'.format(
            strelka=params['config']['strelka'],
            normal_bam_file=params['input_files']['files']['normal_bam_file'],
            tumour_bam_file=params['input_files']['files']['tumour_bam_file'],
            genome_ref=params['input_files']['files']['reference_genome_fasta'],
            call_regions=params['input_files']['files']['regions_bed'],
            out_folder=params['output_files']['folders']['out_folder'])
    )

    print (         
        '{strelka} --normalBam={normal_bam_file} --tumorBam={tumour_bam_file} --outputCallableRegions --referenceFasta={genome_ref} --callRegions={call_regions} --runDir={out_folder} ' 
        '> {strelka_all}'.format(
            strelka=params['config']['strelka'],
            normal_bam_file=params['input_files']['files']['normal_bam_file'],
            tumour_bam_file=params['input_files']['files']['tumour_bam_file'],
            genome_ref=params['input_files']['files']['reference_genome_fasta'],
            call_regions=params['input_files']['files']['regions_bed'],
            out_folder=params['output_files']['folders']['out_folder'])
    )
    print ("Hello")

    # Example:    configureStrelkaSomaticWorkflow.py --normalBam=/genomes/analysis/cancer/new_pipeline/ISO_validation/1_mapping/LP3000417-DNA_H02/LP3000417-DNA_H02.bam --tumorBam=/genomes/analysis/cancer/new_pipeline/ISO_validation/1_mapping/LP3000396-DNA_G02/LP3000396-DNA_G02.bam --outputCallableRegions --referenceFasta=/genomes/analysis/cancer/new_pipeline/Dragen/HASH_tables/GRCh38_full_analysis_set_plus_decoy_hla.01.011.264.3.1.7.fa --callRegions=/genomes/analysis/cancer/new_pipeline/Dragen/HASH_tables/GRCh38_full_analysis_set_plus_decoy_hla.01.011.269.3.2.3.noHLAdecoy.bed.gz --runDir=/genomes/analysis/cancer/new_pipeline/ISO_validation/2_SomaticVC/LP3000396-DNA_G02

    # Do the strelka run. 
    # -m local makes it run as a local machine, rather than as an SGE cluster (no LSF option available)
    # -j 12 uses all 12 nodes you specified in the bsub command

    profile, _ = timed_command(
        '{python} {run_script}  -m local -j 12'
        '> {run_strelka_all}'.format(
            python='/usr/bin/python',
            run_script=params['config']['strelka_runscript'])
    )

     # Example: bsub -P Analysis -J tLP3000396-DNA_F01 -q bio -o tLP3000396-DNA_F01.o -e tLP3000396-DNA_F01.e -n 12 -R "span[hosts=1]" -R "select[largedsk]" -R "hname!='hpc-prod-grid-lsfexec-001' && hname!='hpc-prod-grid-lsfexec-002'" python /genomes/analysis/cancer/new_pipeline/ISO_validation/2_VC/LP3000396-DNA_F01/runWorkflow.py -m local -j 12
     


    # Merge strelka snvs and indels into one sorted VCF file
    merge_strelka_output(params['output_files']['folders']['out_folder'])

    # Tidy up the strelka output folder, and put into sample base structure
    tidy_up_strelka_output(params['output_files']['folders']['out_folder'])



def merge_strelka_output(out_folder):
    """
    Merges the strelka output vcf files and creates index

    """

    # concatenate VCF files
    self._timed_command(
        '{bcftools} concat -a -O gz -o {out_file} {vcf_snv} {vcf_indel} 2>>{bcftools_warning}'.format(
            bcftools=params['config']['bcftools'],
            vcf_snv=params['output_files']['files']['somatic_vcf_snv'],
            vcf_indel=params['output_files']['files']['somatic_vcf_indel'], 
            out_file=params['output_files']['files']['somatic_vcf'],
            bcftools_warning=params['output_files']['files']['bcftools_warning']),
        params['output_files']['files']['bcftools_warning']
    )

    # bcftools concat -a -O gz -o merge.vcf.gz LP3000592-DNA_H09/results/variants/somatic.*.vcf.gz



def tidy_up_strelka_output(out_folder):
    """
    Cleans up the strelka output and keeps only necessary files

    """

    for filetype,file_name in params['output_files']['files'].iteritems():

        if os.path.exists(file_name):
            dest = params['output_files']['folders']['core_folder']
            shutil.move(file_name, dest)
        else:
            missing_files.append(dest)

        # If complete strelka output was not generated, complain
        if missing_files:
            error_message = 'The following files were not generated: {}'.format(missing_files)
            logger.error(error_message)
            raise ProcessingException(error_message)

    shutil.rmtree(out_folder)









# Example function, completely redundant

def make_basedir_writable(vcf_base_dir):
    """
    Makes the parent directory of the VCF file writeable.

    :param vcf_base_dir: the VCF parent dir
    """
    logger.info('Making {0} writable'.format(vcf_base_dir))
    os.chmod(vcf_base_dir, READ_WRITE)
















