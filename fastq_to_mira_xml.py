#!/usr/bin/env python3.3

import argparse
import fastn

parser = argparse.ArgumentParser(
    description = 'Creates an xml file from a fastq file of reads, for use with Mira assembler',
    usage = '%(prog)s [options] <fastq_in> <xml_out>')
parser.add_argument('fastq_in', help='Name of input fastq file')
parser.add_argument('xml_out', help='Name of output xml file')
options = parser.parse_args()
fastn.fastq_to_mira_xml(options.fastq_in, options.xml_out)
