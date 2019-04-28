#!/usr/bin/python3
#!/usr/bin/env python

import os
import sys
import re
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, execute_subprocess, extract_read_list
from bbduk_trimmer import bbduk_trimming

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 28 April 2019
REVISION: 

TODO: Use a list of samples to filter 
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def get_arguments():

    parser = argparse.ArgumentParser(prog = 'snptb.py', description= 'Pipeline to call variants (SNVs) with any non model organism. Specialised in Mycobacterium Tuberculosis')
    
    input_group = parser.add_argument_group('Input', 'Required fiparameters')

    input_group.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='Input directory containing all fast[aq] files')

    output_group = parser.add_argument_group('Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True, help='Output directory to extract all results')
    output_group.add_argument('-s', '--sample_list', type=str, required=False, help='Sample name to handle output files ')

    trimming_group = parser.add_argument_group('Trimming parameters', 'parameters for diferent triming conditions')

    trimming_group.add_argument('-H', '--hdist', type=str, required=False, help='Set hdist parameter, default 2')
    trimming_group.add_argument('-k', '--kmer', type=str, required=False, help='Set k parameter, default 21')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    
    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=4, help='Threads to use')
    params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=8, help='MAx memory to use')



    #argcomplete.autocomplete(parser)
    arguments = parser.parse_args()

    return arguments

args = get_arguments()

#Obtain all R1 and R2 from folder
r1, r2 = extract_read_list(args.input_dir)

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    args.r1_file = r1_file
    args.r2_file = r2_file

    print(WHITE_BG + "STARTING SAMPLE " + sample + END_FORMATTING)

    ##############START PIPELINE#####################
    #################################################


    #INPUT ARGUMENTS
    ################
    check_file_exists(args.r1_file)
    check_file_exists(args.r2_file)

    args.output = os.path.abspath(args.output)
    #QUALITY CHECK
    ##############
    """
    TODO
    """
    
    #QUALITY TRIMMING AND ADAPTER REMOVAL WITH bbduk.sh
    ###################################################
    out_trim_dir = os.path.join(args.output, "Trimmed")
    out_trim_name = sample + "_R1.clean.fastq.gz"
    output_trimming_file = os.path.join(out_trim_dir, out_trim_name)

    if os.path.isfile(output_trimming_file):
        print(YELLOW + output_trimming_file + BOLD + " EXIST, ommiting Trimming fot sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Trimming sample " + sample + END_FORMATTING)
        bbduk_trimming(args)
