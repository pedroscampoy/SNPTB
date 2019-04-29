#!/usr/bin/python3
#!/usr/bin/env python

import os
import sys
import re
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, execute_subprocess, \
    extract_read_list
from bbduk_trimmer import bbduk_trimming
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_recall import picard_dictionary, samtools_faidx, picard_markdup, haplotype_caller, call_variants, \
    select_variants, hard_filter

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
    input_group.add_argument('-r', '--reference', metavar="reference", type=str, required=True, help='File to map against')

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

    print("STARTING SAMPLE " + WHITE_BG + sample + END_FORMATTING)

    ##############START PIPELINE#####################
    #################################################


    #INPUT ARGUMENTS
    ################
    check_file_exists(args.r1_file)
    check_file_exists(args.r2_file)

    args.output = os.path.abspath(args.output)
    check_create_dir(args.output)
    #QUALITY CHECK
    ##############
    """
    TODO: Quality check 
    """
    
    #QUALITY TRIMMING AND ADAPTER REMOVAL WITH bbduk.sh
    ###################################################
    out_trim_dir = os.path.join(args.output, "Trimmed")
    out_trim_name_r1 = sample + "_R1.clean.fastq.gz"
    out_trim_name_r2 = sample + "_R2.clean.fastq.gz"
    output_trimming_file_r1 = os.path.join(out_trim_dir, out_trim_name_r1)
    output_trimming_file_r2 = os.path.join(out_trim_dir, out_trim_name_r2)
    
    if os.path.isfile(output_trimming_file_r1) and os.path.isfile(output_trimming_file_r2):
        print(YELLOW + DIM + output_trimming_file_r1 + BOLD + " EXIST\nOmmiting Trimming for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Trimming sample " + sample + END_FORMATTING)
        #bbduk_trimming(args)

    #MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
    #####################################################
    out_map_dir = os.path.join(args.output, "Bam")
    out_map_name = sample + ".rg.sorted.bam"
    output_map_file = os.path.join(out_map_dir, out_map_name)

    args.r1_file = output_trimming_file_r1
    args.r2_file = output_trimming_file_r2

    if os.path.isfile(output_map_file):
        print(YELLOW + DIM + output_map_file + BOLD + " EXIST\nOmmiting Mapping for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Mapping sample " + sample + END_FORMATTING)
        print("R1: " + output_trimming_file_r1 + "\nR2: " + output_trimming_file_r2 + "\nReference: " + args.reference)
        bwa_mapping(args)
        sam_to_index_bam(args)

    #PREPARE REFERENCE FOR MAPPING + FAI + DICT #########
    #####################################################

    picard_dictionary(args)
    samtools_faidx(args)

    #MARK DUPLICATES WITH PICARDTOOLS ###################
    #####################################################
    #TO DO: remove output_map_file and include markdup in previous step checking for existence of .rg.markdup.sorted.bam
    #out_markdup_dir = os.path.join(args.output, "Bam")
    out_markdup_name = sample + ".rg.markdup.sorted.bam"
    output_markdup_file = os.path.join(out_map_dir, out_markdup_name)

    args.input_bam = output_map_file

    if os.path.isfile(output_markdup_file):
        print(YELLOW + DIM + output_markdup_file + BOLD + " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Marking Dupes in sample " + sample + END_FORMATTING)
        print("Input Bam: " + args.input_bam)
        #picard_markdup(args)
    
    #HPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
    #######################################################
    out_gvcfr_dir = os.path.join(args.output, "GVCF_recal")
    out_markdup_name = sample + ".rg.markdup.sorted.bam"
    output_markdup_file = os.path.join(out_markdup_dir, out_markdup_name)

    args.input_bam = output_map_file

    if os.path.isfile(output_markdup_file):
        print(YELLOW + DIM + output_markdup_file + BOLD + " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Marking Dupes in sample " + sample + END_FORMATTING)
        print("Input Bam: " + args.input_bam)
        #picard_markdup(args)


    #haplotype_caller(args, recalibrate=True, ploidy=1, bamout=False, forceactive=False)
    #call_variants(args, recalibrate=True)

    out = args.output
    sample = args.sample
    #raw_vcf = out + "/VCF_recal/" + sample + ".raw.vcf"
    #select_variants(raw_vcf, select_type='SNP') #select_variants(raw_vcf, select_type='INDEL')

    selected_vcf = out + "/VCF_recal/" + sample + ".snp.vcf"
    hard_filter(selected_vcf, select_type='SNP')