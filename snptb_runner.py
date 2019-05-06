#!/usr/bin/python3
#!/usr/bin/env python

import os
import sys
import re
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, execute_subprocess, \
    extract_read_list, file_to_list
from bbduk_trimmer import bbduk_trimming
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_recall import picard_dictionary, samtools_faidx, picard_markdup, haplotype_caller, call_variants, \
    select_variants, hard_filter, combine_gvcf

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

#Check if there are samples to filter
sample_list_F = []
if args.sample_list == None:
    print("No samples to filter")
    for r1_file, r2_file in zip(r1, r2):
        sample = extract_sample(r1_file, r2_file)
        sample_list_F.append(sample)
else:
    print("samples will be filtered")
    sample_list_F = file_to_list(args.sample_list)
print("%d samples will be analysed: %s" % (len(sample_list_F), ",".join(sample_list_F)))


######################################################################
#####################START PIPELINE###################################
######################################################################

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    if sample in sample_list_F:
        args.r1_file = r1_file
        args.r2_file = r2_file

        print("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + END_FORMATTING)

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
            bbduk_trimming(args)

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
            picard_markdup(args)
        
        #HAPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
        #######################################################
        out_gvcfr_dir = os.path.join(args.output, "GVCF_recal")
        out_gvcfr_name = sample + ".g.vcf"
        output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

        if os.path.isfile(output_gvcfr_file):
            print(YELLOW + DIM + output_gvcfr_file + BOLD + " EXIST\nOmmiting Haplotype Call (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Haplotype Calling (Recall) in sample " + sample + END_FORMATTING)
            haplotype_caller(args, recalibrate=True, ploidy=2, bamout=False, forceactive=False)

        ###############################################################################################################################################
        #############################FOR COMPARING PURPOSE#############################################################################################
        ###############################################################################################################################################

        #CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
        #######################################################
        out_vcfr_dir = os.path.join(args.output, "VCF_recal")
        out_vcfr_name = sample + ".raw.vcf"
        output_vcfr_file = os.path.join(out_vcfr_dir, out_vcfr_name)

        if os.path.isfile(output_vcfr_file):
            print(YELLOW + DIM + output_vcfr_file + BOLD + " EXIST\nOmmiting Variant Calling (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Variant Calling (Recall) in sample " + sample + END_FORMATTING)
            call_variants(args, recalibrate=True, group=False)

        #SELECT VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
        #########################################################
        out_vcfsnpr_name = sample + ".snp.vcf"
        output_vcfsnpr_file = os.path.join(out_vcfr_dir, out_vcfsnpr_name)

        if os.path.isfile(output_vcfsnpr_file):
            print(YELLOW + DIM + output_vcfsnpr_file + BOLD + " EXIST\nOmmiting Variant Selection (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Selecting Variants (Recall) in sample " + sample + END_FORMATTING)
            select_variants(output_vcfr_file, select_type='SNP') #select_variants(output_vcfr_file, select_type='INDEL')
        
        #HARD FILTER VARIANTS 1/2 FOR RECALIBRATION #############
        #########################################################
        out_vcfhfsnpr_name = sample + ".snp.hf.vcf"
        output_vcfhfsnpr_file = os.path.join(out_vcfr_dir, out_vcfhfsnpr_name)

        if os.path.isfile(output_vcfhfsnpr_file):
            print(YELLOW + DIM + output_vcfhfsnpr_file + BOLD + " EXIST\nOmmiting Hard Filtering (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Hard Filtering Variants (Recall) in sample " + sample + END_FORMATTING)
            hard_filter(output_vcfsnpr_file, select_type='SNP')

#ONCE ALL GVCF VARIANTS ARE CALLED, THEY ARE GATHERED AND FILTERED 
# TO RECALIBRATE ORIGINAL MARKDUPPED BAM
######################################################################
##############START GROUP CALLING FOR RECALIBRATION###################
######################################################################
group_name = args.output.split("/")[-1]
print("\n" + WHITE_BG + "STARTING JOINT CALL FOR RECALIBATION IN GROUP: " + group_name + END_FORMATTING)

#CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_gvcfr_dir = os.path.join(args.output, "GVCF_recal")
out_gvcfr_name = group_name + ".cohort.g.vcf"
output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

if os.path.isfile(output_gvcfr_file):
    print(YELLOW + DIM + output_gvcfr_file + BOLD + " EXIST\nOmmiting GVCF Combination (Recall) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "GVCF Combination (Recall) in group " + group_name + END_FORMATTING)
    combine_gvcf(args, recalibrate=True)

#CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_vcfr_dir = os.path.join(args.output, "VCF_recal")
out_vcfr_name = group_name + ".cohort.raw.vcf"
output_vcfr_file = os.path.join(out_vcfr_dir, out_vcfr_name)

if os.path.isfile(output_vcfr_file):
    print(YELLOW + DIM + output_vcfr_file + BOLD + " EXIST\nOmmiting Variant Calling (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Variant Calling (Recall-Group) in group " + group_name + END_FORMATTING)
    call_variants(args, recalibrate=True, group=True)

#SELECT VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#########################################################
out_vcfsnpr_name = group_name + ".cohort.snp.vcf"
out_vcfindelr_name = group_name + ".cohort.indel.vcf"
output_vcfsnpr_file = os.path.join(out_vcfr_dir, out_vcfsnpr_name)
output_vcfindelr_file = os.path.join(out_vcfr_dir, out_vcfindelr_name)

if os.path.isfile(output_vcfsnpr_file) and os.path.isfile(output_vcfindelr_file):
    print(YELLOW + DIM + output_vcfsnpr_file + BOLD + " EXIST\nOmmiting Variant Selection (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Selecting Variants (Recall-Group) in group " + group_name + END_FORMATTING)
    select_variants(output_vcfr_file, select_type='SNP') #select_variants(output_vcfr_file, select_type='INDEL')
    select_variants(output_vcfr_file, select_type='INDEL')

#HARD FILTER VARIANTS 1/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnpr_name = group_name + ".cohort.snp.hf.vcf"
out_vcfhfindelr_name = group_name + ".cohort.indel.hf.vcf"
output_vcfhfsnpr_file = os.path.join(out_vcfr_dir, out_vcfhfsnpr_name)
output_vcfhfindelr_file = os.path.join(out_vcfr_dir, out_vcfhfindelr_name)


if os.path.isfile(output_vcfhfsnpr_file) and os.path.isfile(output_vcfhfindelr_file):
    print(YELLOW + DIM + output_vcfhfsnpr_file + BOLD + " EXIST\nOmmiting Hard Filtering (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Hard Filtering Variants (Recall-Group) in group " + group_name + END_FORMATTING)
    hard_filter(output_vcfsnpr_file, select_type='SNP')
    hard_filter(output_vcfindelr_file, select_type='INDEL')


    ######################################################################
    ##############START RECALIBRATION AND FINAL CALL######################
    ######################################################################

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    if sample in sample_list_F:
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

"""

#./snptb_runner.py -i /home/laura/ANALYSIS/Lofreq/coinfection_designed/raw -r reference/MTB_ancestorII_reference.fasta -o /home/laura/ANALYSIS/Lofreq/coinfection_designed/TEST -s sample_list.txt