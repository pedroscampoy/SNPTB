#!/usr/bin/env python

import os
import sys
import re
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, execute_subprocess, \
    extract_read_list, file_to_list, get_coverage, obtain_group_cov_stats, remove_low_covered
from bbduk_trimmer import bbduk_trimming
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_recall import picard_dictionary, samtools_faidx, picard_markdup, haplotype_caller, call_variants, \
    select_variants, hard_filter, combine_gvcf, select_pass, select_pass_variants, recalibrate_bam, \
    split_vcf_saples
from vcf_process import vcf_consensus_filter

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
    
    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all fast[aq] files')
    input_group.add_argument('-r', '--reference', metavar="reference", type=str, required=True, help='REQUIRED. File to map against')
    input_group.add_argument('-s', '--sample', metavar="sample", type=str, required=False, help='Sample to identify further files')

    output_group = parser.add_argument_group('Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')
    output_group.add_argument('-S', '--sample_list', type=str, required=False, help='Sample names to analyse only in the file supplied')

    trimming_group = parser.add_argument_group('Trimming parameters', 'parameters for diferent triming conditions')

    trimming_group.add_argument('-H', '--hdist', type=str, required=False, help='Set hdist parameter, default 2')
    trimming_group.add_argument('-k', '--kmer', type=str, required=False, help='Set k parameter, default 21')

    gatk_group = parser.add_argument_group('GATK parameters', 'parameters for diferent variant calling')

    gatk_group.add_argument('-p', '--ploidy', type=str, required=False, default=2, help='Set ploidy when HC call, default 2')
    gatk_group.add_argument('-E', '--enrich_gvcf', required=False,  default=False, help='Point a directory with g.vcf files to enrich the analysis')


    vcf_group = parser.add_argument_group('VCF filters', 'parameters for variant filtering')

    vcf_group.add_argument('-m', '--maxnocallfr', type=str, required=False, default=0.2, help='maximun proportion of samples with non genotyped alleles')

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
    print("\n" + "No samples to filter")
    for r1_file, r2_file in zip(r1, r2):
        sample = extract_sample(r1_file, r2_file)
        sample_list_F.append(sample)
else:
    print("samples will be filtered")
    sample_list_F = file_to_list(args.sample_list)
print("\n%d samples will be analysed: %s" % (len(sample_list_F), ",".join(sample_list_F)))


######################################################################
#####################START PIPELINE###################################
######################################################################
output = os.path.abspath(args.output)
group_name = output.split("/")[-1]

print("\n\n" + BLUE + BOLD + "STARTING PIPELINE IN GROUP: " + group_name + END_FORMATTING)


#PREPARE REFERENCE FOR MAPPING + FAI + DICT #########
#####################################################

picard_dictionary(args)
samtools_faidx(args)


for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    args.sample = sample
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
            print(YELLOW + DIM + output_trimming_file_r1 + " EXIST\nOmmiting Trimming for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Trimming sample " + sample + END_FORMATTING)
            bbduk_trimming(args)

        #MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
        #####################################################
        out_map_dir = os.path.join(args.output, "Bam")
        out_map_name = sample + ".rg.sorted.bam"
        output_map_file = os.path.join(out_map_dir, out_map_name)

        out_markdup_name = sample + ".rg.markdup.sorted.bam"
        output_markdup_file = os.path.join(out_map_dir, out_markdup_name)

        args.r1_file = output_trimming_file_r1
        args.r2_file = output_trimming_file_r2

        if os.path.isfile(output_map_file) or os.path.isfile(output_markdup_file):
            print(YELLOW + DIM + output_map_file + " EXIST\nOmmiting Mapping for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Mapping sample " + sample + END_FORMATTING)
            print("R1: " + output_trimming_file_r1 + "\nR2: " + output_trimming_file_r2 + "\nReference: " + args.reference)
            bwa_mapping(args)
            sam_to_index_bam(args)

        #MARK DUPLICATES WITH PICARDTOOLS ###################
        #####################################################
        #TO DO: remove output_map_file and include markdup in previous step checking for existence of .rg.markdup.sorted.bam
        #out_markdup_dir = os.path.join(args.output, "Bam")
        out_markdup_name = sample + ".rg.markdup.sorted.bam"
        output_markdup_file = os.path.join(out_map_dir, out_markdup_name)

        args.input_bam = output_map_file

        if os.path.isfile(output_markdup_file):
            print(YELLOW + DIM + output_markdup_file + " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Marking Dupes in sample " + sample + END_FORMATTING)
            print("Input Bam: " + args.input_bam)
            picard_markdup(args)
        
        #CALCULATE COVERAGE FOR EACH POSITION##################
        #######################################################
        out_cov_dir = os.path.join(args.output, "Coverage")
        out_cov_name = sample + ".cov"
        output_cov_file = os.path.join(out_cov_dir, out_cov_name)

        if os.path.isfile(output_cov_file):
            print(YELLOW + DIM + output_cov_file + " EXIST\nOmmiting coverage calculation for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Calculating coverage in sample " + sample + END_FORMATTING)
            get_coverage(args, output_markdup_file, output_fmt="-d")
        
        #HAPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
        #######################################################
        out_gvcfr_dir = os.path.join(args.output, "GVCF_recal")
        out_gvcfr_name = sample + ".g.vcf"
        output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

        if os.path.isfile(output_gvcfr_file):
            print(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting Haplotype Call (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Haplotype Calling (Recall) in sample " + sample + END_FORMATTING)
            haplotype_caller(args, recalibrate=True, ploidy=args.ploidy, bamout=False, forceactive=False)
"""
        ###############################################################################################################################################
        #############################FOR COMPARING PURPOSE#############################################################################################
        ###############################################################################################################################################

        #CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
        #######################################################
        out_vcfr_dir = os.path.join(args.output, "VCF_recal")
        out_vcfr_name = sample + ".raw.vcf"
        output_vcfr_file = os.path.join(out_vcfr_dir, out_vcfr_name)

        if os.path.isfile(output_vcfr_file):
            print(YELLOW + DIM + output_vcfr_file + " EXIST\nOmmiting Variant Calling (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Variant Calling (Recall) in sample " + sample + END_FORMATTING)
            call_variants(args, recalibrate=True, group=False)

        #SELECT VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
        #########################################################
        out_vcfsnpr_name = sample + ".snp.vcf"
        output_vcfsnpr_file = os.path.join(out_vcfr_dir, out_vcfsnpr_name)

        if os.path.isfile(output_vcfsnpr_file):
            print(YELLOW + DIM + output_vcfsnpr_file + " EXIST\nOmmiting Variant Selection (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Selecting Variants (Recall) in sample " + sample + END_FORMATTING)
            select_variants(output_vcfr_file, select_type='SNP') #select_variants(output_vcfr_file, select_type='INDEL')
        
        #HARD FILTER VARIANTS 1/2 FOR RECALIBRATION #############
        #########################################################
        out_vcfhfsnpr_name = sample + ".snp.hf.vcf"
        output_vcfhfsnpr_file = os.path.join(out_vcfr_dir, out_vcfhfsnpr_name)

        if os.path.isfile(output_vcfhfsnpr_file):
            print(YELLOW + DIM + output_vcfhfsnpr_file + " EXIST\nOmmiting Hard Filtering (Recall) for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Hard Filtering Variants (Recall) in sample " + sample + END_FORMATTING)
            hard_filter(output_vcfsnpr_file, select_type='SNP')
"""



group_name = output.split("/")[-1]
print("\n\n" + BLUE + BOLD + "CHECKING LOW COVERED SAMPLES IN GROUP: " + group_name + END_FORMATTING + "\n")

#GROUP COVERAGE SUMMARY STATS##########################
#######################################################

out_covg_dir = os.path.join(args.output, "Coverage")
out_covg_name = group_name + ".covegare.tab"
output_covg_file = os.path.join(out_covg_dir, out_covg_name)
"""
if os.path.isfile(output_covg_file):
    print(YELLOW + DIM + output_covg_file + " EXIST\nOmmiting group coverage calculation for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Group coverage stats in group " + group_name + END_FORMATTING)
    saples_low_covered = obtain_group_cov_stats(out_covg_dir, low_cov_threshold=20, unmmaped_threshold=20)
"""

saples_low_covered = obtain_group_cov_stats(out_covg_dir, low_cov_threshold=20, unmmaped_threshold=20)

if len(saples_low_covered) > 0:
    print("\n" + YELLOW + BOLD + "There are sample(s) with low coverage that will be removed from the analysis: " + "\n"\
         + ",".join(saples_low_covered) + END_FORMATTING + "\n")
    remove_low_covered(args.output, saples_low_covered)
    #Remove sample from the list of filtered samples
    ################################################
    for samples_to_remove in saples_low_covered:
        sample_list_F.remove(samples_to_remove)
else:
    print("\n" + YELLOW + BOLD + "All samples have a decent depth of coverage according to threshold supplied" + "\n")




#ONCE ALL GVCF VARIANTS ARE CALLED, THEY ARE GATHERED AND FILTERED 
# TO RECALIBRATE ORIGINAL MARKDUPPED BAM
######################################################################
##############START GROUP CALLING FOR RECALIBRATION###################
######################################################################

group_name = output.split("/")[-1]
print("\n\n" + BLUE + BOLD + "STARTING JOINT CALL FOR RECALIBATION IN GROUP: " + group_name + END_FORMATTING + "\n")

#CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_gvcfr_dir = os.path.join(args.output, "GVCF_recal")
out_gvcfr_name = group_name + ".cohort.g.vcf"
output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

if os.path.isfile(output_gvcfr_file):
    print(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting GVCF Combination (Recall) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "GVCF Combination (Recall) in group " + group_name + END_FORMATTING)
    combine_gvcf(args, recalibrate=True, all_gvcf=args.enrich_gvcf)

#CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_vcfr_dir = os.path.join(args.output, "VCF_recal")
out_vcfr_name = group_name + ".cohort.raw.vcf"
output_vcfr_file = os.path.join(out_vcfr_dir, out_vcfr_name)

if os.path.isfile(output_vcfr_file):
    print(YELLOW + DIM + output_vcfr_file + " EXIST\nOmmiting Variant Calling (Recall-Group) for group " + group_name + END_FORMATTING)
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
    print(YELLOW + DIM + output_vcfsnpr_file + " EXIST\nOmmiting Variant Selection (Recall-Group) for group " + group_name + END_FORMATTING)
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
    print(YELLOW + DIM + output_vcfhfsnpr_file + " EXIST\nOmmiting Hard Filtering (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Hard Filtering Variants (Recall-Group) in group " + group_name + END_FORMATTING)
    hard_filter(output_vcfsnpr_file, select_type='SNP')
    hard_filter(output_vcfindelr_file, select_type='INDEL')

#PASS FILTER VARIANTS 1/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnppass_name = group_name + ".cohort.snp.hf.pass.vcf"
out_vcfhfindelpass_name = group_name + ".cohort.indel.hf.pass.vcf"
output_vcfhfsnppass_file = os.path.join(out_vcfr_dir, out_vcfhfsnppass_name)
output_vcfhfindelpass_file = os.path.join(out_vcfr_dir, out_vcfhfindelpass_name)


if os.path.isfile(output_vcfhfsnpr_file) and os.path.isfile(output_vcfhfsnppass_file):
    print(YELLOW + DIM + output_vcfhfsnppass_file + " EXIST\nOmmiting PASS Filtering (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "PASS Filtering Variants (Recall-Group) in group " + group_name + END_FORMATTING)
    select_pass_variants(output_vcfhfsnpr_file, nocall_fr=0.2)
    select_pass_variants(output_vcfhfindelr_file, nocall_fr=0.2)


    ######################################################################
    ##############START RECALIBRATION AND FINAL CALL######################
    ######################################################################

print("\n\n" + BLUE + BOLD + "STARTING RECALIBATION IN GROUP: " + group_name + END_FORMATTING)

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    args.sample = sample
    args.output = os.path.abspath(args.output)

    if sample in sample_list_F:

        print("\n" + WHITE_BG + "RECALIBRATION AND CALL ON SAMPLE: " + sample + END_FORMATTING)

        ##############START BAM RECALIBRATION############
        #################################################

        ################BQSR AND APPLY BQSR##################
        #####################################################
        out_bqsr_name = sample + ".bqsr.bam"
        output_bqsr_file = os.path.join(out_map_dir, out_bqsr_name)

        if os.path.isfile(output_bqsr_file):
            print(YELLOW + DIM + output_bqsr_file + " EXIST\nOmmiting Recalibration for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Recalibration in sample " + sample + END_FORMATTING)
            recalibrate_bam(args, tb=True)

        #HAPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
        #######################################################
        out_gvcf_dir = os.path.join(args.output, "GVCF")
        out_gvcf_name = sample + ".g.vcf"
        output_gvcf_file = os.path.join(out_gvcf_dir, out_gvcf_name)

        #args.input_bam = output_bqsr_file

        if os.path.isfile(output_gvcf_file):
            print(YELLOW + DIM + output_gvcf_file + " EXIST\nOmmiting Haplotype Call for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Haplotype Calling in sample " + sample + END_FORMATTING)
            haplotype_caller(args, recalibrate=False, ploidy=args.ploidy, bamout=False, forceactive=False)
"""
        ###############################################################################################################################################
        #############################FOR COMPARING PURPOSE#############################################################################################
        ###############################################################################################################################################

        #CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
        #######################################################
        out_vcf_dir = os.path.join(args.output, "VCF")
        out_vcf_name = sample + ".raw.vcf"
        output_vcf_file = os.path.join(out_vcf_dir, out_vcf_name)

        if os.path.isfile(output_vcf_file):
            print(YELLOW + DIM + output_vcf_file + " EXIST\nOmmiting Variant Calling for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Variant Calling in sample " + sample + END_FORMATTING)
            call_variants(args, recalibrate=False, group=False)

        #SELECT VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
        #########################################################
        out_vcfsnp_name = sample + ".snp.vcf"
        output_vcfsnp_file = os.path.join(out_vcf_dir, out_vcfsnp_name)

        if os.path.isfile(output_vcfsnp_file):
            print(YELLOW + DIM + output_vcfsnp_file + " EXIST\nOmmiting Variant Selection for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Selecting Variants in sample " + sample + END_FORMATTING)
            select_variants(output_vcf_file, select_type='SNP') #select_variants(output_vcfr_file, select_type='INDEL')

        #HARD FILTER VARIANTS 2/2 FOR RECALIBRATION #############
        #########################################################
        out_vcfhfsnp_name = sample + ".snp.hf.vcf"
        output_vcfhfsnp_file = os.path.join(out_vcf_dir, out_vcfhfsnp_name)

        if os.path.isfile(output_vcfhfsnp_file):
            print(YELLOW + DIM + output_vcfhfsnp_file + " EXIST\nOmmiting Hard Filtering for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Hard Filtering Variants in sample " + sample + END_FORMATTING)
            hard_filter(output_vcfsnp_file, select_type='SNP')
"""
#ONCE ALL GVCF VARIANTS ARE CALLED, THEY ARE GATHERED AND FILTERED 
# FOR FINAL CALLING
######################################################################
##############START GROUP CALLING FOR FINAL CALL######################
######################################################################
group_name = args.output.split("/")[-1]
print("\n\n" + BLUE + BOLD + "STARTING JOINT CALL FOR FINAL CALL IN GROUP: " + group_name + END_FORMATTING + "\n")

#CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_gvcf_dir = os.path.join(args.output, "GVCF")
out_gvcf_name = group_name + ".cohort.g.vcf"
output_gvcf_file = os.path.join(out_gvcf_dir, out_gvcf_name)

if os.path.isfile(output_gvcf_file):
    print(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting GVCF Combination for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "GVCF Combination in group " + group_name + END_FORMATTING)
    combine_gvcf(args, recalibrate=False, all_gvcf=args.enrich_gvcf)

#CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_vcf_dir = os.path.join(args.output, "VCF")
out_vcf_name = group_name + ".cohort.raw.vcf"
output_vcf_file = os.path.join(out_vcf_dir, out_vcf_name)

if os.path.isfile(output_vcf_file):
    print(YELLOW + DIM + output_vcf_file + " EXIST\nOmmiting Variant Calling (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Variant Calling (Group) in group " + group_name + END_FORMATTING)
    call_variants(args, recalibrate=False, group=True)

#SELECT VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#########################################################
out_vcfsnp_name = group_name + ".cohort.snp.vcf"
out_vcfindel_name = group_name + ".cohort.indel.vcf"
output_vcfsnp_file = os.path.join(out_vcf_dir, out_vcfsnp_name)
output_vcfindel_file = os.path.join(out_vcf_dir, out_vcfindel_name)

if os.path.isfile(output_vcfsnp_file) and os.path.isfile(output_vcfindel_file):
    print(YELLOW + DIM + output_vcfsnp_file + " EXIST\nOmmiting Variant Selection (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Selecting Variants (Group) in group " + group_name + END_FORMATTING)
    select_variants(output_vcf_file, select_type='SNP')
    select_variants(output_vcf_file, select_type='INDEL')

#HARD FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnp_name = group_name + ".cohort.snp.hf.vcf"
out_vcfhfindel_name = group_name + ".cohort.indel.hf.vcf"
output_vcfhfsnp_file = os.path.join(out_vcf_dir, out_vcfhfsnp_name)
output_vcfhfindel_file = os.path.join(out_vcf_dir, out_vcfhfindel_name)


if os.path.isfile(output_vcfhfsnp_file) and os.path.isfile(output_vcfhfindel_file):
    print(YELLOW + DIM + output_vcfhfsnp_file + " EXIST\nOmmiting Hard Filtering (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Hard Filtering Variants (Group) in group " + group_name + END_FORMATTING)
    hard_filter(output_vcfsnp_file, select_type='SNP')
    hard_filter(output_vcfindel_file, select_type='INDEL')


#PASS FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnppass_name = group_name + ".cohort.snp.hf.pass.vcf"
out_vcfhfindelpass_name = group_name + ".cohort.indel.hf.pass.vcf"
output_vcfhfsnppass_file = os.path.join(out_vcf_dir, out_vcfhfsnppass_name)
output_vcfhfindelpass_file = os.path.join(out_vcf_dir, out_vcfhfindelpass_name)


if os.path.isfile(output_vcfhfindelpass_file) and os.path.isfile(output_vcfhfsnppass_file):
    print(YELLOW + DIM + output_vcfhfsnppass_file + " EXIST\nOmmiting PASS Filtering (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "PASS Filtering Variants (Group) in group " + group_name + END_FORMATTING)
    select_pass_variants(output_vcfhfsnp_file, nocall_fr=0.2)
    select_pass_variants(output_vcfhfindel_file, nocall_fr=0.2)


split_vcf_saples(output_vcfhfsnppass_file, sample_list=sample_list_F)



for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    args.sample = sample
    args.output = os.path.abspath(args.output)

    if sample in sample_list_F:

        print("\n" + WHITE_BG + "FINAL FILTERING IN SAMPLE " + sample + END_FORMATTING)

        ################FINAL VCF FILTERING##################
        #####################################################
        out_final_name = sample + ".snp.hf.pass.final.vcf"
        in_final_name = sample + ".snp.hf.pass.vcf"
        output_final_vcf = os.path.join(out_vcf_dir, out_final_name)
        in_final_vcf = os.path.join(out_vcf_dir, in_final_name)

        if os.path.isfile(output_final_vcf):
            print(YELLOW + DIM + output_final_vcf + " EXIST\nOmmiting Final filter for sample " + sample + END_FORMATTING)
        else:
            print(GREEN + "Final filter in sample " + sample + END_FORMATTING)
            vcf_consensus_filter(in_final_vcf,  distance=1, AF=0.75, QD=15, window_10=3)


print("\n\n" + MAGENTA + BOLD + "VARIANT CALL FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")

"""
#./snptb_runner.py -i /home/laura/ANALYSIS/Lofreq/coinfection_designed/raw -r reference/MTB_ancestorII_reference.fasta -o /home/laura/ANALYSIS/Lofreq/coinfection_designed/TEST -s sample_list.txt
#/home/laura/DEVELOP/SNPTB/snptb_runner.py -i /home/laura/RAW/Mixtas_Italia/ -r /home/laura/DATABASES/REFERENCES/ancestorII/MTB_ancestorII_reference.fasta -o /home/laura/ANALYSIS/Lofreq/coinfection_italy/
"""