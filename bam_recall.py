# PYTHON_ARGCOMPLETE_OK

import os
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, obtain_output_dir, check_create_dir, get_picard_path, execute_subprocess, check_remove_file



"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 18 April 2019
REVISION:

TODO
    -Handle parameters
    -Explore samtools markdup

================================================================
END_OF_HEADER
================================================================
"""


#COLORS AND AND FORMATTING
"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
"""
END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
DIM = '\033[2m'

def get_arguments():

    parser = argparse.ArgumentParser(prog = 'bam_recalibrate.py', description= 'Call variants, get HQ SNPs and use them to realibrate')
    
    input_group = parser.add_argument_group('Input', 'Required files to map')

    input_group.add_argument('-b', '--bam', dest="input_bam", metavar="input.bam", type=str, required=True, help='Input bam to recalibrate and call')
    input_group.add_argument('-r', '--reference', metavar="ref.fasta", type=str, required=True, help='File to map against')

    output_group = parser.add_argument_group('Output')

    output_group.add_argument('-o', '--output', type=str, required=False, help='Output directory')
    output_group.add_argument('-s', '--sample', type=str, required=False, help='Sample name to handle output files ')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    
    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=1, help='Threads to use')


    #argcomplete.autocomplete(parser)
    arguments = parser.parse_args()

    return arguments

args = get_arguments()

def samtools_markdup(args):
    #http://www.htslib.org/doc/samtools.html
    # Add ms and MC tags for markdup to use later
    #samtools fixmate -m namesort.bam fixmate.bam
    # Markdup needs position order
    #samtools sort -o positionsort.bam fixmate.bam
    # Finally mark duplicates
    #samtools markdup positionsort.bam markdup.bam
    pass

def picard_markdup(args):
    #java -jar picard.jar MarkDuplicates \
    #  I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
    picard_jar = get_picard_path()
    
    input_bam = os.path.abspath(args.input_bam)
    in_param = "I=" + input_bam
    
    path_file_name = input_bam.split(".")[0]
    file_name = path_file_name.split("/")[-1]
    output_markdup = path_file_name + ".rg.markdup.bam"
    output_markdup_sorted = path_file_name + ".rg.markdup.sorted.bam"
    out_param = "O=" + output_markdup

    stat_output_dir = obtain_output_dir(args, "Stats")
    stat_output_file = file_name + ".markdup.metrics.txt"
    stat_output_full = os.path.join(stat_output_dir, stat_output_file)
    stats_param = "M=" + stat_output_full

    check_create_dir(stat_output_dir)

    cmd_markdup = ["java", "-jar", picard_jar, "MarkDuplicates", 
    in_param, out_param, stats_param]
    execute_subprocess(cmd_markdup)
    

    #samtools sort: samtools sort $output_dir/$sample".sorted.bam" -o $output_dir/$sample".sorted.bam"
    cmd_sort = ["samtools", "sort", output_markdup, "-o", output_markdup_sorted]
    execute_subprocess(cmd_sort)

    #Handled in Haplotype Caller function
    #samtools index: samtools index $output_dir/$sample".sorted.bam"
    #subprocess.run(["samtools", "index", output_markdup_sorted], 
    #    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    
    check_remove_file(output_markdup)

def picard_dictionary(args):
    #java -jar picard.jar CreateSequenceDictionary\
    # R=reference.fasta O=reference.dict
    picard_jar = get_picard_path()

    input_reference = os.path.abspath(args.reference)
    ref_param = "R=" + input_reference

    path_file_list = input_reference.split(".")[:-1]
    path_file_name = ".".join(path_file_list)
    dict_file_name = path_file_name + ".dict"
    out_param = "O=" + dict_file_name

    if os.path.exists(dict_file_name):
        print(dict_file_name + " already EXIST")
    else:
        cmd = ["java", "-jar", picard_jar, "CreateSequenceDictionary", 
        ref_param, out_param]
        execute_subprocess(cmd)


def samtools_faidx(args):
    #samtools faidx reference.fa

    input_reference = os.path.abspath(args.reference)
    fai_file_name = input_reference + ".fai"

    if os.path.exists(fai_file_name):
        print(fai_file_name + " already EXIST")
    else:
        cmd = ["samtools", "faidx", input_reference]
        execute_subprocess(cmd)


def haplotype_caller(args, recalibrate=False, ploidy=1, bamout=False, forceactive=False, intervals=False):
    """
    #No excuses
    https://software.broadinstitute.org/gatk/documentation/article?id=11081
    """
    input_bam = os.path.abspath(args.input_bam)
    input_reference = os.path.abspath(args.reference)
    
    path_file_name = input_bam.split(".")[0]
    file_name = path_file_name.split("/")[-1] #sample_name

    if recalibrate:
        output_markdup_sorted = path_file_name + ".rg.markdup.sorted.bam"

        gvcf_output_dir = obtain_output_dir(args, "GVCF_recal")
        gvcf_output_file = file_name + ".g.vcf"
        gvcf_output_full = os.path.join(gvcf_output_dir, gvcf_output_file)
    else:
        output_markdup_sorted = path_file_name + ".bqsr.bam"

        gvcf_output_dir = obtain_output_dir(args, "GVCF")
        gvcf_output_file = file_name + ".g.vcf"
        gvcf_output_full = os.path.join(gvcf_output_dir, gvcf_output_file)

    check_create_dir(gvcf_output_dir)

    hc_args = ["gatk", "HaplotypeCaller",
     "--reference", input_reference,
     "--input", output_markdup_sorted,
     "--output", gvcf_output_full,
     "--emit-ref-confidence", "GVCF",
     "--annotation-group", "AS_StandardAnnotation",
     "--sample-ploidy", str(ploidy)
     ]

    #Create bam index
    cmd_index = ["samtools", "index", output_markdup_sorted]
    execute_subprocess(cmd_index)

    if bamout:
        bamout_output_dir = obtain_output_dir(args, "Bamout")
        bamout_output_file = file_name + ".p" + str(ploidy) + ".out.bam"
        bamout_output_full = os.path.join(bamout_output_dir, bamout_output_file)
        check_create_dir(bamout_output_dir)
        bamout_params = ["--bam-output", bamout_output_full]
        hc_args.extend(bamout_params)
    
    if forceactive:
        force_params = ["--force-active", "--disable-optimizations"]
        hc_args.extend(force_params)

    execute_subprocess(hc_args)

    """
    gatk --java-options "-Xmx4g" HaplotypeCaller --reference ref.fasta --input markdup.bam --output sample.g.vcf" \
   --intervals MTB_anc_PAIR1_50_50:1*-100000 --emit-ref-confidence GVCF --bam-output "Bamout/"$i".out.bam" \
   --force-active --disable-optimizations --dont-trim-active-regions --max-num-haplotypes-in-population 2 \
   --assembly-region-padding 5000 --max-assembly-region-size 10000 --assembly-region-out "/home/laura/ANALYSIS/GATK/coinfection_PAIR15050/Assembly/"$i"_assembly.tsv" \
   --annotation-group AS_StandardAnnotation --annotation-group StandardHCAnnotation
    """

def call_variants(args, recalibrate=False):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
    #Call variants:
    gatk --java-options "-Xmx4g" GenotypeGVCFs -R Homo_sapiens_assembly38.fasta -V input.g.vcf.gz -O output.vcf.gz
    """
    input_bam = os.path.abspath(args.input_bam)
    input_reference = os.path.abspath(args.reference)
    
    path_file_name = input_bam.split(".")[0]
    file_name = path_file_name.split("/")[-1] #sample_name

    if recalibrate:
        gvcf_output_dir = obtain_output_dir(args, "GVCF_recal")
        gvcf_output_file = file_name + ".g.vcf"
        gvcf_output_full = os.path.join(gvcf_output_dir, gvcf_output_file)

        vcf_output_dir = obtain_output_dir(args, "VCF_recal")
        vcf_output_file = file_name + ".raw.vcf"
        vcf_output_full = os.path.join(vcf_output_dir, vcf_output_file)
    else:
        gvcf_output_dir = obtain_output_dir(args, "GVCF")
        gvcf_output_file = file_name + ".g.vcf"
        gvcf_output_full = os.path.join(gvcf_output_dir, gvcf_output_file)

        vcf_output_dir = obtain_output_dir(args, "VCF")
        vcf_output_file = file_name + ".raw.vcf"
        vcf_output_full = os.path.join(vcf_output_dir, vcf_output_file)

    check_create_dir(gvcf_output_dir)
    check_create_dir(vcf_output_dir)

    cmd = ["gatk", "GenotypeGVCFs", 
    "--reference", input_reference,
    "--variant", gvcf_output_full,
    "--output", vcf_output_full]

    execute_subprocess(cmd)
    

def select_variants(raw_vcf, select_type='SNP'):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php
    gatk SelectVariants -V cohort.vcf.gz -select-type SNP -O snps.vcf.gz
    """
    if select_type == "SNP":
        extension = ".snp.vcf"
    elif select_type == "INDEL":
        extension = ".indel.vcf"
    else:
        print(RED + BOLD + "Choose a correct type to filter" + END_FORMATTING)

    input_vcf = os.path.abspath(raw_vcf)
    check_file_exists(input_vcf)
    
    raw_vcf_file_name = input_vcf.split(".")[0]
    #file_name = raw_vcf_file_name.split("/")[-1] #sample_name

    vcf_selected_output_file = raw_vcf_file_name + extension

    cmd = ["gatk", "SelectVariants", 
    "--variant", input_vcf,
    "--select-type-to-include", select_type,
    "--output", vcf_selected_output_file]

    execute_subprocess(cmd)

def hard_filter(selected_vcf, select_type='SNP'):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.php
    https://software.broadinstitute.org/gatk/documentation/article?id=23216
    SNP:
    gatk VariantFiltration -V snps.vcf.gz "--filter-expression", "QD < 2.0", "--filter-name", "QD2" \
    "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30" "--filter-expression", "SOR > 3.0", "--filter-name", "SOR3" "--filter-expression", "FS > 60.0", "--filter-name", "FS60" \
    "--filter-expression", "MQ < 40.0", "--filter-name", "MQ40" "--filter-expression", "MQRankSum < -12.5", "--filter-name", "MQRankSum-12.5" "--filter-expression", "ReadPosRankSum < -8.0" \
   , "--filter-name", "ReadPosRankSum-8" -O snps_filtered.vcf.gz
    INDEL:
    gatk VariantFiltration -V indels.vcf.gz "--filter-expression", "QD < 2.0", "--filter-name", "QD2" "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30" \
    -"--filter-expression", "FS > 200.0", "--filter-name", "FS200" -"--filter-expression", "ReadPosRankSum < -20.0", "--filter-name", "ReadPosRankSum-20" -O indels_filtered.vcf.gz
    #--filterExpression "QD<2.0||FS>60.0||MQ<40.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filterName "my_snp_filter" 
    """

    input_vcf = os.path.abspath(selected_vcf)
    check_file_exists(input_vcf)
    
    selected_vcf_file_name = input_vcf.split(".")[0]

    if select_type == "SNP":
        extension = ".snp.hf.vcf"
        vcf_hard_filtered_output_file = selected_vcf_file_name + extension
        cmd = ["gatk", "VariantFiltration", 
            "--variant", input_vcf,
            "--filter-expression", "QD < 2.0", "--filter-name", "QD2",
            "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30",
            "--filter-expression", "SOR > 3.0", "--filter-name", "SOR3",
            "--filter-expression", "FS > 60.0", "--filter-name", "FS60",
            "--filter-expression", "MQ < 40.0", "--filter-name", "MQ40",
            "--filter-expression", "MQRankSum < -12.5", "--filter-name", "MQRankSum-12.5",
            "--filter-expression", "ReadPosRankSum < -8.0", "--filter-name", "ReadPosRankSum-8",
            "--output", vcf_hard_filtered_output_file]
    elif select_type == "INDEL":
        extension = ".indel.hf.vcf"
        vcf_hard_filtered_output_file = selected_vcf_file_name + extension
        cmd = ["gatk", "VariantFiltration", 
            "--variant", input_vcf,
            "--filter-expression", "QD < 2.0", "--filter-name", "QD2",
            "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30",
            "--filter-expression", "SOR > 10.0", "--filter-name", "SOR10",
            "--filter-expression", "FS > 200.0", "--filter-name", "FS200",
            "--filter-expression", "ReadPosRankSum < -20.0", "--filter-name", "ReadPosRankSum-20",
            "--output", vcf_hard_filtered_output_file]
    else:
        print(RED + BOLD + "Choose a correct type to filter" + END_FORMATTING)

    execute_subprocess(cmd)
    
 



#BaseRecalibrator

#     gatk BaseRecalibrator \
#    --input my_reads.bam \
#    --reference reference.fasta \
#    --known-sites sites_of_variation.vcf \
#    --known-sites another/optional/setOfSitesToMask.vcf \
#    --output recal_data.table

    #ApplyBQSR

# gatk ApplyBQSR \
#    -R reference.fasta \
#    -I input.bam \
#    --bqsr-recal-file recalibration.table \
#    -O output.bam



#print(args)

#picard_markdup(args)
#picard_dictionary(args)
#samtools_faidx(args)
#haplotype_caller(args, recalibrate=True, ploidy=1, bamout=False, forceactive=False)
#call_variants(args, recalibrate=True)

out = args.output
sample = args.sample
#raw_vcf = out + "/VCF_recal/" + sample + ".raw.vcf"
#select_variants(raw_vcf, select_type='SNP') #select_variants(raw_vcf, select_type='INDEL')

selected_vcf = out + "/VCF_recal/" + sample + ".snp.vcf"
hard_filter(selected_vcf, select_type='SNP')


#Filter highest quality
#Use it to recalibrate (BaseRecalibrator and ApplyBQSR in one function)

#python bam_recall.py -b  Bam/AL14621.rg.sorted.bam -r ../references/NC_000962.3.fasta -o .
#python bam_recall.py -b /home/pedro/analysis/Mixed/Bam/1mixed.rg.sorted.bam -r ../references/NC_000962.3.fasta -o /home/pedro/analysis/Mixed