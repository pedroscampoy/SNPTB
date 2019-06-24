#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np
import re
import subprocess
from misc import get_snpeff_path
from vcf_process import calculate_ALT_AD, obtain_output_dir

def picard_markdup(args, vcf_file, database="Mycobacterium_tuberculosis_h37rv VCF_FILE"):
    #http://snpeff.sourceforge.net/SnpEff_manual.html
    # java -jar snpEff.jar -ud 0 -c path/to/snpEff/snpEff.config -stats SAMPLE.starts.annot.html
    #   Mycobacterium_tuberculosis_h37rv VCF_FILE > FILE.annot
   
    vcf_file = os.path.abspath(vcf_file)

    sample = os.path.basename(vcf_file).split(".")[0]
    file_name = os.path.basename(vcf_file).split(".")[:-1]

    snpeff_config = get_snpeff_path()[0]
    snpeff_jar = get_snpeff_path()[1]

    annotate_output_dir = obtain_output_dir(args, "Annotate")
    #output_dir = os.path.abspath(args.output)
    stat_name = sample + ".annot.stat"
    annot_name = file_name + ".annot"

    stat_file = os.path.join(annotate_output_dir, stat_name)
    output_file = os.path.join(annotate_output_dir, annot_name)

    cmd = ["java", "-jar", snpeff_jar, "-ud", "0", "-c", snpeff_config, "-stats", stat_file, database, vcf_file]
    #execute_subprocess(cmd)
    with open(output_file, "w") as outfile:
        #calculate coverage and save it in th eoutput file
        subprocess.run(cmd,
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)


def import_annot_to_pandas(vcf_file, sep='\t'):
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            #print(next_line)
            next_line = f.readline()
    
    if first_line.endswith('VCFv4.2'):
        
        #Use first line as header
        dataframe = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)
        sample = dataframe.columns[-1]
        dataframe.rename(columns={sample:'sample'}, inplace=True)
        
        ann_head = ["Allele","Annotation","Annotation_Impact","Gene_Name",
                    "Gene_ID","Feature_Type","Feature_ID","Transcript_BioType",
                    "Rank","HGVS.c","HGVS.p","cDNA.pos / cDNA.length",
                    "CDS.pos / CDS.length","AA.pos / AA.length","Distance",
                    "ERRORS / WARNINGS / INFO"]
        
        for index, data_row in dataframe.iterrows():
            info_fields = re.findall(r';*([a-zA-Z]{1,20})=', data_row.INFO)
            info_values = re.findall(r'-?\d+\.?\d*e?[+-]?\d{0,2}', data_row.INFO)
            ann_values_re = re.search(r'ANN=(.*)', data_row.INFO)
            all_ann_values = ann_values_re.group(1)
            ann_values = all_ann_values.split(",")[0].split("|")[:16]
            
            format_fields = data_row['FORMAT'].split(":")
            format_values = data_row['sample'].split(":")
                                    
            for ifield, ivalue in zip(info_fields,info_values):
                dataframe.loc[index,ifield] = ivalue
                
            for ffield, fvalue in zip(format_fields,format_values):
                dataframe.loc[index,ffield] = fvalue
            
            dataframe.loc[index,'ANN'] = all_ann_values
            
            for ann_field,ann_value in zip(ann_head, ann_values):
                dataframe.loc[index,ann_field] = ann_value
            
        dataframe.rename(columns={'AF':'af'}, inplace=True)
        
        dataframe['len_AD'] = dataframe['AD'].str.split(",").str.len()
        dataframe['REF_AD'] = dataframe['AD'].str.split(",").str[0]
        #dataframe['ALT_AD'] = dataframe['AD'].str.split(",").str[1]
        dataframe['ALT_AD'] = dataframe.apply(calculate_ALT_AD, axis=1)
        dataframe[['gt0','gt1']] = dataframe['GT'].str.split(r'[/|\|]', expand=True)
                
        to_float = ['QUAL', 'AC', 'af', 'AN', 'BaseQRankSum', 'DP', 'ExcessHet', 'FS',
       'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR','GQ','ALT_AD', 'REF_AD']
        
        to_int = ['POS', 'len_AD', 'gt0', 'gt1']
        
        to_str = ['#CHROM','REF','ALT', 'FILTER']
        
        for column in dataframe.columns:
            if column in to_float:
                dataframe[column] = dataframe[column].astype(float)
                
        for column in dataframe.columns:
            if column in to_int:
                dataframe[column] = dataframe[column].astype(int)
                
        for column in dataframe.columns:
            if column in to_str:
                dataframe[column] = dataframe[column].astype(str)
                
        dataframe['dp'] = (dataframe['REF_AD'] + dataframe['ALT_AD'])
        dataframe['aF'] = dataframe['REF_AD']/dataframe['dp']
        dataframe['AF'] = dataframe['ALT_AD']/dataframe['dp']
        

                
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe