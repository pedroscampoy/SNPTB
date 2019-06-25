#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np
import re
import subprocess
from misc import get_snpeff_path, check_create_dir
from vcf_process import calculate_ALT_AD, obtain_output_dir


def replace_reference(input, ref_old, ref_new, output):
    """
    THis function replace all instances of a reference in a vcf file
    """
    input_file = os.path.abspath(input)
    output_file = os.path.abspath(output)
    output_dir = os.path.dirname(output)

    check_create_dir(output_dir)

    with open(input_file, 'r') as fi:
        with open(output_file, 'w') as fo:
            for line in fi:
                ref = ref_old + "\t"
                new = ref_new + "\t"
                line = line.replace(ref, new)
                fo.write(line)



def snpeff_annotation(args, vcf_file, database="Mycobacterium_tuberculosis_h37rv"):
    #http://snpeff.sourceforge.net/SnpEff_manual.html
    # java -jar snpEff.jar -ud 0 -c path/to/snpEff/snpEff.config -stats SAMPLE.starts.annot.html
    #   Mycobacterium_tuberculosis_h37rv VCF_FILE > FILE.annot
   
    vcf_file = os.path.abspath(vcf_file)

    sample = os.path.basename(vcf_file).split(".")[0]
    file_name = (".").join(os.path.basename(vcf_file).split(".")[:-1])

    snpeff_config = get_snpeff_path()[0]
    snpeff_jar = get_snpeff_path()[1]

    annotate_output_dir = obtain_output_dir(args, "Annotation")
    #output_dir = os.path.abspath(args.output)
    stat_name = sample + ".annot.html"
    annot_name = file_name + ".annot"

    stat_file = os.path.join(annotate_output_dir, stat_name)
    output_file = os.path.join(annotate_output_dir, annot_name)

    cmd = ["java", "-jar", snpeff_jar, "-ud", "0", "-c", snpeff_config, "-stats", stat_file, database, vcf_file]
    #execute_subprocess(cmd)
    with open(output_file, "w+") as outfile:
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
            ann_values_re = re.search(r'ANN=(.*)\|(.*)', data_row.INFO)
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



def add_lineage_Coll(vcf_df):
    dict_lineage_position = {
        '615938' : ['A', '1'],
        '4404247' : ['A', '1.1'],
        '3021283' : ['A', '1.1.1'],
        '3216553' : ['A', '1.1.1.1'],
        '2622402' : ['A', '1.1.2'],
        '1491275' : ['A', '1.1.3'],
        '3479545' : ['A', '1.2.1'],
        '3470377' : ['T', '1.2.2'],
        '497491' : ['A', '2'],
        '1881090' : ['T', '2.1'],
        '2505085' : ['A', '2.2'],
        '797736' : ['T', '2.2.1'],
        '4248115' : ['T', '2.2.1.1'],
        '3836274' : ['A', '2.2.1.2'],
        '346693' : ['T', '2.2.2'],
        '3273107' : ['A', '3'],
        '1084911' : ['A', '3.1.1'],
        '3722702' : ['C', '3.1.2'],
        '1237818' : ['G', '3.1.2.1'],
        '2874344' : ['A', '3.1.2.2'],
        '931123' : ['C', '4**'],
        '62657' : ['A', '4.1'],
        '514245' : ['T', '4.1.1'],
        '1850119' : ['T', '4.1.1.1'],
        '541048' : ['G', '4.1.1.2'],
        '4229087' : ['T', '4.1.1.3'],
        '891756' : ['G', '4.1.2'],
        '107794' : ['T', '4.1.2.1'],
        '2411730' : ['C', '4.2'],
        '783601' : ['C', '4.2.1'],
        '1487796' : ['A', '4.2.2'],
        '1455780' : ['C', '4.2.2.1'],
        '764995' : ['G', '4.3'],
        '615614' : ['A', '4.3.1'],
        '4316114' : ['A', '4.3.2'],
        '3388166' : ['G', '4.3.2.1'],
        '403364' : ['A', '4.3.3'],
        '3977226' : ['A', '4.3.4'],
        '4398141' : ['A', '4.3.4.1'],
        '1132368' : ['T', '4.3.4.2'],
        '1502120' : ['A', '4.3.4.2.1'],
        '4307886' : ['A', '4.4'],
        '4151558' : ['A', '4.4.1'],
        '355181' : ['A', '4.4.1.1'],
        '2694560' : ['C', '4.4.1.2'],
        '4246508' : ['A', '4.4.2'],
        '1719757' : ['T', '4.5'],
        '3466426' : ['A', '4.6'],
        '4260268' : ['C', '4.6.1'],
        '874787' : ['A', '4.6.1.1'],
        '1501468' : ['C', '4.6.1.2'],
        '4125058' : ['C', '4.6.2'],
        '3570528' : ['G', '4.6.2.1'],
        '2875883' : ['T', '4.6.2.2'],
        '4249732' : ['G', '4.7'],
        '3836739' : ['A', '4.8'],
        '1759252' : ['T', '4.9**'],
        '1799921' : ['A', '5'],
        '1816587' : ['G', '6'],
        '1137518' : ['A', '7'],
        '2831482' : ['G', 'BOV'],
        '1882180' : ['T', 'BOV_AFRI']
                }
    list_lineage = []
    
    for index, _ in vcf_df.iterrows():
        position = str(vcf_df.loc[index,'POS'])
        if position in dict_lineage_position.keys():
            if str(vcf_df.loc[index,'ALT']) == dict_lineage_position[str(position)][0]:
                lineage = dict_lineage_position[str(position)][1]
                vcf_df.loc[index,'Lineage'] = lineage
                list_lineage.append(lineage)
                
    if len(list_lineage) > 0:
        list_lineage.sort(reverse=True)
        asterix = ""
        for sublineage_n in range(len(list_lineage)):
            if sublineage_n < (len(list_lineage) - 1):
                if list_lineage[sublineage_n].startswith(list_lineage[sublineage_n + 1]):
                    asterix = asterix + "*"
        final_lineage = list_lineage[0] + " " + asterix
        print("This strain has lineage position(s):\n: " + " ".join([list_lineage[0],asterix]))
        return final_lineage
    else:
        print("No lineage were found\n")


def add_resistance_snp(vcf_df):
    dict_resistance_position = {6575: ['T', 'fluoroquinolones (FQ)'],
                                6620: ['C', 'fluoroquinolones (FQ)', 'A'],
                                6621: ['C', 'fluoroquinolones (FQ)'],
                                6734: ['G', 'fluoroquinolones (FQ)'],
                                6735: ['C', 'fluoroquinolones (FQ)'],
                                6736: ['G', 'fluoroquinolones (FQ)'],
                                6737: ['C', 'fluoroquinolones (FQ)'],
                                6738: ['A', 'fluoroquinolones (FQ)'],
                                6741: ['T', 'fluoroquinolones (FQ)'],
                                6742: ['T', 'fluoroquinolones (FQ)'],
                                6749: ['A', 'fluoroquinolones (FQ)'],
                                6750: ['T', 'fluoroquinolones (FQ)'],
                                7563: ['T', 'fluoroquinolones (FQ)'],
                                7564: ['C', 'fluoroquinolones (FQ)'],
                                7566: ['A', 'fluoroquinolones (FQ)'],
                                7570: ['T', 'fluoroquinolones (FQ)'],
                                7572: ['C', 'fluoroquinolones (FQ)'],
                                7581: ['C', 'fluoroquinolones (FQ)', 'A', 'T'],
                                7582: ['G', 'fluoroquinolones (FQ)', 'C', 'T'],
                                575729: ['T', 'ethionamide (ETH)'],
                                576164: ['T', 'ethionamide (ETH)'],
                                576242: ['T', 'ethionamide (ETH)'],
                                576338: ['T', 'ethionamide (ETH)'],
                                576414: ['A', 'ethionamide (ETH)'],
                                576429: ['C', 'ethionamide (ETH)'],
                                760314: ['T', 'rifampicin (RMP)'],
                                761004: ['G', 'rifampicin (RMP)'],
                                761093: ['C', 'rifampicin (RMP)'],
                                761095: ['C', 'rifampicin (RMP)', 'G'],
                                761098: ['T', 'rifampicin (RMP)', 'C'],
                                761100: ['A', 'rifampicin (RMP)'],
                                761101: ['C', 'rifampicin (RMP)', 'T'],
                                761108: ['T', 'rifampicin (RMP)'],
                                761109: ['T', 'rifampicin (RMP)'],
                                761110: ['G', 'rifampicin (RMP)', 'T'],
                                761111: ['G', 'rifampicin (RMP)'],
                                761120: ['C', 'rifampicin (RMP)', 'G'],
                                761128: ['T', 'rifampicin (RMP)', 'G'],
                                761139: ['A', 'rifampicin (RMP)', 'G', 'T'],
                                761140: ['C', 'rifampicin (RMP)', 'G', 'T'],
                                761141: ['A', 'rifampicin (RMP)'],
                                761154: ['G', 'rifampicin (RMP)'],
                                761155: ['G', 'rifampicin (RMP)', 'T'],
                                761161: ['C', 'rifampicin (RMP)'],
                                761277: ['T', 'rifampicin (RMP)'],
                                781687: ['G', 'streptomycin (SM)'],
                                781821: ['C', 'streptomycin (SM)'],
                                781822: ['G', 'streptomycin (SM)'],
                                801268: ['C', 'linezolid (LZD)'],
                                1472337: ['T', 'streptomycin (SM)'],
                                1472358: ['T', 'streptomycin (SM)'],
                                1472359: ['C', 'streptomycin (SM)'],
                                1472362: ['T', 'streptomycin (SM)'],
                                1472750: ['A', 'streptomycin (SM)'],
                                1472751: ['G', 'streptomycin (SM)'],
                                1472752: ['T', 'streptomycin (SM)'],
                                1473246: ['G', 'amikacin (AMK) kanamycin (KAN) capreomycin (CPR)'],
                                1473247: ['T', 'amikacin (AMK) kanamycin (KAN) capreomycin (CPR)'],
                                1473329: ['T', 'amikacin (AMK) kanamycin (KAN) capreomycin (CPR)'],
                                1475956: ['T', 'linezolid (LZD)'],
                                1476471: ['T', 'linezolid (LZD)'],
                                1673423: ['T', 'isoniazid (INH)'],
                                1673424: ['G', 'isoniazid (INH)'],
                                1673425: ['T', 'isoniazid (INH)'],
                                1673432: ['A', 'isoniazid (INH)', 'C'],
                                1674481: ['G', 'isoniazid (INH) ethionamide (ETH)'],
                                1674782: ['C', 'isoniazid (INH) ethionamide (ETH)'],
                                1833909: ['C', 'pyrazinamide (PZA)'],
                                1834325: ['A', 'pyrazinamide (PZA)'],
                                1834855: ['C', 'pyrazinamide (PZA)'],
                                1917946: ['T', 'capreomycin (CPR)'],
                                1917979: ['T', 'capreomycin (CPR)'],
                                1917991: ['T', 'capreomycin (CPR)'],
                                1918003: ['T', 'capreomycin (CPR)'],
                                1918139: ['A', 'capreomycin (CPR)'],
                                1918144: ['G', 'capreomycin (CPR)'],
                                1918211: ['A', 'capreomycin (CPR)'],
                                1918292: ['C', 'capreomycin (CPR)'],
                                1918322: ['A', 'capreomycin (CPR)'],
                                1918388: ['C', 'capreomycin (CPR)'],
                                1918487: ['T', 'capreomycin (CPR)'],
                                1918489: ['A', 'capreomycin (CPR)', 'T'],
                                1918494: ['G', 'capreomycin (CPR)'],
                                1918651: ['A', 'capreomycin (CPR)'],
                                2102240: ['T', 'isoniazid (INH) ethionamide (ETH)'],
                                2102715: ['C', 'isoniazid (INH) ethionamide (ETH)'],
                                2155167: ['C', 'isoniazid (INH)', 'T'],
                                2155168: ['G', 'isoniazid (INH)', 'T', 'A'],
                                2155169: ['C', 'isoniazid (INH)'],
                                2155206: ['C', 'isoniazid (INH)'],
                                2155212: ['G', 'isoniazid (INH)'],
                                2155214: ['C', 'isoniazid (INH)'],
                                2155222: ['A', 'isoniazid (INH)'],
                                2155289: ['G', 'isoniazid (INH)'],
                                2155699: ['C', 'isoniazid (INH)'],
                                2288683: ['G', 'pyrazinamide (PZA)'],
                                2288697: ['G', 'pyrazinamide (PZA)'],
                                2288703: ['G', 'pyrazinamide (PZA)', 'C'],
                                2288704: ['A', 'pyrazinamide (PZA)'],
                                2288718: ['G', 'pyrazinamide (PZA)'],
                                2288719: ['C', 'pyrazinamide (PZA)'],
                                2288727: ['G', 'pyrazinamide (PZA)'],
                                2288730: ['A', 'pyrazinamide (PZA)'],
                                2288740: ['G', 'pyrazinamide (PZA)'],
                                2288752: ['G', 'pyrazinamide (PZA)'],
                                2288754: ['G', 'pyrazinamide (PZA)'],
                                2288757: ['T', 'pyrazinamide (PZA)'],
                                2288761: ['G', 'pyrazinamide (PZA)'],
                                2288764: ['G', 'pyrazinamide (PZA)'],
                                2288766: ['C', 'pyrazinamide (PZA)'],
                                2288772: ['G', 'pyrazinamide (PZA)', 'C'],
                                2288778: ['C', 'pyrazinamide (PZA)'],
                                2288779: ['T', 'pyrazinamide (PZA)'],
                                2288782: ['C', 'pyrazinamide (PZA)'],
                                2288805: ['T', 'pyrazinamide (PZA)', 'A'],
                                2288806: ['G', 'pyrazinamide (PZA)', 'T'],
                                2288817: ['T', 'pyrazinamide (PZA)', 'A'],
                                2288818: ['C', 'pyrazinamide (PZA)'],
                                2288820: ['G', 'pyrazinamide (PZA)'],
                                2288821: ['A', 'pyrazinamide (PZA)'],
                                2288823: ['G', 'pyrazinamide (PZA)'],
                                2288826: ['C', 'pyrazinamide (PZA)'],
                                2288827: ['G', 'pyrazinamide (PZA)'],
                                2288828: ['C', 'pyrazinamide (PZA)'],
                                2288830: ['G', 'pyrazinamide (PZA)'],
                                2288832: ['C', 'pyrazinamide (PZA)', 'G'],
                                2288833: ['C', 'pyrazinamide (PZA)'],
                                2288836: ['T', 'pyrazinamide (PZA)', 'A'],
                                2288838: ['T', 'pyrazinamide (PZA)'],
                                2288839: ['G', 'pyrazinamide (PZA)'],
                                2288841: ['A', 'pyrazinamide (PZA)'],
                                2288844: ['G', 'pyrazinamide (PZA)'],
                                2288847: ['G', 'pyrazinamide (PZA)', 'T'],
                                2288848: ['A', 'pyrazinamide (PZA)', 'T'],
                                2288853: ['G', 'pyrazinamide (PZA)', 'C'],
                                2288857: ['A', 'pyrazinamide (PZA)'],
                                2288859: ['C', 'pyrazinamide (PZA)'],
                                2288868: ['C', 'pyrazinamide (PZA)'],
                                2288869: ['A', 'pyrazinamide (PZA)'],
                                2288874: ['G', 'pyrazinamide (PZA)'],
                                2288880: ['G', 'pyrazinamide (PZA)'],
                                2288883: ['C', 'pyrazinamide (PZA)', 'T', 'G'],
                                2288885: ['T', 'pyrazinamide (PZA)'],
                                2288886: ['G', 'pyrazinamide (PZA)', 'T'],
                                2288887: ['G', 'pyrazinamide (PZA)', 'C'],
                                2288895: ['C', 'pyrazinamide (PZA)', 'G'],
                                2288902: ['G', 'pyrazinamide (PZA)'],
                                2288920: ['T', 'pyrazinamide (PZA)', 'G'],
                                2288928: ['T', 'pyrazinamide (PZA)'],
                                2288930: ['T', 'pyrazinamide (PZA)'],
                                2288931: ['A', 'pyrazinamide (PZA)'],
                                2288933: ['C', 'pyrazinamide (PZA)'],
                                2288934: ['C', 'pyrazinamide (PZA)', 'G'],
                                2288935: ['C', 'pyrazinamide (PZA)'],
                                2288938: ['T', 'pyrazinamide (PZA)'],
                                2288944: ['C', 'pyrazinamide (PZA)', 'G'],
                                2288945: ['T', 'pyrazinamide (PZA)'],
                                2288952: ['G', 'pyrazinamide (PZA)', 'T'],
                                2288953: ['T', 'pyrazinamide (PZA)'],
                                2288954: ['T', 'pyrazinamide (PZA)'],
                                2288955: ['C', 'pyrazinamide (PZA)', 'G'],
                                2288956: ['G', 'pyrazinamide (PZA)', 'C'],
                                2288957: ['C', 'pyrazinamide (PZA)'],
                                2288960: ['T', 'pyrazinamide (PZA)', 'C'],
                                2288961: ['C', 'pyrazinamide (PZA)', 'G'],
                                2288962: ['G', 'pyrazinamide (PZA)'],
                                2288971: ['A', 'pyrazinamide (PZA)'],
                                2288982: ['A', 'pyrazinamide (PZA)'],
                                2288988: ['C', 'pyrazinamide (PZA)', 'G'],
                                2288997: ['C', 'pyrazinamide (PZA)'],
                                2288998: ['C', 'pyrazinamide (PZA)'],
                                2289000: ['C', 'pyrazinamide (PZA)', 'G'],
                                2289001: ['C', 'pyrazinamide (PZA)'],
                                2289016: ['G', 'pyrazinamide (PZA)'],
                                2289028: ['G', 'pyrazinamide (PZA)'],
                                2289029: ['T', 'pyrazinamide (PZA)'],
                                2289030: ['C', 'pyrazinamide (PZA)'],
                                2289031: ['A', 'pyrazinamide (PZA)'],
                                2289036: ['A', 'pyrazinamide (PZA)'],
                                2289038: ['G', 'pyrazinamide (PZA)', 'A'],
                                2289039: ['G', 'pyrazinamide (PZA)', 'T'],
                                2289040: ['G', 'pyrazinamide (PZA)', 'C'],
                                2289043: ['G', 'pyrazinamide (PZA)'],
                                2289050: ['C', 'pyrazinamide (PZA)'],
                                2289052: ['C', 'pyrazinamide (PZA)'],
                                2289054: ['C', 'pyrazinamide (PZA)'],
                                2289057: ['A', 'pyrazinamide (PZA)'],
                                2289061: ['G', 'pyrazinamide (PZA)'],
                                2289068: ['T', 'pyrazinamide (PZA)', 'C'],
                                2289070: ['G', 'pyrazinamide (PZA)'],
                                2289071: ['C', 'pyrazinamide (PZA)'],
                                2289072: ['C', 'pyrazinamide (PZA)', 'G'],
                                2289073: ['C', 'pyrazinamide (PZA)', 'A'],
                                2289081: ['C', 'pyrazinamide (PZA)', 'T', 'A'],
                                2289082: ['A', 'pyrazinamide (PZA)'],
                                2289089: ['T', 'pyrazinamide (PZA)'],
                                2289090: ['C', 'pyrazinamide (PZA)', 'G'],
                                2289091: ['A', 'pyrazinamide (PZA)'],
                                2289096: ['G', 'pyrazinamide (PZA)', 'C'],
                                2289097: ['T', 'pyrazinamide (PZA)'],
                                2289100: ['C', 'pyrazinamide (PZA)', 'A'],
                                2289103: ['C', 'pyrazinamide (PZA)', 'G'],
                                2289108: ['C', 'pyrazinamide (PZA)'],
                                2289111: ['C', 'pyrazinamide (PZA)'],
                                2289133: ['A', 'pyrazinamide (PZA)'],
                                2289138: ['G', 'pyrazinamide (PZA)'],
                                2289140: ['C', 'pyrazinamide (PZA)'],
                                2289150: ['C', 'pyrazinamide (PZA)'],
                                2289159: ['T', 'pyrazinamide (PZA)'],
                                2289162: ['G', 'pyrazinamide (PZA)'],
                                2289171: ['T', 'pyrazinamide (PZA)'],
                                2289180: ['C', 'pyrazinamide (PZA)'],
                                2289186: ['G', 'pyrazinamide (PZA)'],
                                2289193: ['T', 'pyrazinamide (PZA)'],
                                2289200: ['T', 'pyrazinamide (PZA)'],
                                2289201: ['T', 'pyrazinamide (PZA)'],
                                2289202: ['G', 'pyrazinamide (PZA)'],
                                2289203: ['C', 'pyrazinamide (PZA)'],
                                2289204: ['G', 'pyrazinamide (PZA)'],
                                2289206: ['C', 'pyrazinamide (PZA)'],
                                2289207: ['G', 'pyrazinamide (PZA)'],
                                2289208: ['T', 'pyrazinamide (PZA)'],
                                2289213: ['C', 'pyrazinamide (PZA)', 'G'],
                                2289214: ['T', 'pyrazinamide (PZA)'],
                                2289216: ['G', 'pyrazinamide (PZA)', 'C'],
                                2289218: ['T', 'pyrazinamide (PZA)'],
                                2289219: ['G', 'pyrazinamide (PZA)', 'C'],
                                2289220: ['T', 'pyrazinamide (PZA)'],
                                2289222: ['T', 'pyrazinamide (PZA)', 'C'],
                                2289223: ['A', 'pyrazinamide (PZA)'],
                                2289225: ['G', 'pyrazinamide (PZA)'],
                                2289231: ['G', 'pyrazinamide (PZA)', 'C'],
                                2289234: ['T', 'pyrazinamide (PZA)'],
                                2289235: ['G', 'pyrazinamide (PZA)'],
                                2289239: ['A', 'pyrazinamide (PZA)', 'T'],
                                2289240: ['T', 'pyrazinamide (PZA)', 'G'],
                                2289248: ['G', 'pyrazinamide (PZA)', 'C'],
                                2289252: ['C', 'pyrazinamide (PZA)', 'G', 'A'],
                                2715342: ['T', 'kanamycin (KAN)'],
                                2715346: ['A', 'kanamycin (KAN)'],
                                2726136: ['T', 'isoniazid (INH)'],
                                2726145: ['A', 'isoniazid (INH)'],
                                3073808: ['C', 'para-aminosalicylic acid (PAS)'],
                                4241078: ['G', 'ethambutol (EMB)'],
                                4243221: ['T', 'ethambutol (EMB)'],
                                4243225: ['A', 'ethambutol (EMB)'],
                                4243242: ['A', 'ethambutol (EMB)'],
                                4243245: ['A', 'ethambutol (EMB)'],
                                4243833: ['A', 'ethambutol (EMB)'],
                                4244193: ['A', 'ethambutol (EMB)'],
                                4244281: ['A', 'ethambutol (EMB)'],
                                4244617: ['T', 'ethambutol (EMB)'],
                                4245730: ['C', 'ethambutol (EMB)'],
                                4246734: ['G', 'ethambutol (EMB)'],
                                4247402: ['G', 'ethambutol (EMB)'],
                                4247429: ['G', 'ethambutol (EMB)', 'C'],
                                4247430: ['C', 'ethambutol (EMB)'],
                                4247431: ['A', 'ethambutol (EMB)', 'C', 'T'],
                                4247469: ['C', 'ethambutol (EMB)'],
                                4247495: ['T', 'ethambutol (EMB)'],
                                4247496: ['G', 'ethambutol (EMB)'],
                                4247507: ['C', 'ethambutol (EMB)'],
                                4247513: ['C', 'ethambutol (EMB)'],
                                4247573: ['A', 'ethambutol (EMB)'],
                                4247717: ['G', 'ethambutol (EMB)'],
                                4247723: ['T', 'ethambutol (EMB)'],
                                4247729: ['A', 'ethambutol (EMB)', 'T'],
                                4247730: ['C', 'ethambutol (EMB)', 'A'],
                                4247863: ['G', 'ethambutol (EMB)'],
                                4247873: ['A', 'ethambutol (EMB)'],
                                4248002: ['A', 'ethambutol (EMB)'],
                                4248003: ['G', 'ethambutol (EMB)'],
                                4248747: ['A', 'ethambutol (EMB)'],
                                4249518: ['G', 'ethambutol (EMB)'],
                                4326087: ['T', 'ethionamide (ETH)'],
                                4326236: ['T', 'ethionamide (ETH)'],
                                4326300: ['C', 'ethionamide (ETH)'],
                                4326320: ['T', 'ethionamide (ETH)'],
                                4326333: ['G', 'ethionamide (ETH)'],
                                4326449: ['T', 'ethionamide (ETH)'],
                                4326461: ['C', 'ethionamide (ETH)'],
                                4326738: ['A', 'ethionamide (ETH)'],
                                4326807: ['T', 'ethionamide (ETH)'],
                                4326917: ['T', 'ethionamide (ETH)'],
                                4327224: ['C', 'ethionamide (ETH)'],
                                4327301: ['G', 'ethionamide (ETH)'],
                                4327307: ['G', 'ethionamide (ETH)'],
                                4327322: ['A', 'ethionamide (ETH)'],
                                4327346: ['T', 'ethionamide (ETH)'],
                                4327347: ['A', 'ethionamide (ETH)'],
                                4407604: ['T', 'streptomycin (SM)'],
                                4407790: ['A', 'streptomycin (SM)'],
                                4407824: ['A', 'streptomycin (SM)'],
                                4407931: ['G', 'streptomycin (SM)'],
                                4407940: ['G', 'streptomycin (SM)'],
                                4407992: ['T', 'streptomycin (SM)'],
                                4408009: ['C', 'streptomycin (SM)'],
                                4408102: ['G', 'streptomycin (SM)'],
                                6576: ['A', 'fluoroquinolones_(FQ)'],
                                6579: ['T', 'fluoroquinolones_(FQ)', 'A'],
                                6768: ['C', 'fluoroquinolones_(FQ)'],
                                412339: ['G', 'ethambutol_(EMB)'],
                                413498: ['G', 'ethambutol_(EMB)'],
                                413807: ['T', 'ethambutol_(EMB)'],
                                761106: ['C', 'rifampicin_(RMP)'],
                                761116: ['G', 'rifampicin_(RMP)'],
                                761127: ['CA', 'rifampicin_(RMP)'],
                                1416212: ['C', 'ethambutol_(EMB)'],
                                1417019: ['T', 'ethambutol_(EMB)'],
                                1673393: ['C', 'isoniazid_(INH)'],
                                1673406: ['T', 'isoniazid_(INH)'],
                                1673431: ['A', 'isoniazid_(INH)'],
                                1918125: ['A', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)'],
                                1918135: ['C', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)', 'C'],
                                1918136: ['A', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)'],
                                1918202: ['G', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)'],
                                1918213: ['C', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)'],
                                1918219: ['T', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)', 'A'],
                                1918220: ['A', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)'],
                                1918418: ['G', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)'],
                                1918478: ['G', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)', 'C'],
                                1918661: ['C', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)'],
                                2154445: ['T', 'isoniazid_(INH)'],
                                2154682: ['T', 'isoniazid_(INH)'],
                                2154700: ['C', 'isoniazid_(INH)'],
                                2155211: ['A', 'isoniazid_(INH)'],
                                2155234: ['G', 'isoniazid_(INH)'],
                                2155241: ['T', 'isoniazid_(INH)'],
                                2155306: ['T', 'isoniazid_(INH)'],
                                2155688: ['T', 'isoniazid_(INH)'],
                                2288686: ['A', 'pyrazinamide_(PZA)'],
                                2288717: ['T', 'pyrazinamide_(PZA)'],
                                2288767: ['C', 'pyrazinamide_(PZA)'],
                                2288790: ['G', 'pyrazinamide_(PZA)'],
                                2288815: ['T', 'pyrazinamide_(PZA)'],
                                2288860: ['A', 'pyrazinamide_(PZA)'],
                                2289099: ['C', 'pyrazinamide_(PZA)'],
                                2726100: ['T', 'isoniazid_(INH)'],
                                2726112: ['T', 'isoniazid_(INH)'],
                                3645524: ['T', 'ethambutol_(EMB)'],
                                3646959: ['T', 'ethambutol_(EMB)'],
                                3647041: ['G', 'ethambutol_(EMB)'],
                                4245147: ['T', 'ethambutol_(EMB)'],
                                4245969: ['T', 'ethambutol_(EMB)'],
                                4407904: ['A', 'streptomycin_(SM)_amikacin_(AMK)_kanamycin_(KAN)_capreomycin_(CPR)']
                                }
    dict_high_confidence = {6620: 'A', 6621: 'C', 6734: 'G', 6735: 'C', 6736: 'G', 6737: 'C', 6738: 'A', 6741: 'T', 6742: 'T',
                         6750: 'T', 7563: 'T', 7564: 'C', 7570: 'T', 7572: 'C', 7581: 'T', 7582: 'C', 760314: 'T', 761101: 'T', 
                         761109: 'T', 761110: 'T', 761139: 'T', 761140: 'G', 761155: 'T', 761161: 'C', 761277: 'T', 781687: 'G', 
                         781822: 'G', 1473246: 'G', 1473247: 'T', 1473329: 'T', 1673425: 'T', 1673432: 'C', 1674481: 'G', 2155168: 'A', 
                         2155169: 'C', 2155214: 'C', 2155289: 'G', 2288683: 'G', 2288697: 'G', 2288703: 'C', 2288718: 'G', 2288719: 'C', 
                         2288740: 'G', 2288752: 'G', 2288754: 'G', 2288761: 'G', 2288764: 'G', 2288772: 'C', 2288778: 'C', 2288779: 'T', 
                         2288805: 'A', 2288806: 'G', 2288817: 'A', 2288818: 'C', 2288823: 'G', 2288826: 'C', 2288827: 'G', 2288828: 'C', 
                         2288830: 'G', 2288832: 'G', 2288833: 'C', 2288838: 'T', 2288839: 'G', 2288841: 'A', 2288847: 'T', 2288848: 'A', 
                         2288853: 'C', 2288857: 'A', 2288868: 'C', 2288869: 'A', 2288874: 'G', 2288880: 'G', 2288883: 'G', 2288886: 'T', 
                         2288887: 'C', 2288895: 'G', 2288920: 'G', 2288928: 'T', 2288930: 'T', 2288931: 'A', 2288933: 'C', 2288934: 'C', 
                         2288935: 'C', 2288938: 'T', 2288944: 'G', 2288945: 'T', 2288952: 'T', 2288954: 'T', 2288955: 'G', 2288956: 'G', 
                         2288957: 'C', 2288960: 'C', 2288961: 'G', 2288962: 'G', 2288971: 'A', 2288988: 'C', 2288997: 'C', 2288998: 'C', 
                         2289000: 'G', 2289001: 'C', 2289028: 'G', 2289029: 'T', 2289030: 'C', 2289038: 'A', 2289039: 'T', 2289043: 'G', 
                         2289050: 'C', 2289057: 'A', 2289068: 'C', 2289070: 'G', 2289071: 'C', 2289072: 'G', 2289081: 'T', 2289082: 'A', 
                         2289089: 'T', 2289090: 'G', 2289091: 'A', 2289097: 'T', 2289100: 'A', 2289103: 'G', 2289111: 'C', 2289133: 'A', 
                         2289138: 'G', 2289140: 'C', 2289150: 'C', 2289162: 'G', 2289171: 'T', 2289180: 'C', 2289186: 'G', 2289193: 'T', 
                         2289200: 'T', 2289202: 'G', 2289203: 'C', 2289206: 'C', 2289207: 'G', 2289214: 'T', 2289218: 'T', 2289219: 'G', 
                         2289220: 'T', 2289222: 'C', 2289223: 'A', 2289225: 'G', 2289231: 'G', 2289234: 'T', 2289239: 'T', 2289240: 'G', 
                         2289248: 'C', 2289252: 'A', 4247429: 'C', 4247430: 'C', 4247431: 'T', 4247729: 'T', 4247730: 'A', 4248003: 'G'}

    list_resistance = []
    
    for index, _ in vcf_df.iterrows():
        position = int(vcf_df.loc[index,'POS'])
        alt_nucleotide = str(vcf_df.loc[index,'ALT'])
        nucleotides = []
        
        
        if position in dict_resistance_position.keys():
            #Check position in resistance dict
            #Create a list with all possible nucleotydes in each position
            if len(dict_resistance_position[position]) == 2:
                nucleotides.append(dict_resistance_position[position][0])
            elif len(dict_resistance_position[position]) > 2:
                nucleotides.append(dict_resistance_position[position][0])
                last_nucleotides = dict_resistance_position[position][2:]
                #print(last_nucleotides)
                for nucleotide in last_nucleotides: #Append extra nucleotides
                    nucleotides.append(nucleotide)
            if alt_nucleotide in nucleotides:
                snp_resist = alt_nucleotide #ALT
                resistance = dict_resistance_position[int(position)][1] #Resist name
                list_resistance.append(str(position)) #POS
                list_resistance.append(snp_resist)
                #Evaluate High confidence
                if (int(position) in dict_high_confidence.keys()) and (dict_high_confidence[int(position)] in nucleotides):
                    list_resistance.append("*")
                    
                    vcf_df.loc[index,'Resistance'] = resistance + "*"
                else:
                    vcf_df.loc[index,'Resistance'] = resistance
                    
            list_resistance.append("\t")
    #list_resistance.append(resistance + "\n")
    
    if len(list_resistance) > 0:
        print("This strain has resistance positions:\n:" + ",".join(list_resistance))
        return ",".join(list_resistance)
    else:
        print("No resistance were found\n")


def final_annotation(vcf_file_annot):
    """
    import annotated vcf with snpEff
    add Lineage info -> output final lineage to an external file
    add resistance info -> output final lineage to an external file
    """
    df_vcf = import_annot_to_pandas(vcf_file_annot, sep='\t')

    vcf_path = os.path.abspath(vcf_file_annot)
    output_dir = ("/").join(vcf_path.split("/")[:-1])
    vcf_name = vcf_path.split("/")[-1]

    tab_name = (".").join(vcf_name.split(".")[:-1])
    #extend_raw = ".raw.annot.tab"
    extend_final = ".annot.tsv"

    
    #Add lineage info 
    add_lineage_Coll(df_vcf)

    #Add resistance info
    add_resistance_snp(df_vcf)


    #Retrieve only annotation fields
    df_vcf_annot = df_vcf[['#CHROM', 'POS', 'ID', 'REF', 'ALT','Annotation',
       'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type',
       'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p',
       'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length', 'AA.pos / AA.length','Lineage', 'Resistance']]
    
    #output all raw info into a file
    new_out_file = tab_name + extend_final
    output_raw_tab = os.path.join(output_dir, new_out_file)
    df_vcf_annot.to_csv(output_raw_tab, sep='\t', index=False)
    

    #final_vcf_name = tab_name + extend_final
    #filter_vcf_list(vcf_path, list_positions_to_filter, final_vcf_name)