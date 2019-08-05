#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np
import re
import subprocess
from tabulate import tabulate
from misc import get_snpeff_path, check_create_dir
from vcf_process import calculate_ALT_AD, obtain_output_dir

##Import files containing annotation info and convert them to dictionary
##TO DO: Compensatory mutations
script_dir = os.path.dirname(os.path.realpath(__file__))

annotation_dir = os.path.join(script_dir, "annotation/genes")
annotation_dir_res = os.path.join(script_dir, "annotation/resistance")

essential_file = os.path.join(annotation_dir, "dict_locus_essential.txt")
product_file = os.path.join(annotation_dir, "dict_locus_product.txt")
resistance_file_V1 = os.path.join(annotation_dir_res, "dict_position_resistance_v1.txt")
resistance_file_v2 = os.path.join(annotation_dir_res, "dict_position_resistance_v2_inf.txt") #Crated on 190718


dict_essential = {}
dict_product = {}
dict_res_v1 = {}
dict_res_v2 = {}

#Create dictionary of essential genes
with open(essential_file, 'r') as f:
    for line in f:
        dict_essential[line.split(":")[0]] = line.split(":")[1].strip()

#Create dictionary of gene products
with open(product_file, 'r') as f:
    for line in f:
        dict_product[line.split(":")[0]] = (":").join(line.split(":")[1:]).strip()

#Create dictionary of resistance positions
with open (resistance_file_V1, 'r') as f:
    for line in f:
        dict_res_v1[int(line.split(":")[0])] = line.split(":")[1].strip().split(",")

with open (resistance_file_v2, 'r') as f:
    for line in f:
        dict_res_v2[int(line.split(":")[0])] = line.split(":")[1].strip().split(",")


def annotate_bed_file(dict_position, position):
    """
    Identify a position within a range
    credits: https://stackoverflow.com/questions/6053974/python-efficiently-check-if-integer-is-within-many-ranges
    """
    #dict_position = bed_to_dict(bed_file)
    if any(start <= position <= end for (start, end) in dict_position.items()):
        return True
    else:
        return False


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

        dataframe['ALT_AD'] = dataframe.apply(calculate_ALT_AD, axis=1)
        dataframe[['gt0','gt1']] = dataframe['GT'].str.split(r'[/|\|]', expand=True)

        dataframe['HGVS.c'] = dataframe['HGVS.c'].str.split(".").str[-1]
        dataframe['HGVS.p'] = dataframe['HGVS.p'].str.split(".").str[-1]
        dataframe['Gene length'] = dataframe['CDS.pos / CDS.length'].str.split("/").str[-1]
        dataframe['AA length'] = dataframe['AA.pos / AA.length'].str.split("/").str[-1]
        dataframe['AA length'] = dataframe['AA length'].replace('', 0)
        dataframe['HGVS.p'] = dataframe['HGVS.p'].replace('', 'None')
                
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
        '931123' : ['T', '4'],
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
        '1759252' : ['G', '4.9'],
        '1799921' : ['A', '5'],
        '1816587' : ['G', '6'],
        '1137518' : ['A', '7'],
        '2831482' : ['G', 'BOV'],
        '1882180' : ['T', 'BOV_AFRI']
                }
    list_lineage = []
    
    vcf_df['Lineage'] = np.nan

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
                if str(list_lineage[sublineage_n]).startswith(str(list_lineage[sublineage_n + 1])):
                    asterix = asterix + "*"
        final_lineage = list_lineage[0] + " " + asterix
        print("This strain has lineage position(s):\n: " + " ".join([list_lineage[0],asterix]))
        return final_lineage
    else:
        print("No lineage were found\n")


def add_resistance_snp(vcf_df, dict_resistance_position=dict_res_v1):
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
    
    vcf_df['Resistance'] = np.nan

    for index, _ in vcf_df.iterrows():
        position = int(vcf_df.loc[index,'POS'])
        alt_nucleotide = str(vcf_df.loc[index,'ALT'])
        
        if position in dict_resistance_position.keys():
            #Check position in resistance dict
            #Create a list with all possible nucleotydes in each position
            nucleotides = dict_resistance_position[position][1:]

            if alt_nucleotide in nucleotides:
                snp_resist = alt_nucleotide #ALT
                resistance = dict_resistance_position[int(position)][0] #Resist name
                list_resistance.append(resistance)
                list_resistance.append(str(position)) #POS
                list_resistance.append(snp_resist)
                #Evaluate High confidence
                if (int(position) in dict_high_confidence.keys()) and (dict_high_confidence[int(position)] == alt_nucleotide):
                    list_resistance.append("*")
                    
                    vcf_df.loc[index,'Resistance'] = resistance + "*"
                else:
                    vcf_df.loc[index,'Resistance'] = resistance
                    
            list_resistance.append("\t")
    #list_resistance.append(resistance + "\n")
    
    if len(list_resistance) > 3:
        print("This strain has resistance positions:\n:" + ",".join(list_resistance))
        return ",".join(list_resistance)
    else:
        print("No resistance were found\n")


def add_essential_cateory(row, dict_essential=dict_essential):
    if row.Gene_ID in dict_essential.keys():
        if dict_essential[row.Gene_ID] == "essential":
            return "essential"
        else:
            return "nonessential"

def add_product_cateory(row, dict_product=dict_product):
    if row.Gene_ID in dict_product.keys():
        return dict_product[row.Gene_ID]


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

    #Add essential info
    df_vcf['Is_essential'] = df_vcf.apply(add_essential_cateory, axis=1)

    #Add protein product
    df_vcf['Product'] = df_vcf.apply(add_product_cateory, axis=1)

    #Add lineage info 
    add_lineage_Coll(df_vcf)

    #Add resistance info
    add_resistance_snp(df_vcf)


    #Retrieve only annotation fields
    df_vcf_annot = df_vcf[['#CHROM', 'POS', 'ID', 'REF', 'ALT','Annotation',
       'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type',
       'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p',
       'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length', 'AA.pos / AA.length','Is_essential','Product', 'Lineage', 'Resistance']]
    
    #output all raw info into a file
    new_out_file = tab_name + extend_final
    output_raw_tab = os.path.join(output_dir, new_out_file)
    df_vcf_annot.to_csv(output_raw_tab, sep='\t', index=False)
    

def get_reverse(nucleotyde):
    nucleotyde = str(nucleotyde)
    nucleotyde_rev = {'A' : 'T',
                     'T' : 'A',
                     'C' : 'G',
                     'G': 'C'}
    if len(nucleotyde) > 1:
        nucleotyde_str = nucleotyde[::-1] #Reverse nucleotide
        nucleotyde_str_fin = "".join([nucleotyde_rev[x] for x in nucleotyde_str]) #Complement nucleotide
        return nucleotyde_str_fin
    else:
        return nucleotyde_rev[nucleotyde]

css_report = """

    <style type="text/css">

    body {
        font: normal 20px Verdana, Arial, sans-serif;
    }

    table {
        text-align: center;
        border-color: #000;
        border-spacing: 0px;
        border-style: solid;
        border-width: 1px;
    }

    th, td {
    border-bottom: 1px solid #ddd;
    padding:0 5px 0 5px;
    }

    th {
    background-color: rgb(76, 175, 170);
    }

    tr:nth-child(even) {background-color: #cecccc;}
    tr:hover {background-color:#7c7b7b;}

    </style>

    """

def create_report(tab_annot, css=css_report, species="Mycobacterium tuberculosis", species_report="Main species: <i>Mycobacterium tuberculosis</i><br />"):
    #<div style="position: absolute; bottom: 5px; color: red; background-color: rgb(253, 253, 253)">
    #Text disclaimer 
    #</div>

    script_dir = os.path.dirname(os.path.realpath(__file__))
    annotation_dir = os.path.join(script_dir, "annotation/resistance")

    table_res = os.path.join(annotation_dir, "MTB_Resistance_Mediating.txt")
    df_res = pd.read_csv(table_res, sep="\t", header=0)
    df_res['High Confidence SNP'].fillna("no", inplace=True)


    output = os.path.dirname(tab_annot)
    sample = os.path.basename(tab_annot).split(".")[0]
    report_name = sample + ".annot.report.html"

    output_file = os.path.join(output, report_name)
    
    cummulative_report = ""

    with open(output_file, 'w+') as f:

        f.write(css)

        line_sample = "Sample name: " + sample + "<br /><br />"
        f.write(line_sample)
        cummulative_report = cummulative_report + line_sample

        line_species = "Species: " + "<i>" + str(species) + "</i>" + "<br /><br />"
        f.write(line_species)
        cummulative_report = cummulative_report + species_report + "<br />"

        df_annot = pd.read_csv(tab_annot, sep="\t", header=0)
        
        df_resistance = df_annot[df_annot.Resistance.notnull()]
        df_resistance_F = df_resistance[['POS', 'ALT', 'Annotation', 'Gene_ID', 'Gene_Name', 'HGVS.c', 'HGVS.p', 'Resistance']]
        df_resistance_F.columns = ['Position', 'Alt. base', 'Change', 'Gene ID', 'Gene Name', 'Codon change','AA change', 'Resistance']
        list_resistance = df_annot['Resistance'][df_annot.Resistance.notnull()].tolist()


        list_lineage = df_annot['Lineage'][df_annot.Lineage.notnull()].tolist()
        
        #Output Lineage info
        if len(list_lineage) > 0:
            list_lineage.sort(reverse=True)
            asterix = ""
            for sublineage_n in range(len(list_lineage)):
                if sublineage_n < (len(list_lineage) - 1):
                    if str(list_lineage[sublineage_n]).startswith(str(list_lineage[sublineage_n + 1])):
                        asterix = asterix + "*"
            final_lineage = str(list_lineage[0]) #+ " " + asterix
            line_lineage = "This strain has lineage position(s): " + "<b>" + str(final_lineage) + "</b>" + "<br /><br />"
            f.write(line_lineage)
            cummulative_report = cummulative_report + line_lineage
        else:
            line_lineage = "No lineage positions were found<br /><br />"
            f.write(line_lineage)
            cummulative_report = cummulative_report + line_lineage
        
        #Output Resistance info
        if len(list_resistance) > 0:
            line_res_1 = "This strain has " + str(len(list_resistance)) + " resistance position(s):<br />"
            f.write(line_res_1)
            cummulative_report = cummulative_report + line_res_1
            """
            additional_resistance = []
            final_res_table = pd.DataFrame(columns= df_res.columns.tolist())

            for index, _ in df_annot[df_annot.Resistance.notnull()].iterrows():
                position = str(df_annot.loc[index,'POS'])
                resistance_name = df_annot.loc[index,'Resistance'].strip("*")
                if df_annot.loc[index,'Gene_ID'].endswith("c"):
                    alt_nucleotide = get_reverse(df_annot.loc[index,'ALT']).lower()
                else:
                    alt_nucleotide = df_annot.loc[index,'ALT'].lower()
                    
                    
                if position in df_res['Variant position genome stop'].values.tolist():
                    row = df_res[(df_res['Var. base'] == alt_nucleotide) & (df_res['Variant position genome stop'] == position)]
                    index = row.index[0]
                    final_res_table = final_res_table.append(df_res.iloc[index], ignore_index=True)
                    #df_resistance_F.drop(df_resistance_F.iloc[index])
                else:
                    #df_resistance_F.drop(df_resistance_F.index[index])
                    other_resistances = position + " " + alt_nucleotide + " " + resistance_name
                    additional_resistance.append(other_resistances)
                
            
            final_res_table.reset_index(drop=True, inplace=True)
            final_res_table_F = final_res_table[['Variant position genome stop', 'Var. base',
                                                'Region', 'Gene ID', 'Gene Name',
                                                'AA change', 'Codon change', 'Antibiotic'
                                                ]]
            #df.rename(columns={'oldName1': 'newName1', 'oldName2': 'newName2'}, inplace=True)
            final_res_table_F.columns = ['Position', 'Alt. base',
                                                'Region', 'Gene ID', 'Gene Name',
                                                'AA change', 'Codon change', 'Antibiotic'
                                                ]
            f.write(tabulate(final_res_table_F, headers='keys', tablefmt='html', showindex=False))
            if len(additional_resistance) > 0:
                f.write("<br /><br />")
                f.write("Found other putative resistances:<br />")
                f.write(tabulate(df_resistance_F, headers='keys', tablefmt='html', showindex=False))
                line_other_res = ("<br />").join(additional_resistance)
                f.write(line_other_res)
            
            """
            f.write(tabulate(df_resistance_F, headers='keys', tablefmt='html', showindex=False))
            table_res = tabulate(df_resistance_F, headers='keys', tablefmt='html', showindex=False)
            cummulative_report = cummulative_report + table_res + "\n"

        else:
            f.write("No Resistance positions were found<br />")
            cummulative_report = cummulative_report + "No Resistance positions were found<br />\n"

        f.write("<br /><br />Este informe debe ser utilizado exclusivamente con fines de investigación. No utilizar ningún dato con fines asistenciales.<br />")

    return cummulative_report

if __name__ == '__main__':
    print("#################### ANNOTATION #########################")