import os
import pandas as pd
import argparse
import sys
import subprocess
from sklearn.metrics import jaccard_similarity_score, pairwise_distances
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as ssd #pdist
#import PyQt5
#from PyQt5 import QtGui
#import ete3
#from ete3 import Tree, TreeStyle

from misc import check_file_exists, import_to_pandas
from vcf_process import import_VCF42_cohort_pandas

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

    #Define parser and program
    parser = argparse.ArgumentParser(prog = "VCF_DDTB.py", description= "copare SNP presence in vcf files") 

    #Define subtask/subparsers
    subparsers = parser.add_subparsers( dest = "subtask", help = "new / update / compare / extract commands either add new samples, compare or discard exixting samples")

    new_parser = subparsers.add_parser("new", help = "Create new ddbb with presence/absence of snp")
    update_parser = subparsers.add_parser("update", help = "Add new sample using a list of variants, files supplied or files on folder")
    compare_parser = subparsers.add_parser("compare", help = "Comapare samples supplied or all samples to obtain a pirwise matrix")
    extract_parser = subparsers.add_parser("extract", help = "Remove samples supplied from databse")


    #new_parser
    new_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the newd ddbb")
    new_parser.add_argument("-v", "--vcf", required= False, action='store_true', help="Database will use vcf files instead of ")
    new_parser.add_argument("-s", "--suffix", required= False, type=str, default=".SNP.final.vcf", help="Suffix to filter within vcf with similar suffix")
    new_parser.add_argument("-r", "--recalibrate", required= False, type=str, default=False, help="Folder with tab files to asses discrepancies after comparing")


    new_exclusive = new_parser.add_mutually_exclusive_group(required= True)

    new_exclusive.add_argument("-F", "--folder",  dest = "folder", metavar="folder",required= False, type=str, help="Folder containinig files with snp positions")
    new_exclusive.add_argument("-f", "--file",  dest = "single_file", metavar="file[s]", required= False, type=str, help="individual files with snp positions")
    new_exclusive.add_argument("-l", "--list",  dest = "snp_list", metavar="list", required= False, help="file with a list of positions in a column, not vcf format")


    #update_parser
    update_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the updated ddbb")
    update_parser.add_argument("-s", "--snp-final", required= False, action='store_true', help="Database will use snp.fila instead of gatk vcf ")
    update_parser.add_argument("-d", "--database",  dest = "update_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")

    update_exclusive = update_parser.add_mutually_exclusive_group()

    update_exclusive.add_argument("-F", "--folder",  dest = "folder", metavar="folder",required= False, type=str, help="Folder containinig files with snp positions")
    update_exclusive.add_argument("-f", "--file",  dest = "single_file", metavar="file[s]", required= False, type=str, help="individual files with snp positions")
    update_exclusive.add_argument("-l", "--list",  dest = "snp_list", metavar="list", required= False, help="file with a list of positions in a column, not vcf format")

    update_parser.add_argument("-b", "--backup",  dest = "backup", action="store_true", help="Creates an aditional database with the date as backup")

    #compare_parser
    compare_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")
    compare_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= False, metavar="filename", help="REQUIRED: file name, including PATH, of the matrix comparison")

    compare_exclusive = compare_parser.add_mutually_exclusive_group()

    compare_exclusive.add_argument("-a", "--all",  dest = "all_compare", action="store_true", required= False, help="All files in supplied database will be compared")
    compare_exclusive.add_argument("-s", "--samples",  dest = "samples_compare", metavar="sample_name[s]", nargs="+", required= False, help="Sample names supplied will be compared")

    #extract_parser
    extract_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database with the sample to remove")
    extract_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the updated ddbb")

    extract_parser.add_argument("-s", "--samples",  dest = "samples_extract", metavar="sample name[s]", nargs="+", required= True, help="Sample names supplied will be removed")

    parser.add_argument("--version", action="version", version="%(prog)s 0.1")

    arguments = parser.parse_args()

    return arguments

def blank_database():
    new_pandas_ddtb = pd.DataFrame(columns=['Position','N', 'Samples'])
    return new_pandas_ddtb

def import_VCF4_to_pandas(vcf_file, sep='\t'):
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
        dataframe['POS'] = dataframe['POS'].astype(int)
        
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe

def recheck_variant(format_sample):
    #GT:AD:DP:GQ:PGT:PID:PL:PS
    list_format = format_sample.split(":")
    gt = list_format[0]
    #gt0 = gt[0]
    #gt1 = gt[1]
    ad = list_format[1]
    ref = int(ad.split(',')[0])
    alt = max(int(x) for x in ad.split(',')[0:])
    
    if gt == "0/0":
        value = 0
    elif gt == "1/1":
        value = 1
    else:
        if gt == "./.":
            value = "!"
        elif "2" in gt:
            value = "!"
        elif (ref > alt):
            value = 0
        elif (alt > ref):
            value = 1
        else:
            value = "!"
            
    return value

def recheck_variant_mpileup(reference_file, position, sample, bam_folder):
    #Find reference name
    with open(reference_file) as f:
        reference = f.readline().split(" ")[0].strip(">").strip()
    #Identify correct bam
    for root, _, files in os.walk(bam_folder):
        for name in files:
            filename = os.path.join(root, name)
            if name.startswith(sample) and name.endswith(".bqsr.bam"):
                bam_file = filename
    #format position for mpileuo execution (NC_000962.3:632455-632455)
    position = reference + ":" + str(position) + "-" + str(position)
    
    #Execute command and retrieve output
    cmd = ["samtools", "mpileup", "-f", reference_file, "-aa", "-r", position, bam_file]
    text_mpileup = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
    
    #Extract 5th column to find variants
    variant = text_mpileup.stdout.split()[4]
    var_list = list(variant)
    
    most_freq_var = max(set(var_list), key = var_list.count).upper()
        
    if most_freq_var == "." or most_freq_var == "," or most_freq_var == "*":
        return 0
    else:
        return 1

def identify_nongenotyped_mpileup(reference_file, row_position, sample_list_matrix, list_presence, bam_folder):
    """
    Replace nongenotyped ("!") with the most abundant genotype
    """
    #mode = max(set(list_presence), key = list_presence.count)
    
    count_ng = list_presence.count("!")
    sample_number = len(list_presence)
    
    if "!" not in list_presence:
        return list_presence
    elif count_ng/sample_number > 0.2:
        return 'delete'
    else:
        indices_ng = [i for i, x in enumerate(list_presence) if x == "!"]
        for index in indices_ng:
            #print(reference_file, row_position, sample_list_matrix[index], bam_folder)
            list_presence[index] = recheck_variant_mpileup(reference_file, row_position, sample_list_matrix[index], bam_folder)
        #new_list_presence = [mode if x == "!" else x for x in list_presence]
        return list_presence

def extract_recalibrate_params(pipeline_folder):
    for root, dirs, _ in os.walk(pipeline_folder):
        if root == pipeline_folder:
            for directory in dirs:
                subfolder = os.path.join(root, directory)
                if subfolder.endswith("/VCF"):
                    for file in os.listdir(subfolder):
                        if file.endswith("cohort.combined.hf.vcf"):
                            cohort_file = os.path.join(subfolder, file)
                            
                            with open(cohort_file, 'r') as f:
                                for line in f:
                                    if line.startswith("#"):
                                        if "--reference " in line:
                                            reference_file = line.split("--reference ")[1].strip().split(" ")[0].strip()
                                        
                            
                elif subfolder.endswith("/Bam"):
                    bam_folder = subfolder
                    
    return (cohort_file, bam_folder, reference_file)

def recalibrate_ddbb_vcf(snp_matrix_ddbb, vcf_cohort, bam_folder, reference_file):
    
    vcf_cohort = os.path.abspath(vcf_cohort)
    #snp_matrix_ddbb = os.path.abspath(snp_matrix_ddbb)
    
    df_matrix = snp_matrix_ddbb
    df_cohort = import_VCF42_cohort_pandas(vcf_cohort)
    
    sample_list_matrix = df_matrix.columns[3:]
    n_samples = len(sample_list_matrix)
    #sample_list_cohort = df_cohort.columns.tolist()[9:]
    
    list_index_dropped = []
    #Iterate over non unanimous positions 
    for index, data_row in df_matrix[df_matrix.N < n_samples].iloc[:,3:].iterrows():
        #Extract its position
        row_position = int(df_matrix.loc[index,"Position"])
        #print(data_row.values)
        #Use enumerate to retrieve column index (column ondex + 3)
        presence_row = [recheck_variant(df_cohort.loc[df_cohort.POS == row_position, df_matrix.columns[n + 3]].item()) \
                           for n,x in enumerate(data_row)]
        #print(presence_row, row_position)
        #Resolve non genotyped using gvcf files
        new_presence_row = identify_nongenotyped_mpileup(reference_file, row_position, sample_list_matrix, presence_row, bam_folder)
        
        #find positions with 20% of nongenotyped and delete them OR
        #reasign positions without nongenotyped positions 
        if new_presence_row == 'delete':
            list_index_dropped.append(index)
        else:
            df_matrix.iloc[index, 3:] = new_presence_row
            df_matrix.loc[index, 'N'] = sum(new_presence_row)
        #print(new_presence_row)
        #print("\n")
    #Remove all rows at once to avoid interfering with index during for loop
    df_matrix.drop(index=list_index_dropped, axis=0, inplace=True)
    
    return df_matrix


def ddtb_add(input_folder, output_filename, recalibrate=False):
    directory = os.path.abspath(input_folder)
    output_filename = os.path.abspath(output_filename)

    final_ddbb = blank_database()

    #Make sure output exist to force change name
    if os.path.isfile(output_filename):
        print(YELLOW + "ERROR: " + BOLD + "output database EXIST, choose a different name or manually delete" + END_FORMATTING)
        sys.exit(1)

    print("Previous final database contains %s rows and %s columns\n" % final_ddbb.shape)
    print("The directory selected is: %s" % directory)
    

    all_samples = 0
    new_samples = 0
    for filename in os.listdir(directory):
        if not filename.startswith('.') and filename.endswith(".combined.hf.SNP.final.vcf"):
            print("\nThe file is: %s" % filename)
            
            all_samples = all_samples + 1
            positions_shared = []
            positions_added = []
            
            sample = filename.split(".")[0] #Manage sample name
            
            file = os.path.join(directory, filename) #Whole file path
            check_file_exists(file) #Manage file[s]. Check if file exist and is greater than 0

            new_sample = import_VCF4_to_pandas(file) #Import files in annotated vcf format

            #Handle each new_sample
            #print("This file contains %s SNPs" % len(new_sample.index))
            
            #Check if sample exist
            ######################
            if sample not in final_ddbb.columns.tolist():
                print("Adding new sample %s to %s" % (sample, os.path.basename(output_filename)))
                new_samples = new_samples + 1
                new_colum_index = len(final_ddbb.columns) #extract the number of columns to insert a new one
                #final_ddbb[sample] = sample #adds a new column but fills all blanks with the value sample
                final_ddbb.insert(new_colum_index, sample, 0) #add a new column with defauls values = 0
                
                #Check if position exist
                ########################
                for position in new_sample['POS'].unique(): #extract first column in file
                    
                    if position not in final_ddbb["Position"].values:
                        positions_added.append(position) #Count new positions for stats
                        
                        new_row = len(final_ddbb.index)
                        final_ddbb.loc[new_row,'Position'] = position
                        final_ddbb.loc[new_row,'Samples'] = sample
                        final_ddbb.loc[new_row,'N'] = int(1)
                        final_ddbb.loc[new_row,sample] = str(1)
                    else:
                        positions_shared.append(position) #Count shared positions for stats
                        
                        #Check whether the column matches the value and retrieve the first position [0]
                        #of the object index generated
                        index_position = final_ddbb.index[final_ddbb["Position"] == position][0]
                        #Add sample to corresponding cell [position, samples]
                        number_samples_with_position = final_ddbb.loc[index_position,'N']
                        names_samples_with_position = final_ddbb.loc[index_position,'Samples']
                        new_names_samples = names_samples_with_position + "," + sample
                        #Sum 1 to the numbes of samples containing the position
                        final_ddbb.loc[index_position,'N'] = number_samples_with_position + 1
                        final_ddbb.loc[index_position,'Samples'] = new_names_samples
                        final_ddbb.loc[index_position,sample] = str(1) #Add "1" in cell with correct position vs sample (indicate present)

                print("\nSAMPLE:\t%s\nTOTAL Variants:\t%s\nShared Variants:\t%s\nNew Variants:\t%s\n"
                % (sample, len(new_sample.index), len(positions_shared), len(positions_added)))
            else:
                print(YELLOW + "The sample " + sample + " ALREADY exist" + END_FORMATTING)

    final_ddbb = final_ddbb.fillna(0).sort_values("Position") #final_ddbb = final_ddbb["Position"].astype(int)
    final_ddbb['N'] = final_ddbb['N'].astype(int)
    final_ddbb = final_ddbb.reset_index(drop=True)

    print("Final database now contains %s rows and %s columns" % final_ddbb.shape)
    if recalibrate == False:
        output_filename = output_filename + ".tsv"
        final_ddbb.to_csv(output_filename, sep='\t', index=False)
    else:
        recalibrate = os.path.abspath(recalibrate)
        if os.path.exists(recalibrate):
            recalibrate_params = extract_recalibrate_params(recalibrate)
            print("\n" + MAGENTA + "Recalibration selected" + END_FORMATTING)
            print(output_filename)
            output_filename = output_filename + ".revised.tsv"
            final_ddbb_revised = recalibrate_ddbb_vcf(final_ddbb, recalibrate_params[0], recalibrate_params[1], recalibrate_params[2])
            final_ddbb_revised.to_csv(output_filename, sep='\t', index=False)
        else:
            print("The directory supplied for recalculation does not exixt")
            sys.exit(1)
    print(output_filename)

    #Create small report with basic count
    #####################################
            
    print("\n" + GREEN + "Position check Finished" + END_FORMATTING)
    print(GREEN + "Added " + str(new_samples) + " samples out of " + str(all_samples) + END_FORMATTING + "\n")


    ###########################COMPARE FUNCTIONS#####################################################################

def compare_snp_columns(sample1, sample2, df):
    jaccard_similarity = jaccard_similarity_score(df[sample1], df[sample2]) #similarities between colums
    hamming_similarity = 1 - jaccard_similarity #disagreements between colums
    snp_distance = int(hamming_similarity * (len(df.index)+1))
    return snp_distance

def snp_distance_pairwise(dataframe, output_file):
    if os.path.exists(output_file):
        os.remove(output_file)
    with open(output_file, "a") as f:
        for sample1 in dataframe.iloc[:,3:].columns: #remove first 3 colums
            for sample2 in dataframe.iloc[:,3:].columns:
                if sample1 != sample2:
                    snp_distance = compare_snp_columns(sample1, sample2, dataframe)
                    line_distance = "%s\t%s\t%s\n" % (sample1, sample2, snp_distance)
                    f.write(line_distance)

def snp_distance_matrix(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    snp_distance_df = pd.DataFrame(hamming_distance * len(dataframe_only_samples.index), index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    snp_distance_df = snp_distance_df.astype(int)
    snp_distance_df.to_csv(output_file, sep='\t', index=True)

def hamming_distance_matrix(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    hamming_distance_df = pd.DataFrame(hamming_distance, index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    hamming_distance_df.to_csv(output_file, sep='\t', index=True)

def clustermap_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    sns.clustermap(dataframe_only_samples, annot=False, cmap="YlGnBu", figsize=(13, 13))
    plt.savefig(output_file, format="png")

def dendogram_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average') #method='single'

    plt.rcParams['lines.linewidth'] = 8 #Dendrogram line with
    plt.rcParams['xtick.major.size'] = 10 #Only affect to tick (line) size
    plt.rcParams.update({'font.size': 30}) #Increase x tick label size
    #plt.tick_params(labelsize=30)
    plt.figure(figsize=(30, 50))
    plt.ylabel('samples', fontsize=30)
    plt.xlabel('snp distance', fontsize=30)

    shc.dendrogram(Z, labels=labelList, orientation='left', distance_sort='descending', show_leaf_counts=True, color_threshold=10, leaf_font_size=20)

    
    plt.savefig(output_file, format="png")

# Convert dendrogram to Newick
def linkage_to_newick(dataframe, output_file):
    """
    Thanks to https://github.com/biocore/scikit-bio/issues/1579
    Input :  Z = linkage matrix, labels = leaf labels
    Output:  Newick formatted tree string
    """
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average')

    tree = shc.to_tree(Z, False)
    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            #print("%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick))
            return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = f"):{(parentdist - node.dist)/2}{newick}"
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = buildNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            #print(newick)
            return newick

    with open(output_file, 'w') as f:
        f.write(buildNewick(tree, "", tree.dist, labelList))
    return buildNewick(tree, "", tree.dist, labelList)


def ddtb_compare(final_database):

    database_file = os.path.abspath(final_database)
    check_file_exists(database_file)
    presence_ddbb = import_to_pandas(database_file, header=True)

    output_path = database_file.split(".")[0]

    print("Output path is: " + output_path)


    print(BLUE + BOLD + "Comparing all samples in " + database_file + END_FORMATTING)
    prior_pairwise = datetime.datetime.now()

    #Calculate pairwise snp distance for all and save file
    print(CYAN + "Pairwise distance" + END_FORMATTING)
    pairwise_file = output_path + ".snp.pairwise.tsv"
    snp_distance_pairwise(presence_ddbb, pairwise_file)
    after_pairwise = datetime.datetime.now()
    print("Done with pairwise in: %s" % (after_pairwise - prior_pairwise))

    #Calculate snp distance for all and save file
    print(CYAN + "SNP distance" + END_FORMATTING)
    snp_dist_file = output_path + ".snp.tsv"
    snp_distance_matrix(presence_ddbb, snp_dist_file)

    #Calculate hamming distance for all and save file
    print(CYAN + "Hamming distance" + END_FORMATTING)
    hmm_dist_file = output_path + ".hamming.tsv"
    hamming_distance_matrix(presence_ddbb, hmm_dist_file)
    """
    #Represent pairwise snp distance for all and save file
    print(CYAN + "Drawing distance" + END_FORMATTING)
    prior_represent = datetime.datetime.now()
    png_dist_file = output_path + ".snp.distance.png"
    #clustermap_dataframe(presence_ddbb, png_dist_file)
    after_represent = datetime.datetime.now()
    print("Done with distance drawing in: %s" % (after_represent - prior_represent))
    """
    #Represent dendrogram snp distance for all and save file
    print(CYAN + "Drawing dendrogram" + END_FORMATTING)
    png_dend_file = output_path + ".snp.dendrogram.png"
    dendogram_dataframe(presence_ddbb, png_dend_file)


    #Output a Newick file distance for all and save file
    print(CYAN + "Newick dendrogram" + END_FORMATTING)
    newick_file = output_path + ".nwk"
    linkage_to_newick(presence_ddbb, newick_file)


if __name__ == '__main__':
    print("#################### COMPARE SNPS #########################")