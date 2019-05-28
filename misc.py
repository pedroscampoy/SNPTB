
import os
import sys
import re
import subprocess
import pandas as pd
import numpy as np


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

def check_file_exists(file_name):
    """
        Check file exist and is not 0 Kb, if not program exit.
    """
    file_info = os.stat(file_name) #Retrieve the file info to check if has size > 0

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        print(RED + BOLD + "File: %s not found or empty\n" % file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)


def check_remove_file(file_name):
    """
    Check file exist and remove it.
    """
    if os.path.exists(file_name):
        os.remove(file_name)
    

def extract_sample(R1_file, R2_file):
    """
    Extract sample from R1, R2 files.
    """
    basename_R1 = os.path.basename(R1_file)
    basename_R2 = os.path.basename(R2_file)

    sample_name_R = os.path.commonprefix([basename_R1, basename_R2])
  
    long_suffix = re.search('_S.*', sample_name_R)
    short_suffix = re.search('_R.*', sample_name_R)
    bar_suffix = re.search('_$', sample_name_R)
    
    if long_suffix:
        match = long_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif short_suffix:
        match = short_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif bar_suffix:
        match = bar_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    else:
        sample_name = sample_name_R

    return sample_name


def obtain_output_dir(args, subfolder=None):
    """
    Returns output folder and output file depending on the output supplied.
    """
    if args.output != None:
        output_dir_arg = os.path.abspath(args.output)
        output_dir = os.path.join(output_dir_arg, subfolder)
    elif args.r1_file:
        r1 = os.path.abspath(args.r1_file)
        output_dir_arg = os.path.dirname(r1)
        output_dir = os.path.join(output_dir_arg, subfolder)
    elif args.input_bam:
        bam = os.path.abspath(args.input_bam)
        output_dir_arg = os.path.dirname(bam)
        output_dir = os.path.join(output_dir_arg, subfolder)
    return output_dir
    
def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def get_picard_path():
    type_route = subprocess.run(["whereis", "picard.jar"],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
    regex = re.compile(r'\/.*\.jar')
    picard_route = re.search(regex, type_route.stdout)

    return picard_route.group()

def execute_subprocess(cmd):
    """
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """
    if cmd[0] == "java":
        prog = cmd[2].split("/")[-1] + " " + cmd[3]
        param = cmd[4:]
    elif cmd[0] == "samtools" or cmd[0] == "bwa" or cmd[0] == "gatk":
        prog = " ".join(cmd[0:2])
        param = cmd[3:]
    else:
        prog = cmd[0]
        param = cmd[1:]
    
    try:
        command = subprocess.run(cmd , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if command.returncode == 0:
            print(GREEN + "Program %s successfully executed" % prog + END_FORMATTING)
        else:
            print (RED + BOLD + "Command %s FAILED\n" % prog + END_FORMATTING
                + BOLD + "WITH PARAMETERS: " + END_FORMATTING + " ".join(param) + "\n"
                + BOLD + "EXIT-CODE: %d\n" % command.returncode +
                "ERROR:\n" + END_FORMATTING + command.stderr.decode().strip())
    except OSError as e:
        sys.exit(RED + BOLD + "failed to execute program '%s': %s" % (prog, str(e)) + END_FORMATTING)


def extract_read_list(input_dir):
    """
    Search files in a directory sort by name and extract comon name of R1 and R2
    with extract_sample() function
    """
    input_dir = os.path.abspath(input_dir)
    r1_list = []
    r2_list = []
    for root, _, files in os.walk(input_dir):
        for name in files:
            filename = os.path.join(root, name)
            is_fasta = re.match(r'.*\.fast[aq](\.gz)*',filename)
            r1 = re.match(r'.*(_R1_|_1|_1_|_R1).*\.fast[aq](\.gz)*',filename)
            r2 = re.match(r'.*(_R2_|_2|_2_|_R2).*\.fast[aq](\.gz)*',filename)

            if is_fasta:
                if r1:
                    r1_list.append(r1.group())
                elif r2:
                    r2_list.append(r2.group())
                else:
                    print(RED + "ERROR, file is not R1 nor R2" + END_FORMATTING)
                    sys.exit(1)
    r1_list = sorted(r1_list)
    r2_list = sorted(r2_list)
    return r1_list, r2_list

    # for filename in sorted(os.listdir(input_dir)):
    #     filename = os.path.abspath(filename)
    #     is_fasta = re.match('.*\.fast[aq](\.gz)*',filename)
    #     r1 = re.match('.*(_R1_|_1_).*\.fast[aq](\.gz)*',filename)
    #     r2 = re.match('.*(_R2_|_2_).*\.fast[aq](\.gz)*',filename)
    #     if is_fasta:
    #         if r1:
    #             r1_list.append(r1.group())
    #         elif r2:
    #             r2_list.append(r2.group())
    #         else:
    #             print(RED + "ERROR, file is not R1 nor R2" + END_FORMATTING)
    #             sys.exit(1)
    
    # return r1_list, r2_list


def extract_sample_list():
    #sample_list = []
    # for r1, r2 in zip(r1_list, r2_list):
    #     sample = extract_sample(r1, r2)
    #     sample_list.append(sample)
    pass

def return_codon_position(number):
    position = number % 3
    if position == 0:
        position = 3
    print("number= %s, pos= %s" % (number,position))

def file_to_list(file_name):
    list_F = []
    file_name_abs = os.path.abspath(file_name)
    with open(file_name_abs, "r") as f:
        for line in f:
            list_F.append(line.strip())
    return list_F


def get_coverage(args, input_bam, output_fmt="-d"):
    """
    #Calculate genome coverage at each position using bedtools and an input bam
    https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
    """
    #reference = os.path.abspath(args.reference)

    input_bam = os.path.abspath(input_bam)
    input_bam_base = os.path.basename(input_bam)

    sample = input_bam_base.split(".")[0]
    output_dir = obtain_output_dir(args, "Coverage")
    sample_name = sample + ".cov"
    output_file = os.path.join(output_dir, sample_name)

    check_create_dir(output_dir)

    #execute_subprocess(cmd)
    with open(output_file, "w") as outfile:
        #calculate coverage and save it in th eoutput file
        subprocess.run(["genomeCoverageBed", "-ibam", input_bam, output_fmt], 
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)

def calculate_cov_stats(file_cov):
    df = pd.read_csv(file_cov, sep="\t", names=["#CHROM", "POS", "COV" ])
    unmmaped_pos = len(df.POS[df.COV == 0].tolist())
    total_pos = df.shape[0]
    unmmaped_prop = "%.2f" % ((unmmaped_pos/total_pos)*100)
    mean_cov = "%.2f" % (df.COV.mean())
    return mean_cov, unmmaped_prop

def obtain_group_cov_stats(directory):
    directory_path = os.path.abspath(directory)
    if directory_path.endswith("Coverage"):
        file_name = directory_path.split("/")[-2]
    else:
        file_name = "samples"

    output_file_name = file_name + ".covegare.tab"
    output_file = os.path.join(directory_path,output_file_name)


    with open(output_file, "w") as outfile:
        outfile.write("#SAMPLE" + "\t" + "MEAN_COV" + "\t" + "UNMMAPED_PROP" + "\n")
        for root, _, files in os.walk(directory_path):
            for name in files:
                filename = os.path.join(root, name)
                file_name_cov = os.path.basename(filename)
                sample = file_name_cov.split(".")[0]
                if filename.endswith(".cov") and (os.path.getsize(filename) > 0):
                    mean_cov, unmmaped_prop = calculate_cov_stats(filename)
                    outfile.write("%s\t%s\t%s\n" % (sample, mean_cov, unmmaped_prop))