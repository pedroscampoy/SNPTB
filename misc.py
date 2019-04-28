
import os
import sys
import re
import subprocess


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

def check_file_exists(file_name):
    """
    Check file exist, if not program exit.
    """
    if not os.path.isfile(file_name):
        print(RED + BOLD + "File: %s not found\n" % file_name + END_FORMATTING)
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
        sample_name = sample_name_R.rstrip(match)
    elif short_suffix:
        match = short_suffix.group()
        sample_name = sample_name_R.rstrip(match)
    elif bar_suffix:
        match = bar_suffix.group()
        sample_name = sample_name_R.rstrip(match)
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
            print(GREEN + BOLD + "Program %s successfully executed" % prog + END_FORMATTING)
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
            r1 = re.match(r'.*(_R1_|_1_).*\.fast[aq](\.gz)*',filename)
            r2 = re.match(r'.*(_R2_|_2_).*\.fast[aq](\.gz)*',filename)
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


extract_read_list("/home/pedro/analysis/Mixed/RAW/")