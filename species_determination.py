#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np
import re
import subprocess
from misc import check_create_dir, obtain_output_dir, extract_sample

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 09 July 2019
REVISION: 

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

    input_group.add_argument('-r', '--r1_file', metavar="reference", type=str, required=True, help='REQUIRED. File to map against')
    input_group.add_argument('-R', '--r2_file', metavar="sample", type=str, required=True, help='Sample to identify further files')
    
    output_group = parser.add_argument_group('Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    params_group.add_argument('-c', '--mincov', type=int, required=False, default=20, help='Minimun coverage to add samples into analysis')
    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=4, help='Threads to use')
    params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=8, help='MAx memory to use')

    arguments = parser.parse_args()

    return arguments


def zcat_concat_reads(args):
    """
    decompress gz r1 and r2 into a combined .fastq with the common sample in the same directory
    """

    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)

    sample = extract_sample(r1, r2)

    output_dir = ("/").join(r1.split("/")[:-1])
    output_name = sample + ".fastq"
    output_file = os.path.join(output_dir, output_name)

    cmd = ["zcat", r1, r2]
    #execute_subprocess(cmd)
    with open(output_file, "w+") as outfile:
        #calculate coverage and save it in th eoutput file
        subprocess.run(cmd,
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)
    
    return output_file


def mash_screen(args, winner=True, mash_database="/home/laura/DATABASES/Mash/refseq.genomes.k21s1000.msh"):
    #https://mash.readthedocs.io/en/latest/index.html
    #https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh #MASH refseq database
    # mash screen -w -p 4 ../refseq.genomes.k21s1000.msh 4_R1.fastq.gz 4_R2.fastq.gz > 4.winner.screen.tab
    #identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment

    if not os.path.isfile(mash_database):
        print(RED + BOLD + "Mash database can't be found\n" + END_FORMATTING + "You can download it typing:\n\
            wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh")
        sys.exit(1)

    threads = args.threads

    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)

    sample = extract_sample(r1, r2)

    species_output_dir = obtain_output_dir(args, "Species")
    check_create_dir(species_output_dir)
    species_output_name = sample + ".screen.tab"
    species_output_file = os.path.join(species_output_dir, species_output_name)

    cmd = ["mash", "screen", "-p", str(threads), mash_database, r1, r2]

    if winner == True:
        cmd.insert(2,"-w")

    #cmd.extend([mash_database, r1, r2])

    prog = cmd[0]
    param = cmd[1:]

    print(" ".join(cmd))

    try:
    #execute_subprocess(cmd)
        with open(species_output_file, "w+") as outfile:
            #calculate mash distance and save it in output file
            command = subprocess.run(cmd,
            stdout=outfile, stderr=subprocess.PIPE, universal_newlines=True)
        if command.returncode == 0:
            print(GREEN + "Program %s successfully executed" % prog + END_FORMATTING)
        else:
            print (RED + BOLD + "Command %s FAILED\n" % prog + END_FORMATTING
                + BOLD + "WITH PARAMETERS: " + END_FORMATTING + " ".join(param) + "\n"
                + BOLD + "EXIT-CODE: %d\n" % command.returncode +
                "ERROR:\n" + END_FORMATTING + command.stderr)

    except OSError as e:
        sys.exit(RED + BOLD + "failed to execute program '%s': %s" % (prog, str(e)) + END_FORMATTING)



if __name__ == '__main__':
    print("#################### SPECIES #########################")
    args = get_arguments()
    #zcat_concat_reads(args)
    mash_screen(args, winner=True, mash_database="/home/laura/DATABASES/Mash/refseq.genomes.k21s1000.msh")