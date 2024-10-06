#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 4/1/2024
# Purpose: Run Rfam to detect non-coding RNA (ncRNAs)

# %% 
import argparse
import os
import subprocess
from checkd_software import cmscan, esl_seqstat, Rfam_cm, Rfam_clanin

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-fas_file", "--fas_file", metavar="path to assembly fasta file", help="fasta file", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)

args = parser.parse_args()
fas_file = args.fas_file
output_dir = args.output_dir
mag_name = args.assembly_name

seqstat_file_out = f'{output_dir}/{mag_name}.seqstat'
cmscan_file_out = f'{output_dir}/{mag_name}.cmscan'
tblout_file_out = f'{output_dir}/{mag_name}.tblout'
tblout_deovlp_file_out = f'{output_dir}/{mag_name}.deoverlapped.tblout'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 
# # Test Vars
# mag_name = 'bin.76'
# output_dir = '/home/kmorten/humanO1_CheckD/hifiasm_MAGS/rfam_out'
# fas_file = f'/home/kmorten/hifi/NM2022/humanO1/MAGs//{mag_name}.fa'
# seqstat_file_out = f'{output_dir}/{mag_name}.seqstat'
# cmscan_file_out = f'{output_dir}/{mag_name}.cmscan'
# tblout_file_out = f'{output_dir}/{mag_name}.tblout'
# tblout_deovlp_file_out = f'{output_dir}/{mag_name}.deoverlapped.tblout'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# %% 
'''Rfam Sequence Stats'''

cmd = f'{esl_seqstat} {fas_file} \
    > {seqstat_file_out}'
subprocess.check_output(cmd, shell=True) 

def find_totMB(seqstat_file) :
    with open(seqstat_file, 'r') as fp:
        totMB = 'NA'
        MB = 10**6
        for line in fp:
            if line.startswith('Total # residues:'):
                print(line)
                linep = line.strip().split(':')
                residues = int(linep[1].strip(' '))
                totMB = float((residues*2)/MB)
                print(f'Total Megabases: {totMB}')
    return totMB

totMB = find_totMB(seqstat_file_out)

'''Rfam cmscan'''

cmd = f'{cmscan} -Z {totMB} \
    --cut_ga --rfam --nohmmonly \
    --tblout {tblout_file_out} \
    --fmt 2 --clanin {Rfam_clanin} \
    {Rfam_cm} {fas_file} \
    > {cmscan_file_out}'
subprocess.check_output(cmd, shell=True) 


'''Removing lower-scoring overlaps from a tblout file'''

cmd = f'grep -v " = " {tblout_file_out} > {tblout_deovlp_file_out}'
subprocess.check_output(cmd, shell=True) 


