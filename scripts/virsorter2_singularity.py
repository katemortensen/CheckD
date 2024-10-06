#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 8/30/2024
# Purpose: Run VirSorter2 with singularity to detect viral DNA in metagenomes 


# %% 
import argparse
import os
import subprocess

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-fasta_file", "--fasta_file", metavar="path to assembly fasta file", help="fasta file", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)
parser.add_argument("-min_length", "--min_length", metavar="min_length", help="min_length flag set for virsorter2", required=False)
parser.add_argument("-threads", "--threads", metavar="threads", help="threads", required=False)

args = parser.parse_args()

if args.threads:
    threads = int(args.threads)
else:
    threads = 1
if args.min_length:
    min_length = int(args.min_length)
else:
    min_length = 1500
    
fasta_file = args.fasta_file
fasta_file_name = fasta_file.split('/')[-1]
output_dir = args.output_dir
assembly_name = args.assembly_name
mag_dir = '/'.join(fasta_file.strip().split('/')[:-1])

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
if not os.path.exists(f'{output_dir}/virsorter2.sif'):
    cmd = f'''
    cd {output_dir}
    module load apptainer
    which singularity
    apptainer build virsorter2.sif docker://jiarong/virsorter:latest
    '''
    subprocess.run(cmd, shell=True, executable='/bin/bash')
    
# delete this hardcoded comment before running


# %%
# humanO1

# assembly_name = 'bin.0'
# wkdir =  f'/home/kmorten/humanO1_CheckD'
# output_dir = f'{wkdir}/virsorter2_out/'
# fasta_file = f'/home/kmorten/hifi/NM2022/humanO1/MAGs//{assembly_name}.fa'
# fasta_file_name = fasta_file.split('/')[-1]
# mag_dir = '/'.join(fasta_file.strip().split('/')[:-1])
# min_length = 1500
# threads = 8 

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# if not os.path.exists(f'{output_dir}/virsorter2.sif'):
#     cmd = f'''
#     cd {output_dir}
#     module load apptainer
#     which singularity
#     apptainer build virsorter2.sif docker://jiarong/virsorter:latest
#     '''
#     subprocess.run(cmd, shell=True, executable='/bin/bash')
    
# %% 
cmd = f'''
cd {output_dir}
module load apptainer
which singularity
singularity run -B {mag_dir}:/mnt/host_files {output_dir}/virsorter2.sif run -w {output_dir}/{assembly_name} -i /mnt/host_files/{fasta_file_name} --min-length {min_length} -j {threads} all
'''

print(cmd)

try:
    subprocess.run([
        f'{cmd}'
    ], shell=True, executable='/bin/bash', check=True, capture_output=True, text=True)
    print(f'done')
except subprocess.CalledProcessError as e:
    print(f'Error: ')
    print(f"Command failed with error: {e.stderr}")

