#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/14/2024
# Purpose: rename contig names in assembly fasta file to avoid parsing issues 

'''NOTE: mag names cannot have underscores - compromises parsing scripts.'''
'''NOTE: contig names cannot have underscores - compromises parsing scripts.'''

# %% 
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("-mag_dir", "--mag_dir", metavar="mag_dir", help="specify the path what directory has the mags", required=True)
# parser.add_argument("-renamed_mag_dir", "--renamed_mag_dir", metavar="rename_mag_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)

args = parser.parse_args()
mag_dir = args.mag_dir
# renamed_mag_dir = args.renamed_mag_dir
output_dir  = args.output_dir
renamed_mag_dir = f'{output_dir}/renamed_mags'
naming_key_path = f'{output_dir}/naming_key.tsv'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(renamed_mag_dir):
    os.makedirs(renamed_mag_dir)

# delete this hardcoded comment before running

# %% 

# illumina example

# wkdir =  f'/home/kmorten/CheckD/illumina_example/illumina_example_checkd'
# output_dir =  f'/home/kmorten/CheckD/illumina_example/illumina_example_checkd'
# mag_dir = f'/home/kmorten/CheckD/illumina_example/mags'
# renamed_mag_dir = f'{output_dir}/renamed_mags'
# naming_key_path = f'{output_dir}/naming_key.tsv'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# if not os.path.exists(renamed_mag_dir):
#     os.makedirs(renamed_mag_dir)
    
# %% 

    
mags = os.listdir(mag_dir)
with open (naming_key_path, 'w') as fp_key_out:
    header = f'# Note: Assembly and contigs are renamed to avoid interference with parsing scripts.\
        \n# Key Structure:\
        \n# old_assembly_name\tnew_assembly_name\
        \n# >old_contig_name\t>new_contig_name\
        \n#'
    fp_key_out.write(header)
    mag_count = 0
    for m in mags:            
        extension = m.strip().split('.')[-1]
        old_mag_name = m.strip(f'.{extension}')
        new_mag_name = f'bin.{mag_count}'
        fp_key_out.write(f'\n{old_mag_name}\t{new_mag_name}')
        with open(f'{renamed_mag_dir}/{new_mag_name}.fa', 'w') as fp_mag_out:
            pmag = f'{mag_dir}/{m}'
            with open(pmag, 'r') as fp_in:
                contig_count = 1
                for line in fp_in:
                    if line.startswith('>'):
                        old_contig_name = line.strip().strip('>')
                        new_contig_name = f'contig.{contig_count}'
                        fp_key_out.write(f'\n>{old_contig_name}\t>{new_contig_name}')
                        fp_mag_out.write(f'>{new_contig_name}\n')
                        contig_count += 1
                    else:
                        linep = line.strip()
                        fp_mag_out.write(f'{linep}\n')
        mag_count += 1

# %%
