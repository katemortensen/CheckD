#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/16/2024
# Purpose: nucleotide accounting: CDS + nonCDS = all nucleotides in assembly

# %% 
import argparse
import os
from checkd_functions import create_ad_norm_dict
from checkd_functions import sort_list_by_mag

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_dir", "--ad_norm_dir", metavar="path-to-ad_norm-files", help="path to <bin>.ad_norm", required=True)

args = parser.parse_args()
ad_norm_dir = args.ad_norm_dir

# delete this hardcoded comment before running

# %% 
# temp vars for testing 

# illumina testing
# ad_norm_dir = f'/home/kmorten/CheckD/illumina_example/illumina_example_checkd/ad_norm'

# # humanO1
# ad_norm_dir = '/home/kmorten/humanO1_CheckD/ad_norm'

            
# %% 
all_count = 0
cds_count = 0
noncds_count = 0

mag_names = []
for file in os.listdir(ad_norm_dir):
    if file.endswith('.ad_norm'):
        name = file.strip().split('/')[-1].split('.ad_norm')[0]
        mag_names.append(name)
mag_names_sorted = sort_list_by_mag(mag_names)

for mag_name in mag_names_sorted:
    print(mag_name)
    ad_norm_file = f'{ad_norm_dir}/{mag_name}.ad_norm'
    cds_ad_norm_file = f'{ad_norm_dir}/{mag_name}.cds_ad_norm'
    noncds_ad_norm_file = f'{ad_norm_dir}/{mag_name}.noncds_ad_norm'

    all_count_tmp = len(create_ad_norm_dict(ad_norm_file)[1])
    cds_count_tmp = len(create_ad_norm_dict(cds_ad_norm_file)[1])
    noncds_count_tmp = len(create_ad_norm_dict(noncds_ad_norm_file)[1])

    if all_count_tmp != cds_count_tmp + noncds_count_tmp:
        print(f'all nucls: {all_count_tmp}')
        print(f'CDS nucls: {cds_count_tmp}')
        print(f'nonCDS nucls: {noncds_count_tmp}')
        
        all_dict = create_ad_norm_dict(ad_norm_file)[0]
        cds_dict = create_ad_norm_dict(cds_ad_norm_file)[0]
        noncds_dict = create_ad_norm_dict(noncds_ad_norm_file)[0]

        for seq in cds_dict.keys():
            contig = seq.split('_')[0]
            for pos in cds_dict[seq].keys():
                if all_dict[contig].get(pos) !=  None:
                    all_dict[contig].pop(pos)
                else:
                    print(mag_name)
                    print(seq)
                    print(contig)
                    print(pos)
        
    all_count += all_count_tmp
    cds_count += cds_count_tmp
    noncds_count += noncds_count_tmp
    
    
print(f'Cohort')   
print(f'all nucls: {all_count}')
print(f'CDS nucls: {cds_count}')
print(f'nonCDS nucls: {noncds_count}')
# assert all_count == cds_count + noncds_count 

if all_count == cds_count + noncds_count :
        print(f'Warning: Discrepency between FragGeneScan and mpileup (VCF) INDEL identification.')
        print(f'Likely FragGeneScan identified more INDELS than mpileup.')

# %%
