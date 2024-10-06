#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/2/2024
# Purpose: Viral DNA, exploratory data analysis of mapping diversity stats wtihin a mag

# %% 
import argparse
import os
from checkd_functions import create_desc_stats_dict
import statistics


# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_dir", "--ad_norm_dir", metavar="ad_norm_dir", help="specify the path where input files are (<mag>.vir_ad_norm)", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-basename", "--basename", metavar="basename", help="output basename", required=True)

args = parser.parse_args()
ad_norm_dir = args.ad_norm_dir
output_dir = args.output_dir
basename = args.basename
vir_cohort_desc_stats_out = f'{output_dir}/{basename}.vir_cohort_desc_stats'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# delete this hardcoded comment before running

# %% 
# # Test Vars

# # humanO1
# wkdir = '/home/kmorten/humanO1_CheckD/'
# basename = 'humanO1'

# ################## 

# ad_norm_dir = f'{wkdir}/ad_norm/'
# output_dir = f'{wkdir}/ad_norm_cohort_stats'
# vir_cohort_desc_stats_out = f'{output_dir}/{basename}.vir_cohort_desc_stats'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)


# %% 

with open(vir_cohort_desc_stats_out, 'w') as fp_out:
    vir_header = '# Viral stats for MAG. Normalized allelic depth is "ad_norm".\
        \n# (1) Viral description\
        \n# (2) count, number of viral type occurrences in assembly\
        \n# (3) length, total nucleotides across all occurrences\
        \n# (4) mean ad_norm for viral positions\
        \n# (5) median ad_norm for viral positions\
        \n# (6) standard deviation of ad_norm values for viral positions\
        \n# (7) percent of viral ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm\
        \n# cols: desc, count, length, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1 \
        \n#'
    col_len = 7
    fp_out.write(vir_header)
    fp_out.write('\n# desc\tcount\tlength\tmean_ad_norm\tmedian_ad_norm\tstdev_ad_norm\tprcnt_lessthan1')
    
    mag_names = []
    for file in os.listdir(ad_norm_dir):
        if file.endswith('.ad_norm'):
            name = file.strip().split('/')[-1].split('.ad_norm')[0]
            mag_names.append(name)
            
    vir_ad_norm_files = []
    for mag_name in mag_names:
        vir_ad_norm_file = f'{ad_norm_dir}/{mag_name}.vir_ad_norm'
        vir_ad_norm_files.append(vir_ad_norm_file)

    # virDict = create_desc_stats_dict([vir_ad_norm_files[3], vir_ad_norm_files[13]])
    virDict = create_desc_stats_dict(vir_ad_norm_files)
    
    for desc in virDict.keys():   
        count = virDict[desc].get('count')
        print(f'{desc}: {count}' )
        length = virDict[desc].get('length')
        mean_ad_norm = virDict[desc].get('mean')
        assert mean_ad_norm != None
        median_ad_norm = virDict[desc].get('median')
        stdev_ad_norm = virDict[desc].get('stdev')
        prcnt_lessthan1 = virDict[desc].get('prcnt_lessthan1')          
        
        new_line_list = [desc,count,length,mean_ad_norm,median_ad_norm,stdev_ad_norm,prcnt_lessthan1]

        assert len(new_line_list) == col_len
        new_line = '\t'.join(str(i) for i in new_line_list)
        fp_out.write(f'\n{new_line}')

# %%
