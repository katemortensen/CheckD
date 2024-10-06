#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 8/17/2024
# Purpose: Compile stats for allelic depth (normalized) across CRISPR regions and artifacts of MAGs in a cohort of MAGs.


# %% 
import argparse
import os
from checkd_functions import create_ad_norm_dict
from checkd_functions import mean_med_stdev
from checkd_functions import sort_list_by_mag

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_dir", "--ad_norm_dir", metavar="ad_norm_dir", help="specify the path where input files are (<mag>-fgs.ffn, <mag>-fgs.out, <mag>.ref_pos_req)", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-basename", "--basename", metavar="basename", help="output basename", required=True)

args = parser.parse_args()
ad_norm_dir = args.ad_norm_dir
output_dir = args.output_dir
basename = args.basename
mag_name = basename 
crispr_mag_stats_out = f'{output_dir}/{basename}.crispr_mag_stats'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 
# temp vars for testing 

# illumina example
# wkdir = '/home/kmorten/CheckD/illumina_example/illumina_example_checkd/'
# ad_norm_dir = f'{wkdir}/ad_norm'
# output_dir = '/home/kmorten/humanO1_CheckD/ad_norm_cohort_stats/'
# basename = 'humanO1'
# crispr_mag_stats_out = f'{output_dir}/{basename}.crispr_mag_stats'



# humanO1

# ad_norm_dir = '/home/kmorten/humanO1_CheckD/ad_norm/'
# output_dir = '/home/kmorten/humanO1_CheckD/ad_norm_cohort_stats/'
# basename = 'humanO1'
# crispr_mag_stats_out = f'{output_dir}/{basename}.crispr_mag_stats'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)


# %%
'''Position Stats for CRISPR regions'''

with open(crispr_mag_stats_out, 'w') as fp_out:
    # print('CRISPR Stats.')
    header = '# Allelic depth for any position within a CRISPR artifact or region.\
        \n# Note that the normalized allelic depth of a reference nucleotide\
        \n# at a position will henceforth be referred to as "ad_norm".\
        \n# If a nucleotide at a position is accounted for multiple times \
        \n# for sequences of the same domain, then the ad_norm value for this \
        \n# position is only accounted for once in order to prevent unwanted \
        \n# dilution of ad_norm values for a domain. \
        \n# (1) MAG name\
        \n# (2) total length of assembly (nucleotides)\
        \n# (3) percent CRISPR (nucleotides)\
        \n# (4) CRISPR artifact instances (count)\
        \n# (5) mean ad_norm for CRISPR positions\
        \n# (6) median ad_norm for CRISPR positions\
        \n# (7) standard deviation of ad_norm values for CRISPR positions\
        \n# (8) percent of CRISPR ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm'
    col_len = 8
    fp_out.write(header)
    fp_out.write('\n# cols: mag, tot_len, prcnt_crispr, crispr_count, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1')
    fp_out.write('\n#')
    fp_out.write('\n# mag\ttot_len\tprcnt_crispr\tcrispr_count\tmean_ad_norm\tmedian_ad_norm\tstdev_ad_norm\tprcnt_lessthan1')

    mag_names = []
    for file in os.listdir(ad_norm_dir):
        if file.endswith('.ad_norm'):
            name = file.strip().split('/')[-1].split('.ad_norm')[0]
            mag_names.append(name)
    mag_names_sorted = sort_list_by_mag(mag_names)

    for mag_name in mag_names_sorted:
        print(mag_name)
        ad_norm_file = f'{ad_norm_dir}/{mag_name}.ad_norm'
        crispr_file = f'{ad_norm_dir}/{mag_name}.crispr_ad_norm'

        crispr_count = 0
        with open(crispr_file, 'r') as fp_in:
            for line in fp_in:
                if line.startswith('>'):
                    crispr_count += 1

        compiled_crispr_ad_norms = create_ad_norm_dict(crispr_file)
        crispr_dict = compiled_crispr_ad_norms[0]
        crispr_ad_norms = compiled_crispr_ad_norms[1]
        compiled_ad_norms = create_ad_norm_dict(ad_norm_file)
        ad_norm_dict = compiled_ad_norms[0]
        all_ad_norms = compiled_ad_norms[1]
        
        tot_len = len(all_ad_norms)
        prcnt_crispr = round((len(crispr_ad_norms)/tot_len)*100,4)
        general_stats = mean_med_stdev(crispr_ad_norms)
        mean_ad_norm = general_stats[0]
        median_ad_norm = general_stats[1]
        stdev_ad_norm = general_stats[2]
        prcnt_lessthan1 = general_stats[3]
        new_line_list = [mag_name,tot_len,prcnt_crispr,crispr_count,mean_ad_norm,median_ad_norm,stdev_ad_norm,prcnt_lessthan1]
        
        assert len(new_line_list) == col_len
        new_line = '\t'.join(str(i) for i in new_line_list)
        fp_out.write(f'\n{new_line}')

# %%
