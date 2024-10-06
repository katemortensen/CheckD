#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/2/2024
# Purpose: Compile stats for allelic depth (normalized) across all MAGs in a cohort of MAGs.


# %% 
import argparse
import os
from checkd_functions import create_ad_norm_dict
from checkd_functions import mean_med_stdev

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
domain_stats_out = f'{output_dir}/{basename}.cohort_domain_stats'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 
# temp vars for testing 

# humanO1

# ad_norm_dir = '/home/kmorten/humanO1_CheckD/ad_norm'
# output_dir = '/home/kmorten/humanO1_CheckD/ad_norm_cohort_stats/'
# basename = 'humanO1'
# domain_stats_out = f'{output_dir}/{basename}.cohort_domain_stats'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)


# %% 

'''Position Stats (all positions)'''

with open(domain_stats_out, 'w') as fp_out:
    header = '# Allelic depth for positions by domain type across all assemblies in cohort.\
        \n# Note that the normalized allelic depth of a reference nucleotide\
        \n# at a position will henceforth be referred to as "ad_norm".\
        \n# If a nucleotide at a position is accounted for multiple times \
        \n# for sequences of the same domain, then the ad_norm value for this \
        \n# position is only accounted for once in order to prevent unwanted \
        \n# dilution of ad_norm values for a domain. However, if a nucleotide \
        \n# at a position occurs in more than one domain, then this nucleotide  \
        \n# will be counted once for each domain in which it occurs. Therefore \
        \n# the nucleotides for each domain may sum to more than the total length \
        \n# of the total nucleotides across all MAGs in a cohort.\
        \n# (1) Domain type\
        \n# (2) nucleotides in domain\
        \n# (3) mean ad_norm\
        \n# (4) median ad_norm\
        \n# (5) standard deviation of ad_norm values\
        \n# (6) percent of ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm\
        \n# cols: domain, tot_len, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1\
        \n#\
        \n# domain\ttot_len\tmean_ad_norm\tmedian_ad_norm\tstdev_ad_norm\tprcnt_lessthan1'
    fp_out.write(header)


    extensions = ['.ad_norm','.cds_ad_norm','.ncrna_ad_norm','.dgr_ad_norm','.crispr_ad_norm','.noncds_ad_norm', 'vir_ad_norm']
    domains = ['all','CDS','ncRNA','DGR','CRISPR','non-CDS', 'VIR']
    extDict = dict(zip(extensions,domains))
    for ext in extDict.keys():
        domain = extDict.get(ext)
        ad_norms  = []
        for file in os.listdir(ad_norm_dir):
            if file.endswith(ext):
                ad_norm_file = f'{ad_norm_dir}/{file}'
                new_ad_norms = create_ad_norm_dict(ad_norm_file)[1]
                ad_norms = new_ad_norms + ad_norms
        tot_len = len(ad_norms)
        general_stats = mean_med_stdev(ad_norms)
        mean_ad_norm = general_stats[0]
        median_ad_norm = general_stats[1]
        stdev_ad_norm = general_stats[2]
        prcnt_lessthan1 = general_stats[3]
        new_line = f'\n{domain}\t{tot_len}\t{mean_ad_norm}\t{median_ad_norm}\t{stdev_ad_norm}\t{prcnt_lessthan1}' 
        fp_out.write(new_line)   

                                               
                
# %%
