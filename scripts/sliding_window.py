#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 8/6/2024
# Purpose: EDA cohort_mag_desc_stats output


# %% 
import argparse
import os
import pandas as pd 
from checkd_functions import create_ad_norm_dict
from checkd_functions import create_general_stats_dict

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="path-to-reference-positional-frequency", help="path to <bin>.ad_norm", required=True)
parser.add_argument("-window", "--window", metavar="window size", help="window size", required=True)
parser.add_argument("-step", "--step", metavar="step size", help="step size", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)

args = parser.parse_args()
window = int(args.window)
step = int(args.step)
ad_norm_file = args.ad_norm_file
output_dir = args.output_dir
assembly_name = args.assembly_name

output_file = f'{output_dir}/{assembly_name}.w{window}_s{step}'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# delete this hardcoded comment before running

# %% 
# Test vars

# humanO1
# wkdir = '/home/kmorten/humanO1_CheckD/'
# assembly_name = 'bin.13'
# ad_norm_file =  f'{wkdir}/ad_norm/{assembly_name}.ad_norm'
# window = 50000
# step = 10000
# output_dir = f'{wkdir}/sliding_window/'
# output_file = f'{output_dir}/{assembly_name}.w{window}_s{step}'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)


# %% 

ad_norm_Dict = create_ad_norm_dict(ad_norm_file)[0]

with open(output_file, 'w') as fp_out:

    header = f'# General stats for window size {window} and step size {step}.\
        \n# Note that the normalized allelic depth of a reference nucleotide\
        \n# at a position will henceforth be reffered to as "ad_norm".\
        \n# (1) contig\
        \n# (2) mean ad_norm for positions in window\
        \n# (3) median ad_norm for positions in window\
        \n# (4) standard deviation of ad_norm values for positions in window\
        \n# (5) percent of ref positions in window with < 100% of mapped reads equal to ref nucleotide, aka ad_norm\
        \n# cols: contig, pos, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1\
        \n#\
        \n# contig\tpos\tmean_ad_norm\tmedian_ad_norm\tstdev_ad_norm\tprcnt_lessthan1'
        
    fp_out.write(header)
    
    for contig in ad_norm_Dict.keys():     
        first = int(list(ad_norm_Dict[contig].keys())[0])
        last = int(list(ad_norm_Dict[contig].keys())[-1])
        start = first
        count = 0
        while start <= last:
            end = start + window -1
            ad_norm_list = []
            for pos in range(start,end,1):
                pos = str(pos)
                if ad_norm_Dict[contig].get(pos) != None:
                    new_ad_norm = float(ad_norm_Dict[contig][pos].get('ad_norm'))
                    ad_norm_list.append(new_ad_norm)
            general_stats_dict = create_general_stats_dict(ad_norm_list)
            assert int(general_stats_dict.get('length')) <= window
            mean_ad_norm = general_stats_dict.get('mean')
            median_ad_norm = general_stats_dict.get('median')
            stdev_ad_norm = general_stats_dict.get('stdev')
            prcnt_lessthan1 = general_stats_dict.get('prcnt_lessthan1')
            contig_name = contig.strip('>')
            new_line = f'\n{contig_name}\t{pos}\t{mean_ad_norm}\t{median_ad_norm}\t{stdev_ad_norm}\t{prcnt_lessthan1}'
            fp_out.write(new_line)  
            start += step
            count += 1
            # print(new_line)
        assert int(pos) >= last
        assert count >= (last-first)/step
        
            

# %%
