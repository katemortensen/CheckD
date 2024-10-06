#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/2/2024
# Purpose: Compile cohort description stats across all regions types 

# %% 
import argparse
import os

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_cohort_stats_dir", "--ad_norm_cohort_stats_dir", metavar="ad_norm_cohort_stats_dir", help="specify the path where input files are (<basename>.<type>_cohort_desc_stats)", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-basename", "--basename", metavar="basename", help="output basename", required=True)

args = parser.parse_args()
ad_norm_cohort_stats_dir = args.ad_norm_cohort_stats_dir
output_dir = args.output_dir
basename = args.basename
compiled_cohort_desc_stats_out = f'{output_dir}/{basename}.cohort_desc_stats'

types = ['ncrna', 'crispr', 'dgr', 'vir']

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# delete this hardcoded comment before running

# %% 
# Test Vars

# # humanO1
# wkdir = '/home/kmorten/humanO1_CheckD/'
# basename = 'humanO1'

# ################### 

# ad_norm_cohort_stats_dir = f'{wkdir}/ad_norm_cohort_stats'
# output_dir = f'{wkdir}/ad_norm_cohort_stats'
# compiled_cohort_desc_stats_out = f'{output_dir}/{basename}.cohort_desc_stats'
# types = ['ncrna', 'crispr', 'dgr', 'vir]

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# %% 

compile_desc_dict = {}
for t in types:
    file = f'{ad_norm_cohort_stats_dir}/{basename}.{t}_cohort_desc_stats' 
    with open(file, 'r') as fp_in:
        for line in fp_in:
            if line.startswith('# cols:'):
                info = line.strip('# cols: ').strip('\n').strip(' ').split(', ')
            elif line.startswith('#') is False:
                linep = line.strip().split('\t')
                tmp_dict = dict(zip(info,linep))
                desc = tmp_dict.get('desc')
                assert compile_desc_dict.get(desc) == None
                tmp_dict.pop('desc')
                compile_desc_dict[desc] = tmp_dict
                compile_desc_dict[desc]['type'] = t

compile_desc_dict = dict(sorted(compile_desc_dict.items(), key=lambda item: item[1]['mean_ad_norm'], reverse=False))

with open(compiled_cohort_desc_stats_out,'w') as fp_out:
    fp_out.write('# Sequence and region stats for cohort of MAGs. Normalized allelic depth is "ad_norm".\
        \n# (1) description\
        \n# (2) region type\
        \n# (3) count, number of occurrences across all MAGs in cohort\
        \n# (4) length, total nucleotides across all occurrences\
        \n# (5) mean ad_norm for positions\
        \n# (6) median ad_norm for positions\
        \n# (7) standard deviation of ad_norm values for positions\
        \n# (8) percent of ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm\
        \n# cols: desc, type, count, length, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1\
        \n#\
        \n# desc\ttype\tcount\tlength\tmean_ad_norm\tmedian_ad_norm\tstdev_ad_norm\tprcnt_lessthan1')
    col_len = 8
    for desc in compile_desc_dict.keys():
        t = compile_desc_dict[desc].get('type')
        count = compile_desc_dict[desc].get('count')
        length = compile_desc_dict[desc].get('length')
        mean_ad_norm = compile_desc_dict[desc].get('mean_ad_norm')
        median_ad_norm = compile_desc_dict[desc].get('median_ad_norm')
        stdev_ad_norm = compile_desc_dict[desc].get('stdev_ad_norm')
        prcnt_lessthan1 = compile_desc_dict[desc].get('prcnt_lessthan1')          
        
        new_line_list = [desc,t,count,length,mean_ad_norm,median_ad_norm,stdev_ad_norm,prcnt_lessthan1]

        assert len(new_line_list) == col_len
        new_line = '\t'.join(str(i) for i in new_line_list)
        fp_out.write(f'\n{new_line}')     
    
# %%
