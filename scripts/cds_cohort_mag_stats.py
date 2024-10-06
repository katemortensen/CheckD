#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 7/1/2024
# Purpose: Compile stats for allelic depth (normalized) across coding domain sequence (CDS)  regions of MAGs in a cohort of MAGs.


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
cds_mag_stats_out = f'{output_dir}/{basename}.cds_mag_stats'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 
# temp vars for testing 

# humanO1

# ad_norm_dir = '/home/kmorten/humanO1_CheckD/ad_norm/'
# output_dir = '/home/kmorten/humanO1_CheckD/ad_norm_cohort_stats/'
# basename = 'humanO1'
# cds_mag_stats_out = f'{output_dir}/{basename}.cds_mag_stats'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)


# %% 
'''Position Stats for CDS only'''

with open(cds_mag_stats_out, 'w') as fp_out:
    # print('CDS Stats')
    header = '# Allelic depth stats for CDS positions.\
        \n# Note that the normalized allelic depth of a reference nucleotide\
        \n# at a position will henceforth be referred to as "ad_norm".\
        \n# If a nucleotide at a position is accounted for multiple times \
        \n# for sequences of the same domain, then the ad_norm value for this \
        \n# position is only accounted for once in order to prevent unwanted \
        \n# dilution of ad_norm values for a domain. \
        \n# (1) MAG name\
        \n# (2) total length of assembly (nucleotides)\
        \n# (3) percent CDS (nucleotides)\
        \n# (4) mean ad_norm for CDS positions\
        \n# (5) median ad_norm for CDS positions\
        \n# (6) standard deviation of ad_norm values for CDS positions\
        \n# (7) percent of CDS ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm.\
        \n# (8) codon1 - mean ad_norm\
        \n# (9) codon1 - median ad_norm\
        \n# (10) codon1 - standard deviation of ad_norm\
        \n# (11) codon1 - percent of ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm\
        \n# (12) codon2 - mean ad_norm\
        \n# (13) codon2 - median ad_norm\
        \n# (14) codon2 - standard deviation of ad_norm\
        \n# (15) codon2 - percent of ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm\
        \n# (16) codon3 - mean ad_norm\
        \n# (17) codon3 - median ad_norm\
        \n# (18) codon3 - standard deviation of ad_norm\
        \n# (19) codon3 - percent of ref positions with < 100% of mapped reads equal to ref nucleotide, aka ad_norm'
    col_len = 19
    fp_out.write(header)
    
    fp_out.write('\n# cols: mag, tot_len, prcnt_cds, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1')
    fp_out.write(', codon1_mean_ad_norm, codon1_median_ad_norm, codon1_stdev_ad_norm, codon1_prcnt_lessthan1')
    fp_out.write(', codon2_mean_ad_norm, codon2_median_ad_norm, codon2_stdev_ad_norm, codon2_prcnt_lessthan1')
    fp_out.write(', codon3_mean_ad_norm, codon3_median_ad_norm, codon3_stdev_ad_norm, codon3_prcnt_lessthan1')
    fp_out.write('\n#')
    fp_out.write('\n# mag\ttot_len\tprcnt_cds\tmean_ad_norm\tmedian_ad_norm\tstdev_ad_norm\tprcnt_lessthan1')
    fp_out.write('\tcodon1_mean_ad_norm\tcodon1_median_ad_norm\tcodon1_stdev_ad_norm\tcodon1_prcnt_lessthan1')
    fp_out.write('\tcodon2_mean_ad_norm\tcodon2_median_ad_norm\tcodon2_stdev_ad_norm\tcodon2_prcnt_lessthan1')
    fp_out.write('\tcodon3_mean_ad_norm\tcodon3_median_ad_norm\tcodon3_stdev_ad_norm\tcodon3_prcnt_lessthan1')

    mag_names = []
    for file in os.listdir(ad_norm_dir):
        if file.endswith('.ad_norm'):
            name = file.strip().split('/')[-1].split('.ad_norm')[0]
            mag_names.append(name)
    mag_names_sorted = sort_list_by_mag(mag_names)

    for mag_name in mag_names_sorted:
        print(mag_name)
        ad_norm_file = f'{ad_norm_dir}/{mag_name}.ad_norm'
        cds_file = f'{ad_norm_dir}/{mag_name}.cds_ad_norm'

        compiled_cds_ad_norms = create_ad_norm_dict(cds_file)
        cds_dict = compiled_cds_ad_norms[0]
        cds_ad_norms = compiled_cds_ad_norms[1]
        compiled_ad_norms = create_ad_norm_dict(ad_norm_file)
        ad_norm_dict = compiled_ad_norms[0]
        all_ad_norms = compiled_ad_norms[1]
        
        tot_len = len(all_ad_norms)
        prcnt_cds = round((len(cds_ad_norms)/tot_len)*100,4)
        general_stats = mean_med_stdev(cds_ad_norms)
        mean_ad_norm = general_stats[0]
        median_ad_norm = general_stats[1]
        stdev_ad_norm = general_stats[2]
        prcnt_lessthan1 = general_stats[3]
        new_line_list = [mag_name,tot_len,prcnt_cds,mean_ad_norm,median_ad_norm,stdev_ad_norm,prcnt_lessthan1]

        '''codon_dict'''
        '''codon > [ad_norms]'''
        
        codon_dict = {'codon1':[],'codon2':[],'codon3':[]}

        for i in cds_dict.keys():
            contig = i.split('_')[0]
            for pos in cds_dict[i].keys():
                if ad_norm_dict[contig].get(pos) != None:
                    ref_nuc_check = ad_norm_dict[contig][pos].get('ref_nuc')
                    ref_nuc = cds_dict[i][pos].get('ref_nuc')
                    assert ref_nuc == ref_nuc_check
                    ad_norm_check = ad_norm_dict[contig][pos].get('ad_norm')
                    ad_norm = cds_dict[i][pos].get('ad_norm')
                    assert ad_norm == ad_norm_check
                    codon_pos = cds_dict[i][pos].get('codon_pos')
                    codon_key = 'codon%s' % (codon_pos)
                    tmp_list = codon_dict.get(codon_key)
                    tmp_list.append(ad_norm)
                    codon_dict[codon_pos] = tmp_list
                else:
                    assert cds_dict[i][pos].get('ref_nuc') == 'NA'
                    assert cds_dict[i][pos].get('ad_norm') == 'NA'
                
        codon1_ad_norms = codon_dict.get('codon1')
        codon2_ad_norms = codon_dict.get('codon2')
        codon3_ad_norms = codon_dict.get('codon3')

        codon1_general_stats = mean_med_stdev(codon1_ad_norms)
        codon1_mean_ad_norm = codon1_general_stats[0]
        codon1_median_ad_norm = codon1_general_stats[1]
        codon1_stdev_ad_norm = codon1_general_stats[2]
        codon1_prcnt_lessthan1 = codon1_general_stats[3]
        new_line_list = new_line_list + [codon1_mean_ad_norm,codon1_median_ad_norm,codon1_stdev_ad_norm,codon1_prcnt_lessthan1]

        codon2_general_stats = mean_med_stdev(codon2_ad_norms)
        codon2_mean_ad_norm = codon2_general_stats[0]
        codon2_median_ad_norm = codon2_general_stats[1]
        codon2_stdev_ad_norm = codon2_general_stats[2]
        codon2_prcnt_lessthan1 = codon2_general_stats[3]
        new_line_list = new_line_list + [codon2_mean_ad_norm,codon2_median_ad_norm,codon2_stdev_ad_norm,codon2_prcnt_lessthan1]
        
        codon3_general_stats = mean_med_stdev(codon3_ad_norms)
        codon3_mean_ad_norm = codon3_general_stats[0]
        codon3_median_ad_norm = codon3_general_stats[1]
        codon3_stdev_ad_norm = codon3_general_stats[2]
        codon3_prcnt_lessthan1 = codon3_general_stats[3]
        new_line_list = new_line_list + [codon3_mean_ad_norm,codon3_median_ad_norm,codon3_stdev_ad_norm,codon3_prcnt_lessthan1]

        assert len(new_line_list) == col_len
        new_line = '\t'.join(str(i) for i in new_line_list)
        fp_out.write(f'\n{new_line}')



# %%
