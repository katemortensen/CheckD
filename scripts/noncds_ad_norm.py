#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/14/2024
# Purpose: FragGeneScan CDS region, exploratory data analysis of mapping diversity stats

# %% 
import argparse
import os
from checkd_functions import create_ad_norm_dict

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="path-to-reference-positional-frequency", help="path to <bin>.ad_norm", required=True)
parser.add_argument("-cds_ad_norm_file", "--cds_ad_norm_file", metavar="path-to-reference-positional-frequency", help="path to <bin>.cds_ad_norm", required=True)
parser.add_argument("-fgs_file", "--fgs_file", metavar="path-to-mapping-diversity-stats", help="path to <bin>-fgs.out", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)

args = parser.parse_args()
fgs_file = args.fgs_file
ad_norm_file = args.ad_norm_file
output_dir = args.output_dir
assembly_name = args.assembly_name
mag_name = assembly_name 
cds_ad_norm_file = args.cds_ad_norm_file

noncds_output_file = f'{output_dir}/{mag_name}.noncds_ad_norm'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# delete this hardcoded comment before running

# %% 
# temp vars for testing 

# mag_name = 'bin.195.strict'
# # mag_name = 'bin.311.strict'
# # mag_name = 'bin.213.permissive'
# # mag_name = 'bin.339.strict'
# # mag_name = 'bin.388.strict'
# # mag_name = 'bin.1'
# # mag_name = 'bin.405.strict'
# # test = '>NODE_324_length_541_cov_0.784483'

# # illumina example
# wkdir = f'/home/kmorten/CheckD/illumina_example/illumina_example_checkd/'
# output_dir = f'{wkdir}/ad_norm'
# fgs_file =  f'{wkdir}/fgs_out/{mag_name}-fgs.out'
# ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'
# cds_ad_norm_file = f'{output_dir}/{mag_name}.cds_ad_norm'
# noncds_output_file = f'{output_dir}/{mag_name}.noncds_ad_norm'

# humanO1
# mag_name = 'bin.0'
# output_dir = f'/home/kmorten/humanO1_CheckD/ad_norm'
# fgs_file =  f'/home/kmorten/humanO1_CheckD/fgs_out/{mag_name}-fgs.out'
# ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'
# cds_ad_norm_file = f'{output_dir}/{mag_name}.cds_ad_norm'
# noncds_output_file = f'{output_dir}/{mag_name}.noncds_ad_norm'


# %% 

'''ad_norm_dict'''
'''contig > pos > ref_nuc, ad_norm, etc.'''
def create_ad_norm_dict(ad_norm_file):
    result_dict = {}
    with open(ad_norm_file,'r') as fp_in:
        contig_count = 0
        for line in fp_in:
            if line.startswith('# cols:'):
                info = line.strip().split('# cols: ')[-1].split(', ')
            elif line.startswith('>'):
                contig_name = line.strip()
                contig_count += 1
            elif contig_count >= 1:
                linep = line.strip('\n').split('\t')
                temp_dict = dict(zip(info,linep))
                pos = temp_dict.get('pos')
                temp_dict.pop('pos')
                ad_norm = temp_dict.get('ad_norm')
                if result_dict.get(contig_name) == None:
                    result_dict[contig_name] = {pos:temp_dict}
                else:
                    result_dict[contig_name][pos] = temp_dict
        if contig_count >= 1:
            assert len(result_dict[contig_name].keys()) > 1
            assert len(result_dict.keys()) == contig_count
    result_list = [] # ad_norms list
    if len(result_dict.keys()) > 0:
        if '_' in list(result_dict.keys())[0]: 
            ddupl_dict = {}
            for seq in result_dict.keys():
                contig = seq.split('_')[0]
                for pos in result_dict[seq].keys():
                    ad_norm = result_dict[seq][pos].get('ad_norm')
                    if ddupl_dict.get(contig) == None:
                        ddupl_dict[contig] = {pos:ad_norm}
                        result_list.append(ad_norm)
                    elif ddupl_dict[contig].get(pos) == None:
                        ddupl_dict[contig][pos]=ad_norm
                        result_list.append(ad_norm)
                    if ad_norm != ddupl_dict[contig].get(pos):
                        print(f'{contig}, {pos}')
                        print(f'{seq}')
                        print(f'ad_norm: {ad_norm}')
                        print(f'ddupl_dict: {ddupl_dict[contig].get(pos)}')
                    assert ad_norm == ddupl_dict[contig].get(pos)
        else:
            for contig in result_dict.keys():
                for pos in result_dict[contig].keys():
                    ad_norm = result_dict[contig][pos].get('ad_norm')
                    result_list.append(ad_norm)
    result_list = list(filter(lambda x: x not in ['N', 'NA', 'NaN', None], result_list))
    return result_dict, result_list


# %% 

print(f'{mag_name}')

'''ad_norm_dict'''
'''contig > pos > ref_nuc, ad_norm, etc.'''

ad_norm_dict = create_ad_norm_dict(ad_norm_file)[0]


'noncds_pos_dict'
'''contig > pos > ref_nuc, ad_norm'''

noncds_dict = create_ad_norm_dict(ad_norm_file)[0]

with open(fgs_file,'r') as fgs_in:
        # next(fgs_in)
        old_line = None
        for line in fgs_in:
            if line.startswith('>'):
                contig_name = line.strip()
            else:
                linep = line.strip().split('\t')
                start = int(linep[0])
                end = int(linep[1])
                
                # exclude positions in CDS regions
                for i in range(start, end+1):
                    pos = str(i)
                    if noncds_dict.get(contig_name) != None:
                        if noncds_dict[contig_name].get(pos) != None:
                            noncds_dict[contig_name].pop(pos)
                            assert noncds_dict[contig_name].get(pos) == None

# %% 

'''cds_dict'''
'''contig > pos > ref_nuc, ad_norm, etc.'''

cds_dict = create_ad_norm_dict(cds_ad_norm_file)[0]

for seq in cds_dict.keys():
    contig = seq.split('_')[0]
    if noncds_dict.get(contig)!= None:
        for pos in cds_dict[seq].keys():
            if noncds_dict[contig].get(pos) != None:
                print(f'DUPLICATE: {mag_name}, {contig}, {pos}')

# %%           
                    
'non-CDS positions'

'''contig > pos > frequency of reference at position'''

with open(noncds_output_file, 'w') as noncds_out:
    noncds_out_header = '# Non-coding domain sequences (non-CDS) info.\
        \n# The below information is separated by ">" for each contig. INDELS are excluded.\
        \n# (1) nucleotide position in contig\
        \n# (2) reference nucleotide\
        \n# (3) allelic depth, normalized [0,1] (high-quality bases only)\
        \n# (4) depth (high-quality bases only)\
        \n# cols: pos, ref_nuc, ad_norm, depth\n'
    noncds_out.write(noncds_out_header)
    for contig_name in noncds_dict.keys():
        noncds_out.write(f'\n{contig_name}')
        for pos in noncds_dict[contig_name].keys():
            ref_nuc = noncds_dict[contig_name][pos].get('ref_nuc')
            ad_norm = noncds_dict[contig_name][pos].get('ad_norm')
            depth = noncds_dict[contig_name][pos].get('depth')                                             
            new_line = f'\n{pos}\t{ref_nuc}\t{ad_norm}\t{depth}'
            noncds_out.write(new_line)




