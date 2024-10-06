#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 6/24/2024
# Purpose: write allelic depth of non-coding RNA (ncRNAs) to file

# %% 
import argparse
import os

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="path-to-reference-positional-frequency", help="path to <bin>.ad_norm", required=True)
parser.add_argument("-tblout_deovlp_file", "--tblout_deovlp_file", metavar="path to Rfam cmscan tblout_deovlp_file (deoverlapped)", help="if file not present, run rfam.py 1st", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)

args = parser.parse_args()
ad_norm_file= args.ad_norm_file
tblout_deovlp_file = args.tblout_deovlp_file
output_dir = args.output_dir
mag_name = args.assembly_name

ncrna_ad_norm_out = f'{output_dir}/{mag_name}.ncrna_ad_norm'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 
# # Test Vars
# mag_name = 'bin.76'
# output_dir = '/home/kmorten/humanO1_CheckD/ad_norm'
# tblout_deovlp_file = f'/home/kmorten/humanO1_CheckD/rfam_out/{mag_name}.deoverlapped.tblout'
# ad_norm_file= f'{output_dir}/{mag_name}.ad_norm'
# ncrna_ad_norm_out = f'{output_dir}/{mag_name}.ncrna_ad_norm'


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
                if result_dict.get(contig_name) == None:
                    result_dict[contig_name] = {pos:temp_dict}
                else:
                    result_dict[contig_name][pos] = temp_dict        
        assert len(result_dict[contig_name].keys()) > 1
        assert len(result_dict.keys()) == contig_count 
    return result_dict

ad_norm_dict = create_ad_norm_dict(ad_norm_file)
    
'''ncRNA freq(pos)'''
ncrna_seqs = []
with open(ncrna_ad_norm_out, 'w') as fp_out:
    fp_out_header = '# Non-coding RNA info.\
        \n# The below information is separated by ">" for each contig. INDELS are excluded\
        \n# Contig headers reads as follows: >contig_start_end_stand;description\
        \n# (1) nucleotide position in contig\
        \n# (2) reference nucleotide\
        \n# (3) allelic depth, normalized [0,1] (high-quality bases only)\
        \n# (4) depth (high-quality bases only)\
        \n# cols: pos, ref_nuc, ad_norm, depth\
        \n'
    fp_out.write(fp_out_header)
    with open(tblout_deovlp_file) as fp_in:
        for line in fp_in:
            if line.startswith('#') == False:
                linep = line.strip().split(' ')
                entry = list(filter(None, linep))
                contig_name = f'>{entry[3]}'
                start = int(entry[9])
                end = int(entry[10])
                strand = entry[11]  
                desc = entry[1]
                # desc = line.strip()[216:]
                seq_header = f'{contig_name}_{start}_{end}_{strand};{desc}'
                ncrna_seqs.append(seq_header)
                fp_out.write(f'\n{seq_header}') 
                if strand == '+':
                    step = 1               
                if strand == '-':
                    step = -1      
                count = 0
                for i in range(start, end + step, step):
                    pos = str(i)
                    '''Read Count Acceptability'''
                    # Note: Quality Scores are only given for positions that meet mpileup parameter reqs
                    if ad_norm_dict.get(contig_name) != None:
                        if ad_norm_dict[contig_name].get(pos) != None:
                            ref_nuc = ad_norm_dict[contig_name][pos].get('ref_nuc')                            
                            ad_norm = ad_norm_dict[contig_name][pos].get('ad_norm')                    
                            depth = ad_norm_dict[contig_name][pos].get('depth')                            
                    else:
                        ref_nuc = 'N'
                        ad_norm = 'NA'
                        depth = 'NA'
                    new_line = f'\n{pos}\t{ref_nuc}\t{ad_norm}\t{depth}'
                    fp_out.write(new_line)
                    count += 1
                    assert count > 0

if len(ncrna_seqs) != len(set(ncrna_seqs)):
    print('Issue: Sequence recorded in duplicate')

assert len(ncrna_seqs) == len(set(ncrna_seqs))
assert os.stat(ncrna_ad_norm_out).st_size > 0



# %%
