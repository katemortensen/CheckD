#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/17/2024
# Purpose: write allelic depth of CRISPR artifacts to file

# %% 
import argparse
import os
import time 
from checkd_functions import create_ad_norm_dict

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="path-to-reference-positional-frequency", help="path to <bin>.ad_norm", required=True)
parser.add_argument("-crt_file", "--crt_file", metavar="path to CRISPRone crt file", help="if file not present, run CRISPRone", required=True)
parser.add_argument("-crisprone_gff_file", "--crisprone_gff_file", metavar="path to CRISPRone crt file", help="if file not present, run CRISPRone", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)

args = parser.parse_args()
ad_norm_file = args.ad_norm_file
crt_file = args.crt_file
crisprone_gff_file = args.crisprone_gff_file
output_dir = args.output_dir
mag_name = args.assembly_name

crispr_ad_norm_out = f'{output_dir}/{mag_name}.crispr_ad_norm'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 
# # illumina example 
# mag_name = 'bin.4'
# wkdir = '/home/kmorten/CheckD/illumina_example/illumina_example_checkd'
# output_dir = f'{wkdir}/ad_norm'
# crt_file = f'{wkdir}/crisprone_out/{mag_name}/{mag_name}.crt'
# crisprone_gff_file = f'{wkdir}/crisprone_out/{mag_name}/{mag_name}-sm.gff'
# ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'
# crispr_ad_norm_out = f'{output_dir}/{mag_name}.crispr_ad_norm'


# humanO1
# mag_name = 'bin.6'
# output_dir = '/home/kmorten/humanO1_CheckD/ad_norm'
# crt_file = f'/home/kmorten/humanO1_CheckD/crisprone_out/{mag_name}/{mag_name}.crt'
# crisprone_gff_file = f'/home/kmorten/humanO1_CheckD/crisprone_out/{mag_name}/{mag_name}-sm.gff'
# ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'
# crispr_ad_norm_out = f'{output_dir}/{mag_name}.crispr_ad_norm'


# %% 

with open(crispr_ad_norm_out, 'w') as fp_out:
    try:
        file_size = os.stat(crisprone_gff_file).st_size
        if file_size == 0:
            print(f'CRISPRone gff file is empty for {mag_name}')
            fp_out.write("No CRISPR artifacts found - CRISPRone gff file empty.\n")
            fp_out.flush()
            os.fsync(fp_out.fileno())
            # assert os.stat(crispr_ad_norm_out).st_size > 0, "Output file is empty after writing an error message."
    except FileNotFoundError as e:
        print(f'CRISPRone gff file NOT found for {mag_name}')
        fp_out.write("No CRISPR artifacts found - CRISPRone gff file not found.\n")
        fp_out.flush()
        os.fsync(fp_out.fileno())
        # assert os.stat(crispr_ad_norm_out).st_size > 0, "Output file is empty after writing a FileNotFound error message."
    else:
        
        ad_norm_dict = create_ad_norm_dict(ad_norm_file)[0]

        def find_strand(start, end):
            if start < end:
                strand = '+'
            elif start > end:
                strand = '-'
            return strand

        def query_ad_norm_dict(contig, start, end, strand, seq_header, crispr_seqs):
            if seq_header not in crispr_seqs:
                fp_out.write(f'\n{seq_header}')
                crispr_seqs.append(seq_header)
                if contig.startswith('>') == False:
                    contig_name = f'>{contig}'
                if strand == '-':
                    step = -1
                else:
                    strand = "+"
                    step = 1                    
                for i in range(start, end, step):
                    pos = str(i)
                    if ad_norm_dict.get(contig_name) != None:
                        if ad_norm_dict[contig_name].get(pos) != None:
                            ref_nuc = ad_norm_dict[contig_name][pos].get('ref_nuc')                            
                            ad_norm = ad_norm_dict[contig_name][pos].get('ad_norm')                    
                            depth = ad_norm_dict[contig_name][pos].get('depth')                            
                        else:
                            ref_nuc = 'N'
                            ad_norm = 'NA'
                            depth = 'NA'                       
                    else:
                        ref_nuc = 'N'
                        ad_norm = 'NA'
                        depth = 'NA'
                    new_line = f'\n{pos}\t{ref_nuc}\t{ad_norm}\t{depth}'
                    fp_out.write(new_line)
                
        
        '''Parse CRISPRone crt file'''        
        fp_out_header = '# CRISPR artifact info.\
            \n# The below information is separated by ">" for each contig. INDELS are excluded\
            \n# Contig headers reads as follows: >contig_start_end_stand;description\
            \n# (1) nucleotide position in contig\
            \n# (2) reference nucleotide\
            \n# (3) allelic depth, normalized [0,1] (high-quality bases only)\
            \n# (4) depth (high-quality bases only)\
            \n# cols: pos, ref_nuc, ad_norm, depth\
            \n'
        fp_out.write(fp_out_header)
        
        crispr_seqs = []
        with open(crt_file, 'r') as crt_in:
            switch = 0
            for line in crt_in:
                if line.startswith('---') and switch == 3:
                    switch = 0
                if switch == 3:
                    linep = line.strip().split('\t')
                    if line.strip().endswith(']') == True:            
                        repeat_length = int(line.split(',')[0].split(' ')[1])
                        spacer_length = int(line.split(',')[1].split(' ')[1])
                        # repeat_start = int(line.split('\t')[0])
                        repeat_start = int(line.split('\t')[0].strip('-'))
                        repeat_end = repeat_start + repeat_length - 1 
                        spacer_start = repeat_start + repeat_length
                        spacer_end = spacer_start + spacer_length - 1
                        strand = find_strand(repeat_start, repeat_end)
                        
                        seq_header = f'>{contig}_{repeat_start}_{repeat_end}_{strand};repeat'
                        query_ad_norm_dict(contig, repeat_start, repeat_end, strand, seq_header, crispr_seqs)
                        
                        seq_header = f'>{contig}_{spacer_start}_{spacer_end}_{strand};spacer'
                        query_ad_norm_dict(contig, spacer_start, spacer_end, strand, seq_header, crispr_seqs)
                        
                    else:
                        repeat_length = len(linep[-1])
                        # repeat_start = int(linep[0])
                        repeat_start = int(line.split('\t')[0].strip('-'))                    
                        repeat_end = repeat_start + repeat_length - 1 
                        strand = find_strand(repeat_start, repeat_end)

                        seq_header = f'>{contig}_{repeat_start}_{repeat_end}_{strand};repeat'
                        query_ad_norm_dict(contig, repeat_start, repeat_end, strand, seq_header, crispr_seqs)
                if line.startswith('SEQ'):
                    contig = line.strip().split(' ')[1]
                # if line.startswith('CRISPR'):
                    # linep = line.strip().split(' ')
                    # start = int(linep[3])
                    # end = int(linep[5])
                    # strand = find_strand(start, end)
                if line.startswith('POSITION'):
                    linep = line.split('\t')
                    switch += 1
                    if linep[1] == 'REPEAT':
                        switch += 1
                if line.startswith('---') and switch == 2 :
                    switch += 1
            
        '''Parse CRISPRone gff file'''
        with open(crisprone_gff_file, 'r') as gff_in:
            for line in gff_in:
                if 'what' in line:
                    linep = line.strip().split('\t')
                    contig = linep[0]
                    what = linep[-1].split('what=')[1].split(';')[0]
                    if what != 'CRISPR':
                        start = int(linep[3])
                        end = int(linep[4])
                        strand = find_strand(start, end)
                        assert start < end 
                        seq_header = f'>{contig}_{start}_{end}_{strand};{what}'
                        query_ad_norm_dict(contig, start, end, strand, seq_header, crispr_seqs)
        assert len(crispr_seqs) == len(set(crispr_seqs))

assert os.stat(crispr_ad_norm_out).st_size > 0

            
            
        

# %%
