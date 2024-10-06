#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 6/22/2024
# Purpose: Summarize mapping diversity, contamination, and completion for all mags

# %% 
import argparse
import os
import time 
import csv as csv

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-vcf_file", "--vcf_file", metavar="path-to-vcf-file", help="path to MAG vcf file", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_directory", help="specify the path where the outputs will go", required=True)
parser.add_argument("-basename", "--basename", metavar="basename", help="prefix of output file", required=False, default=False)
# parser.add_argument("-threads", "--threads", metavar="threads", help="threads", required=False, default=False)
# parser.add_argument("-fasta_file", "--fasta_file", metavar="path-to-fasta-file", help="path to MAG fasta file", required=False, default=None)
# parser.add_argument("-qs_score", "--qs_score", help="(boolean flag) quick run but less safe, pulls QS scores from VCF file", action='store_true')


args = parser.parse_args()
vcf_file  = args.vcf_file 
output_dir = args.output_dir
basename = 'out'
if args.basename != None:
    basename = args.basename
output_file = f'{output_dir}/{basename}.ad_norm'

# if fasta_file == None:
#     qs_score == True

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_file = f'{output_dir}/{basename}.ad_norm'
   
# delete this hardcoded comment before running

# %% 
# temp stats for testing

# mag_name = 'bin.339.strict'
# mag_name = 'bin.195.strict'
# mag_name = 'bin.311.strict'
# mag_name = 'bin.213.permissive'
# mag_name = 'bin.76'
# mag_name = 'bin.1'

# ADMB
# vcf_file = f'/home/kmorten/CheckD_old/READS2MAGS/{mag_name}.vcf'
# vcf_file = f'/home/kmorten/ADMB_CheckD/METAWRAP_MAGS/reads2mags_minimap/{mag_name}.vcf'
# fasta_file = f'/home/kmorten/ADMB/METAWRAP_MAGS/BIN_REASSEMBLY/reassembled_bins/{mag_name}.fa'
# output_file = f'/home/kmorten/ADMB_CheckD/METAWRAP_MAGS/allelic_depth/{mag_name}.ad_norm'
# qs_score = False

# humanO1
# vcf_file = f'/home/kmorten/humanO1_CheckD/reads2mags_minimap/{mag_name}.vcf'
# output_dir = f'/home/kmorten/humanO1_CheckD/ad_norm'
# output_file = f'{output_dir}/{mag_name}.ad_norm'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)


# illumina example
# wkdir = f'/home/kmorten/CheckD/illumina_example/illumina_example_checkd'
# vcf_file = f'{wkdir}/reads2mags_minimap/{mag_name}.vcf'
# output_dir = f'/home/kmorten/humanO1_CheckD/ad_norm'
# output_file = f'{output_dir}/{mag_name}.ad_norm'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# %% 
'''

VCF Method

INDELS are skipped, only high quality bases are considered, and a min depth is required.

'''
# start_time = time.time()

qs_score = True
if qs_score == True:
    with open(vcf_file, 'r') as vcf_in, open(f'{output_file}', 'w') as fp_out:
        fp_out_header = '# The below information is separated by ">" for each contig. INDELS identified in the VCF are excluded.\
            \n# (1) nucleotide position in contig\
            \n# (2) reference nucleotide\
            \n# (3) allelic depth, normalized [0,1] (high-quality bases only)\
            \n# (4) depth (high-quality bases only)\
            \n# cols: pos, ref_nuc, ad_norm, depth\
            \n'
        fp_out.write(fp_out_header)
        old_contig_name = None
        old_pos = 0
        for line in vcf_in:
            temp_dict = {}
            if line.startswith('#CHROM') is True:
                header = line.strip('#').strip().split('\t')
            
            if line.startswith('#') is False:
                linep = line.strip().split('\t')
                temp_dict = dict(zip(header,linep))
                AD_key = list(temp_dict.keys())[-1]
                
                # ignore indels 
                if 'INDEL' not in temp_dict.get('INFO'):
                    contig_name = temp_dict.get('CHROM')
                    if contig_name != old_contig_name:
                        fp_out.write(f'\n>{contig_name}')
                        old_pos = 0
                    pos = int(temp_dict.get('POS'))
                    assembly_ref_nuc = temp_dict.get('REF')
                    alt_nucs = temp_dict.get('ALT').split(',')
                    nucs = [assembly_ref_nuc] + alt_nucs
                    freqs = [float(i) for i in temp_dict.get(AD_key).split('\t')[-1].split(':')[-2].split(',')]
                    assert len(freqs) == len(nucs)
                    ndict = dict(zip(nucs,freqs))
                    ndict_sorted = dict(ndict.items(), key=lambda item: item[1])
                    ref_nuc = list(ndict_sorted.keys())[0]
                    ad = int(ndict_sorted.get(ref_nuc))
                    assert ad == ndict.get(ref_nuc)
                    depth = int(temp_dict.get(AD_key).split('\t')[-1].split(':')[1])
                    
                    if depth == 0:
                        ref_nuc = 'N'
                        ad_norm = 0
                        depth = 0
                    else:
                        ad_norm = round(ad/int(depth), 6)
                    
                    
                    if old_pos != pos-1:
                        assert old_pos < pos
                
                    new_line = f'\n{pos}\t{ref_nuc}\t{ad_norm}\t{depth}'
                    fp_out.write(new_line)    

                    
                elif 'INDEL' in temp_dict.get('INFO'):
                    pos = int(temp_dict.get('POS'))
                    assert old_pos == pos
                    
                old_pos = pos 
                old_contig_name = contig_name
                
# end_time = time.time()
# print("total time: ", end_time - start_time)


# %%
