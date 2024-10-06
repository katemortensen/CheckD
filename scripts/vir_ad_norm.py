#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/4/2024
# Purpose: write allelic depth of viral DNA

# %% 
import argparse
import os
import time 
from checkd_functions import create_ad_norm_dict

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)
parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="path to ad_norm file", help="<mag>.ad_norm file", required=True)
parser.add_argument("-final_viral_boundary_file", "--final_viral_boundary_file", metavar="path to virsorter output file (final_viral_boundary_file.tsv file)", help="final_viral_boundary_file.tsv file", required=True)
parser.add_argument("-final_viral_combined_file", "--final_viral_combined_file", metavar="path to virsorter output file (final_viral_combined_file.fa file)", help="final_viral_combined_file.fa file", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)

args = parser.parse_args()
assembly_name = args.assembly_name
ad_norm_file = args.ad_norm_file
final_viral_boundary_file = args.final_viral_boundary_file
final_viral_combined_file = args.final_viral_combined_file
output_dir = args.output_dir

vir_ad_norm_out = f'{output_dir}/{assembly_name}.vir_ad_norm'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 

# humanO1
# assembly_name = 'bin.1'
# wkdir = '/home/kmorten/humanO1_CheckD/'

# # illumina example
# assembly_name = 'bin.1'
# wkdir = '/home/kmorten/CheckD/illumina_example/illumina_example_checkd/'

# output_dir = f'{wkdir}/ad_norm'
# ad_norm_file = f'{output_dir}/{assembly_name}.ad_norm'
# final_viral_boundary_file = f'{wkdir}/virsorter2_out/{assembly_name}/final-viral-boundary.tsv'
# final_viral_combined_file = f'{wkdir}/virsorter2_out/{assembly_name}/final-viral-combined.fa'
# vir_ad_norm_out = f'{output_dir}/{assembly_name}.vir_ad_norm'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)
    
# %% 

with open(vir_ad_norm_out, 'w') as fp_out:
    
    virsorter_run = True
    if os.path.exists(final_viral_combined_file):
        file_size = os.stat(final_viral_combined_file).st_size
        if file_size == 0:
            virsorter_run = False
            print(f'VirSorter2 {final_viral_combined_file} file is empty for {assembly_name}')
            fp_out.write(f'No Viral DNA found by VirSorter2 - {final_viral_combined_file} empty.\n')
            fp_out.flush()
            os.fsync(fp_out.fileno())
            # assert os.stat(vir_ad_norm_out).st_size > 0, "Output file is empty after writing an error message."
    else:
        virsorter_run = False
        print(f'VirSorter2 {final_viral_combined_file} file is empty for {assembly_name}')
        fp_out.write(f'No Viral DNA found by VirSorter2 - {final_viral_combined_file} empty.\n')    
    
    if virsorter_run is True:
        vir_seqs = []
        ad_norm_dict = create_ad_norm_dict(ad_norm_file)[0]

        '''Virus normalized allelic depth'''

        start_time = time.time()

        with open(vir_ad_norm_out, 'w') as fp_out:
            fp_out_header = '# Viral DNA info.\
                \n# The below information is separated by ">" for each contig. INDELS are excluded\
                \n# Contig headers reads as follows: >contig_start_end_stand;description\
                \n# (1) nucleotide position in contig\
                \n# (2) reference nucleotide\
                \n# (3) allelic depth, normalized [0,1] (high-quality bases only)\
                \n# (4) depth (high-quality bases only)\
                \n# cols: pos, ref_nuc, ad_norm, depth\
                \n'
            fp_out.write(fp_out_header)
            
            faDict = {}
            with open(final_viral_combined_file, 'r') as fa_fp:
                count = 0
                for line in fa_fp:
                    linep = line.strip()
                    if line.startswith('>'):
                        seqname = linep.strip('>')
                        count += 1 
                    else:
                        seq = linep
                        faDict[seqname] = seq
                        assert count == len(faDict)
            bdDict = {}
            with open(final_viral_boundary_file, 'r') as bound_fp:
                count = 0
                for line in bound_fp:
                    linep = line.strip()
                    if line.startswith('seqname'):
                        info = linep.split('\t')
                    else:
                        count += 1
                        instance = linep.split('\t') 
                        tmpDict = dict(zip(info, instance))
                        seqname = tmpDict.get('seqname_new')
                        if bdDict.get(seqname) != None:
                            print(seqname)
                        assert bdDict.get(seqname) == None
                        bdDict[seqname] = tmpDict
                        assert len(bdDict) == count
            # NOTE: final-viral-boundary.tsv and final-viral-combined.tsv don't always sum zero
            # assert len(faDict) == len(bdDict)
            for seqname in bdDict.keys():
            # for seqname in faDict:            
                tmp_start = int(bdDict[seqname].get('trim_bp_start'))
                tmp_end = int(bdDict[seqname].get('trim_bp_end'))
                if tmp_start < tmp_end:
                    start = tmp_start
                    end = tmp_end
                    strand = '+'
                    step = 1
                else:
                    start = tmp_end
                    end = tmp_start
                    strand = '-'
                    step = -1
                contig_name = bdDict[seqname].get('seqname')
                what = bdDict[seqname].get('group')
                seq_header = f'{contig_name}_{tmp_start}_{tmp_end}_{strand}'  
                vir_seq = f'{seq_header};{what}'
                if vir_seq not in vir_seqs:
                    vir_seqs.append(vir_seq)
                    fp_out.write(f'\n>{seq_header};{what}')  
                bp_count = 0
                for i in range(start, end + step, step):
                    pos = str(i)
                    # Note: Quality Scores are only given for positions that meet mpileup parameter reqs
                    if ad_norm_dict.get(f'>{contig_name}') != None:
                        if ad_norm_dict[f'>{contig_name}'].get(pos) != None:
                            ref_nuc = ad_norm_dict[f'>{contig_name}'][pos].get('ref_nuc') 
                            if faDict.get(seqname) != None:
                                seq = faDict.get(seqname)
                                fa_nucl = seq[bp_count]                    
                                if fa_nucl != ref_nuc:
                                    print(vir_seq)
                                    print(pos)
                                    print(f'ref_nuc: {ref_nuc}')
                                    print(f'fa_nucl: {fa_nucl}')
                                    assert ref_nuc == 'N'
                            ad_norm = ad_norm_dict[f'>{contig_name}'][pos].get('ad_norm')                    
                            depth = ad_norm_dict[f'>{contig_name}'][pos].get('depth')                            
                    else:
                        ref_nuc = 'N'
                        ad_norm = 'NA'
                        depth = 'NA'
                    new_line = f'\n{pos}\t{ref_nuc}\t{ad_norm}\t{depth}'
                    fp_out.write(new_line)
                    bp_count += 1
                    
        assert len(vir_seqs) == len(set(vir_seqs))               
                        
                        
                        
            
  
# %%
