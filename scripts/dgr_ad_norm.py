#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 7/31/2024
# Purpose: write allelic depth of diversity generating retroelements (DGRs) to file 

# %% 
import argparse
import os
import time 
from checkd_functions import create_ad_norm_dict

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="path-to-reference-positional-frequency", help="path to <bin>.ad_norm", required=True)
parser.add_argument("-dgr_gff_file", "--dgr_gff_file", metavar="path to DGRscan GFF file", help="if file not present, run DGRscan first", required=True)
# parser.add_argument("-dgr_summary_file", "--dgr_summary_file", metavar="path to DGRscan summary file", help="if file not present, DGRs are not found in assembly", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)

args = parser.parse_args()
ad_norm_file = args.ad_norm_file
# dgr_summary_file = args.dgr_summary_file
dgr_gff_file = args.dgr_gff_file
output_dir = args.output_dir
mag_name = args.assembly_name

dgr_ad_norm_out = f'{output_dir}/{mag_name}.dgr_ad_norm'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 

# humanO1

# # mag_name = 'bin.41'
# mag_name = 'bin.62'
# # mag_name = 'bin.75'

# wkdir = '/home/kmorten/humanO1_CheckD/'
# output_dir = f'{wkdir}/ad_norm'
# dgr_summary_file = f'{wkdir}/dgrscan_out/{mag_name}/{mag_name}-DGRscan.summary'
# dgr_gff_file = f'{wkdir}/dgrscan_out/{mag_name}/{mag_name}-DGRscan.gff'
# ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'
# dgr_ad_norm_out = f'{output_dir}/{mag_name}.dgr_ad_norm'

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# %% 

with open(dgr_ad_norm_out, 'w') as fp_out:
    try:
        file_size = os.stat(dgr_gff_file).st_size
        if file_size == 0:
            print(f'DGRscan gff file is empty for {mag_name}')
            fp_out.write("No DGR systems found by DGRscan - DGRscan gff file empty.\n")
            fp_out.flush()
            os.fsync(fp_out.fileno())
            # assert os.stat(dgr_ad_norm_out).st_size > 0, "Output file is empty after writing an error message."
    except FileNotFoundError as e:
        print(f'DGRscan gff file NOT found for {mag_name}')
        fp_out.write("No DGR systems found - DGRscan gff file not found.\n")
        fp_out.flush()
        os.fsync(fp_out.fileno())
        # assert os.stat(dgr_ad_norm_out).st_size > 0, "Output file is empty after writing a FileNotFound error message."
    else:
        
        ad_norm_dict = create_ad_norm_dict(ad_norm_file)[0]
        
        '''DGR freq(pos)'''

        start_time = time.time()

        dgr_seqs = []
        with open(dgr_ad_norm_out, 'w') as fp_out:
            fp_out_header = '# Diversity Generating Retroelement (DGR) info.\
                \n# The below information is separated by ">" for each contig. INDELS are excluded\
                \n# Contig headers reads as follows: >contig_start_end_stand;description\
                \n# (1) nucleotide position in contig\
                \n# (2) reference nucleotide\
                \n# (3) allelic depth, normalized [0,1] (high-quality bases only)\
                \n# (4) depth (high-quality bases only)\
                \n# cols: pos, ref_nuc, ad_norm, depth\
                \n'
            fp_out.write(fp_out_header)
            with open(dgr_gff_file) as fp_in:
                next(fp_in)
                next(fp_in)
                for line in fp_in:
                    if 'DGRscan' in line:
                        linep = line.strip().split('\t')
                        seq_header = linep[-1].strip('ID=').split(';')[0]
                        seq_header = f'>{seq_header}'
                        what = line.split('what=')[1].split(';')[0].strip()
                        seq_headerp = seq_header.split('_')
                        contig_name = seq_headerp[0]
                        if len(seq_headerp)!=4:
                            tmp_start = linep[3]
                            tmp_end = linep[4]
                            strand = linep[6]                        
                            seq_header = f'{contig_name}_{tmp_start}_{tmp_end}_{strand}'
                            seq_headerp = seq_header.split('_')
                        strand = seq_header[-1]
                        dgr_seq = f'{seq_header};{what}'
                        if dgr_seq not in dgr_seqs:
                            dgr_seqs.append(dgr_seq)
                            fp_out.write(f'\n{seq_header};{what}')      
                            if strand == '-':
                                step = -1
                                end = int(linep[3])
                                start = int(linep[4])
                                if '_' in seq_header:
                                    assert end == int(seq_headerp[1])
                                    assert start == int(seq_headerp[2])
                            else:
                                strand = "+"
                                step = 1  
                                start = int(linep[3])
                                end = int(linep[4])
                                if '_' in seq_header:
                                    start == int(seq_headerp[1])
                                    end == int(seq_headerp[2])    
                            for i in range(start, end + step, step):
                                pos = str(i)
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
                                
                assert len(dgr_seqs) == len(set(dgr_seqs))
                
assert os.stat(dgr_ad_norm_out).st_size > 0


# %%
