#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/17/2024
# Purpose: FragGeneScan CDS region, exploratory data analysis of mapping diversity stats

# %% 
import argparse
import os

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="path-to-reference-positional-frequency", help="path to <bin>.ad_norm", required=True)
parser.add_argument("-fgs_file", "--fgs_file", metavar="path-to-mapping-diversity-stats", help="path to <bin>-fgs.out", required=True)
parser.add_argument("-ffn_file", "--ffn_file", metavar="path-to-fgs-ffn", help="path to <bin>-fgs.ffn", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="assembly name as an output prefix", required=True)

args = parser.parse_args()
fgs_file = args.fgs_file
ffn_file = args.ffn_file
ad_norm_file = args.ad_norm_file
output_dir = args.output_dir
assembly_name = args.assembly_name
mag_name = assembly_name 

cds_output_file = f'{output_dir}/{mag_name}.cds_ad_norm'
noncds_output_file = f'{output_dir}/{mag_name}.noncds_ad_norm'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# delete this hardcoded comment before running

# %% 
# temp vars for testing 

# mag_name = 'bin.195.strict'
# mag_name = 'bin.311.strict'
# mag_name = 'bin.213.permissive'
# mag_name = 'bin.339.strict'
# mag_name = 'bin.388.strict'
# mag_name = 'bin.39'
# mag_name = 'bin.405.strict'
# test = '>NODE_324_length_541_cov_0.784483'
# mag_name = 'bin.0'

# illumina example
# assembly_name = mag_name
# wkdir = f'/home/kmorten/CheckD/illumina_example/illumina_example_checkd/'
# fgs_file = f'{wkdir}/fgs_out/{mag_name}-fgs.out'
# ffn_file = f'{wkdir}/fgs_out/{mag_name}-fgs.ffn'
# ad_norm_file = f'{wkdir}/ad_norm/{mag_name}.ad_norm'
# cds_output_file = f'{wkdir}/ad_norm/{mag_name}.cds_ad_norm'

# # ADMB
# fgs_file =  f'/home/kmorten/ADMB_CheckD/METAWRAP_MAGS/fgs_out/{mag_name}-fgs.out'
# ffn_file = f'/home/kmorten/ADMB_CheckD/METAWRAP_MAGS/fgs_out/{mag_name}-fgs.ffn'
# ad_norm_file = f'/home/kmorten/ADMB_CheckD/METAWRAP_MAGS/pos_freq/{mag_name}.ad_norm'
# cds_output_file = f'/home/kmorten/ADMB_CheckD/METAWRAP_MAGS/pos_freq/{mag_name}.cds_ad_norm'
# noncds_output_file = f'/home/kmorten/ADMB_CheckD/METAWRAP_MAGS/pos_freq/{mag_name}.noncds_ad_norm'

# # humanO1
# output_dir = f'/home/kmorten/humanO1_CheckD/ad_norm'
# fgs_file =  f'/home/kmorten/humanO1_CheckD/fgs_out/{mag_name}-fgs.out'
# ffn_file = f'/home/kmorten/humanO1_CheckD/fgs_out/{mag_name}-fgs.ffn'
# ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'
# cds_output_file = f'{output_dir}/{mag_name}.cds_ad_norm'


# %% 

'''cds_ffn_dict'''
'''contig > seq'''

cds_ffn_dict = {}

with open(ffn_file,'r') as file:
    for line in file:
        if line.startswith('>'):
            contig_name = line.strip()
        else:
            seq = line.strip()
            cds_ffn_dict[contig_name] = seq

'''indel_dict'''
'''contig > inserts positions, deletion positions, indel positions'''

indel_dict = {}
with open(fgs_file, 'r') as fp_in:
    for line in fp_in:
        if line.startswith('>'):
            contig_name = line.strip()
        else:
            linep = line.strip().split('\t')
            I = linep[-2]
            ins = list(filter(None, I.strip('I:').split(',')))
            D = linep[-1]
            dels = list(filter(None, D.strip('D:').split(',')))
            indels = ins + dels 
            if indel_dict.get(contig_name) == None:
                indel_dict[contig_name] = {'inserts':ins,'deletions':dels,'indels':indels}
            else:
                new_ins = ins + indel_dict[contig_name].get('inserts')
                new_dels = dels + indel_dict[contig_name].get('deletions')
                new_indels = indels + indel_dict[contig_name].get('indels')
                indel_dict[contig_name] = {'inserts':new_ins,'deletions':new_dels,'indels':new_indels}
    

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

'Convert Nucleotide'

def convert_nucl(nucleotide, direction):
    nucl_dict = {'A':'T', 'T':'A','C':'G','G':'C','N':'N'}
    if direction == '+':
        result = nucleotide
    elif direction == '-':
        result = nucl_dict.get(nucleotide)
    else:
        result = f'Issue with nucleotide of direction of strand ({assembly_name})'
    return(result)


'Stop Codon Check'

def stop_check(seq):
    last_codon = seq[-3:]
    stop_codons = ['TAA','TAG','TGA']
    result = last_codon in stop_codons
    return(result)


'CDS & Codon positions'

with open(cds_output_file, 'w') as cds_out:
    cds_out_header = '# Coding domain sequence (CDS) info.\
        \n# Note: If ref_nuc, ad_norm, and depth are "NA", then pos was not mapped in VCF file but accounted for in CDS position by FragGeneScan. \
        \n# The below information is separated by ">contigID_start_end_strand;desc" where strand is forward (+) or reverse (-) and the description (desc) is included if available. \
        \n# INDELS are identified by FragGeneScan are excluded\
        \n# Positions shared by multiple CDS regions are entered redundantly to account for the difference in nucleotide and codon position.\
        \n# (1) nucleotide position in contig\
        \n# (2) reference nucleotide\
        \n# (3) allelic depth, normalized [0,1] (high-quality bases only)\
        \n# (4) depth (high-quality bases only)\
        \n# (5) codon position (1, 2, or 3)\
        \n# cols: pos, ref_nuc, ad_norm, depth, codon_pos\n'
    cds_out.write(cds_out_header)
    with open(ffn_file,'r') as ffn_in:
        for line in ffn_in:
            if line.startswith('>'):
                contig_header = line.strip()
                l = len(contig_header.split('_'))
                contig_name = '_'.join(contig_header.strip().split('_')[0:l-3])
                linep = line.strip().split('_')
                '''Directionality'''
                direction = linep[-1]
                if direction == '-':
                    start = int(linep[-2])
                    end = int(linep[-3])
                    step = -1
                elif direction == '+':
                    start = int(linep[-3])
                    end = int(linep[-2])
                    step = 1     
                cds_out.write(f'\n{contig_header}')
            else:
                seq = line.strip()                
                codon_pos = 0
                count = 0
                for i in range(start, end + step, step):
                    pos = str(i)
                    nucl = seq[count]
                    '''Read Count Acceptability'''
                    # Note: Quality Scores are only given for positions that meet mpileup parameter reqs
                    if ad_norm_dict.get(contig_name) != None:
                        if ad_norm_dict[contig_name].get(pos) != None:
                            ad_norm = ad_norm_dict[contig_name][pos].get('ad_norm')    
                            depth = ad_norm_dict[contig_name][pos].get('depth')    
                            ref_nuc = ad_norm_dict[contig_name][pos].get('ref_nuc')
                            ref_nuc_check = convert_nucl(nucl,direction)
                            if ref_nuc != 'N':                                 
                                assert convert_nucl(nucl,direction) == ad_norm_dict[contig_name][pos].get('ref_nuc')
                            else:
                                assert float(ad_norm) == float('0')
                        else:
                            ad_norm = 'NA'
                            ref_nuc = 'NA'
                            depth = 'NA'
                            # print(f'Position was not mapped in VCF file: {assembly_name},{contig_name},{pos}')
                    '''Codon Positions'''
                    if pos not in indel_dict[contig_name].get('indels'):
                        if codon_pos == 3:
                            codon_pos = 1
                        else:
                            codon_pos += 1
                        new_line = f'\n{pos}\t{ref_nuc}\t{ad_norm}\t{depth}\t{codon_pos}'
                        cds_out.write(new_line)    
                    '''Deletions & Codon Positions'''
                    # Note: Inserts detected by fgs are skipped
                    # Note: positions marked as deletions by fgs are adjusted by an extra position
                    if pos in indel_dict[contig_name].get('deletions'):
                        if codon_pos == 1:
                            codon_pos = 3
                        else:
                            codon_pos += -1 
                        new_line = f'\n{pos}\t{ref_nuc}\t{ad_norm}\t{depth}\t{codon_pos}'
                        cds_out.write(new_line)
                    count += 1
                assert codon_pos == 3



# %%
