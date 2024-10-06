#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/2/2024
# Purpose: EDA cohort_mag_desc_stats output


# %% 
import argparse
import os
import pandas as pd 
import numpy as np
from checkd_functions import create_ad_norm_dict
from checkd_functions import create_general_stats_dict
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-sliding_window_dir", "--sliding_window_dir", metavar="path-to-sliding_window_file", help="path to sliding window director", required=True)
parser.add_argument("-metric", "--metric", metavar="metric", help="mead_ad_norm, stdev, or prcnt_lessthan1", required=True)
parser.add_argument("-window", "--window", metavar="window size", help="window size", required=True)
parser.add_argument("-step", "--step", metavar="step size", help="step size", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="ex: bin.62", required=True)
parser.add_argument("-ad_norm_dir", "--ad_norm_dir", metavar="path-to-ad_norm-files", help="path to <bin>.ad_norm", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="path-to-output_dir", help="output dir", required=True)

args = parser.parse_args()
sliding_window_dir = args.sliding_window_dir
window = int(args.window)
step = int(args.step)
metric = args.metric
assembly_name = args.assembly_name
ad_norm_dir = args.ad_norm_dir
output_dir = args.output_dir

if 'mean' in metric:
    inverse_metric = True
if 'median' in metric:
    inverse_metric = True
if 'stdev' in metric:
    inverse_metric = False
if 'prcnt' in metric:
    inverse_metric = False 

window_file = f'{sliding_window_dir}/{assembly_name}.w{window}_s{step}'
dgr_ad_norm_file = f'{ad_norm_dir}/{assembly_name}.dgr_ad_norm'
crispr_ad_norm_file = f'{ad_norm_dir}/{assembly_name}.crispr_ad_norm'
ncrna_ad_norm_file = f'{ad_norm_dir}/{assembly_name}.ncrna_ad_norm'
vir_ad_norm_file = f'{ad_norm_dir}/{assembly_name}.vir_ad_norm'

ad_norm_files = [dgr_ad_norm_file, 
                 crispr_ad_norm_file, 
                 ncrna_ad_norm_file, 
                 vir_ad_norm_file]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
    
# delete this hardcoded comment before running

# %% 
# Test vars

# human01
# wkdir = f'/home/kmorten/humanO1_CheckD/'
# sliding_window_dir = f'{wkdir}/sliding_window/'
# assembly_name = 'bin.62'
# metric = 'mean_ad_norm'
# # metric = 'stdev_ad_norm'
# # metric = 'prcnt_lessthan1'
# window = 50000
# step = 40000
# window_file = f'{sliding_window_dir}/{assembly_name}.w{window}_s{step}'

# # domains: CDS, nonCDS, DGR, crispr, ncRNA
# # cds_ad_norm_file = f'{wkdir}/ad_norm/{assembly_name}.cds_ad_norm'
# # noncds_ad_norm_file = f'{wkdir}/ad_norm/{assembly_name}.noncds_ad_norm'
# dgr_ad_norm_file = f'{wkdir}/ad_norm/{assembly_name}.dgr_ad_norm'
# crispr_ad_norm_file = f'{wkdir}/ad_norm/{assembly_name}.crispr_ad_norm'
# ncrna_ad_norm_file = f'{wkdir}/ad_norm/{assembly_name}.ncrna_ad_norm'
# vir_ad_norm_file = f'{ad_norm_dir}/{assembly_name}.vir_ad_norm'

# ad_norm_files = [dgr_ad_norm_file, 
#                  crispr_ad_norm_file, 
#                  ncrna_ad_norm_file, 
#                  vir_ad_norm_file]

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)
    
# %% 

'''winDict'''
'''contig > pos > mean, median, stdev, prcnt_lessthan1'''

winDict = {}
contig_pos_Dict = {}
with open(window_file, 'r') as fp_in:
    for line in fp_in:
        if line.startswith('# cols: '):
            info = line.strip().split('# cols: ')[-1].split(', ')
        if line.startswith('#') is False:
            linep = line.strip('\n').split('\t')
            tmp_dict = dict(zip(info,linep))    
            contig = tmp_dict.get('contig')
            tmp_dict.pop('contig')
            pos = tmp_dict.get('pos')
            tmp_dict.pop('pos')
            if winDict.get(contig) == None:
                winDict[contig] = {pos:tmp_dict}
            else: 
                assert winDict[contig].get(pos) == None
                winDict[contig][pos] = tmp_dict
                

'''domainDict'''
'''contig > seq > start, end, label'''


key_map = {'dgr_ad_norm': 'DGR',
           'crispr_ad_norm': 'CRISPR',
           'ncrna_ad_norm': 'ncRNA', 
           'vir_ad_norm': 'VIR'}

domainDict = {}
for i in ad_norm_files:
    domain_ad_norm = i.strip().split('.')[-1]
    domain = key_map.get(domain_ad_norm)
    with open(i, "r") as fp_in:
        for line in fp_in:
            if line.startswith('>'):
                seq_header = line.strip()
                linep = line.strip().strip('>').split('_')
                contig = linep[0]
                strand = line.split(';')[0].split('_')[-1]
                label = seq_header.split(';')[-1]
                assert strand in ['-','+']
                if strand == '-':
                    start = max(int(linep[1]),int(linep[2]))
                    end = min(int(linep[1]),int(linep[2]))
                elif strand == '+':
                    start = min(int(linep[1]),int(linep[2]))
                    end = max(int(linep[1]),int(linep[2]))             
                tmp_dict = {'strand':strand,
                            'start':start,
                            'end':end,
                            'label':label}
                if domainDict.get(domain) is None:
                    domainDict[domain] = {contig:{seq_header:tmp_dict}}
                elif domainDict[domain].get(contig) is None:
                    domainDict[domain][contig] = {seq_header:tmp_dict}
                else:
                    domainDict[domain][contig][seq_header] = tmp_dict
                    

# %% 

def create_plot_by_metric(winDict, domainDict, metric, inverse_metric):
    data = []
    contigs = []
    for contig in winDict:
        contigs.append(contig)
        for pos in winDict[contig].keys():
            measurement = winDict[contig][pos].get(metric)
            tmp_list = [contig, pos, measurement]
            data.append(tmp_list)

    df = pd.DataFrame(data, columns=["Contig", "Position", metric])
    df['Position'] = pd.to_numeric(df['Position'])
    df[metric] = df[metric].astype(float)

    if inverse_metric is True:
        df['diversity'] = 1-df[metric]
        ylabel = f'Diversity\n(1-{metric})'
    else:
        df['diversity'] = df[metric]
        ylabel = f'Diversity\n({metric})'


    domains = ['DGR', 'CRISPR', 'ncRNA', 'VIR']
    domain_colors = ['green', 'red', 'purple', 'goldenrod']
    custom_palette = dict(zip(domains, domain_colors))

    for contig in winDict.keys():

        plot_name = f'{assembly_name}_{contig}'
        fig, ax = plt.subplots(figsize=(14, 7))

        contig_data = df[df['Contig'] == contig]
        ax.plot(contig_data['Position'], contig_data['diversity'], 
                label=contig, color='gray', linewidth = 2.5)

        # Shade specific position ranges for domains and manage legend entries
        legend_handled = {}  
        for domain in domainDict:
            if domain in custom_palette:
                domain_color = custom_palette[domain]
                if domain in domainDict and domainDict[domain].get(contig):
                    for j in domainDict[domain][contig]:
                        start = int(domainDict[domain][contig][j]['start'])
                        end = int(domainDict[domain][contig][j]['end'])
                        if end > int(contig_data['Position'].iloc[-1]):
                            print(int(contig_data['Position'].iloc[-1]))
                            print(start)
                            print(end)
                            print(j)
                        ax.axvspan(start, end, color=domain_color, alpha=0.3)
                        
                        # Check if the domain is already added to the legend
                        if domain not in legend_handled:
                            legend_handled[domain] = mpatches.Patch(color=domain_color, label=domain)

        # set global font sizes
        plt.rcParams['font.size'] = 20  # Sets the default font size
        plt.rcParams['axes.titlesize'] = 20  # Font size for titles
        plt.rcParams['axes.labelsize'] = 15  # Font size for x and y labels
        plt.rcParams['xtick.labelsize'] = 15  # Font size for x-tick labels
        plt.rcParams['ytick.labelsize'] = 15  # Font size for y-tick labels
        plt.rcParams['legend.fontsize'] = 15  # Font size for legend
        
        # Assembling the legend from handled entries
        window_text = window/1000
        step_text = step/1000
        info_patch = mpatches.Patch(color='none', label=f'Window = {window_text}kb, Step = {step_text}kb')
        ax.legend(handles=list(legend_handled.values()) + 
                [mpatches.Patch(color='gray', label=f'Diversity in {contig}'), info_patch],
                loc='upper left', bbox_to_anchor=(1, 1))
        
        ax.set_title(f'{assembly_name}; {contig}')
        ax.set_xlabel(f'Position in {contig} (kb)')
        x_min, x_max = df['Position'].min(), df['Position'].max()
        ax.set_xticklabels([f"{int(tick/1000)} kb" for tick in ax.get_xticks()])
        ax.set_ylabel(ylabel)
        y_min, y_max = df['diversity'].min(), df['diversity'].max()
        ax.set_yticks(np.linspace(y_min, y_max, num=10))

        plot_out = f'{output_dir}/{assembly_name}_{contig}_w{window}_s{step}_{metric}'
        plt.savefig(f'{plot_out}.png', format="png", bbox_inches="tight")
        plt.savefig(f'{plot_out}.pdf', format="pdf", bbox_inches="tight")
        plt.show(block=False)
        plt.pause(2)  # Pause for 2 seconds before showing the next plot
        plt.close(fig)  # Close the plot window to avoid overlapping
            
        plt.show()

# %%
create_plot_by_metric(winDict, domainDict, metric, inverse_metric)


# %%
