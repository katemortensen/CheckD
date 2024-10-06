#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/2/2024
# Purpose: Codons, exploratory data analysis of mapping diversity stats

# %% 
import argparse
import os
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_dir", "--ad_norm_dir", metavar="ad_norm_dir", help="specify the path where input files are (<mag>-fgs.ffn, <mag>-fgs.out, <mag>.ad_norm)", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-basename", "--basename", metavar="basename", help="output basename", required=True)

args = parser.parse_args()
ad_norm_dir = args.ad_norm_dir
output_dir = args.output_dir
basename = args.basename
mag_name = basename 
codon_dist_plot_out = f'{output_dir}/{basename}.codon_dist_plot'
codon_boxplot_plot_out = f'{output_dir}/{basename}.codon_boxplot_plot'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 
# Test Vars

# # illumina example
# basename = 'illumina_example'
# wkdir = '/home/kmorten/CheckD/illumina_example/illumina_example_checkd/'
# ad_norm_dir = f'{wkdir}/ad_norm'
# output_dir = f'{wkdir}/plots/cohort_plots/'

# # humanO1
# basename = 'humanO1'
# wkdir = f'/home/kmorten/{basename}_CheckD/'
# ad_norm_dir = f'{wkdir}/ad_norm/'
# output_dir = f'{wkdir}/plots/cohort_plots/'

#################### 

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

codon_dist_plot_out = f'{output_dir}/{basename}.codon_dist_plot'
codon_boxplot_plot_out = f'{output_dir}/{basename}.codon_boxplot_plot'

# %% 
file_names = []
for file in os.listdir(ad_norm_dir):
    if file.endswith('.cds_ad_norm'):
        file_names.append(file)

# data = []
tmp_file = f'{output_dir}/.cds_plot_tmp'
with open(tmp_file, 'w') as tmp:
    tmp.write(f'codon\tad_norm')
    for file in file_names:
        mag_name = file.strip().split('/')[-1].split('.cds_ad_norm')[0]  
        print(mag_name)
        cds_ad_norm_file = f'{ad_norm_dir}/{mag_name}.cds_ad_norm'
        with open(cds_ad_norm_file, 'r') as fp_in:
            next(fp_in)
            next(fp_in)
            for line in fp_in:
                if (line.startswith('>') == False) and (line.startswith('#') == False) and (line.startswith('\n') == False) and (line.startswith('N')==False):
                    linep = line.strip().split('\t')
                    if linep[2] not in ['N', 'NA', 'NaN', None]:
                        ad_norm = float(linep[2])
                        codon = int(linep[3])
                        # data.append({'codon':codon,'ad_norm':ad_norm})  
                        tmp.write(f'\n{codon}\t{ad_norm}')

# df = pd.DataFrame(data)
# df['ad_norm'] = df['ad_norm'].astype(float)

tmp_file_path = f'{output_dir}/.cds_plot_tmp'

df = pd.read_csv(tmp_file_path, delimiter='\t')

# %% 
sns.set_theme()
codon_order = df.codon.value_counts().index
sea = sns.FacetGrid(df, row = "codon", 
                    row_order = codon_order,
                    height = 1.7, 
                    aspect = 4)
# sea.map(sns.histplot, "ad_norm", kde=True, stat="density")
sea.map(sns.kdeplot, "mean_ad_norm")
sea.set_xlabels("mean_ad_norm")
sea.figure.subplots_adjust(top=0.85)
sea.figure.suptitle(f'{basename}')
plt.savefig(f'{codon_dist_plot_out}.pdf', format="pdf", bbox_inches="tight")
plt.savefig(f'{codon_dist_plot_out}.png', format="png", bbox_inches="tight")
plt.show()



# %%
# sns.set_theme()

# sns.violinplot(y=df['mean_ad_norm'], x=df['codon'], data=df)
# plt.title(f'{basename}', fontsize=20)
# plt.ylabel('ad_norm', fontsize = 15)
# plt.xticks(rotation=45,horizontalalignment='right')
# # plt.xlabel(f'{}', fontsize = 15) 
# plt.show()
# %%
