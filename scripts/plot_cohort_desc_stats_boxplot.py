#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 10/26/2024
# Purpose: box & whiskers plot of allelic depth stats 

# %% 
import argparse
import os
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator

# %% 
# # arguments
parser = argparse.ArgumentParser()
parser.add_argument("-cohort_desc_stats_file", "--cohort_desc_stats_file", metavar="<basename>.cohort_desc_stats", help="specify the path where input files are (<basename>.cohort_desc_stats)", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-basename", "--basename", metavar="basename", help="output basename", required=True)

args = parser.parse_args()
cohort_desc_stats_file = args.cohort_desc_stats_file 
output_dir = args.output_dir
basename = args.basename
in_extension = cohort_desc_stats_file.split('.')[-1]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# delete this hardcoded comment before running

# %% 

# Test vars

# humanO1
# wkdir = '/home/kmorten/humanO1_CheckD/'
# basename = 'humanO1'
# cohort_desc_stats_file = f'{wkdir}/ad_norm_cohort_stats/{basename}.cohort_desc_stats'


# humanO1
# wkdir = '/home/kmorten/sheepA_checkd/'
# basename = 'sheepA'
# cohort_desc_stats_file = f'{wkdir}/ad_norm_cohort_stats/{basename}.cohort_desc_stats'


# ##################
# output_dir = f'{wkdir}/plots/ad_norm_cohort_stats_plots'
# in_extension = cohort_desc_stats_file.split('.')[-1]

# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)
    
# %%
# Collect data
with open(cohort_desc_stats_file, 'r') as fp_in:
    count = 0
    for line in fp_in:
        if line.startswith('# cols:'):
            info = line.strip().split('# cols: ')[-1].split(', ')
            list_of_data = []
            count += 1
        if (count == 1) & (line.startswith('#') == False):
            linep = line.strip('\n').split('\t')
            if 'NA' not in linep:
                list_of_data.append(linep)

df = pd.DataFrame(list_of_data, columns=info)

df.loc[df['type'] == 'ncrna', 'type'] = 'ncRNA'
df.loc[df['type'] == 'crispr', 'type'] = 'CRISPR'
df.loc[df['type'] == 'dgr', 'type'] = 'DGR'
df.loc[df['type'] == 'vir', 'type'] = 'VIR'



# %% 
def create_plot_by_metric(df, metric, plot_out, inverse_metric=True):
    
    df[metric] = df[metric].astype(float)
    
    if inverse_metric is True:
        df['diversity'] = 1-df[metric]
        ylabel = f'Diversity\n(1-{metric})'
    else:
        df['diversity'] = df[metric]
        ylabel = f'Diversity\n({metric})'

    # Plotting params
    sns.set_theme(style='whitegrid')
    sns.boxplot(x = 'type', y = f'diversity', 
                    data = df,
                    linewidth=2.5)
    # adding data points 
    sns.stripplot(x='type', y='diversity', data=df, color="orange", alpha = 0.6) 
    plt.title(f'{basename}', fontsize=20)
    plt.ylabel(f'{ylabel}', fontsize = 15)
    plt.xticks(rotation=45,horizontalalignment='right', )
    plt.xlabel(f'Sequence Type', fontsize = 15)
    ax = plt.gca()
    ax.yaxis.set_major_locator(MaxNLocator(integer=False))

    plt.savefig(f'{plot_out}.png', format="png", bbox_inches="tight")
    plt.savefig(f'{plot_out}.pdf', format="pdf", bbox_inches="tight")
    plt.show()
    plt.close()

    return


# %%

metric = 'mean_ad_norm'
plot_out = f'{output_dir}/{basename}.{in_extension}.{metric}.boxplot'
create_plot_by_metric(df, metric, plot_out, inverse_metric=True)

metric = 'stdev_ad_norm'
plot_out = f'{output_dir}/{basename}.{in_extension}.{metric}.boxplot'
create_plot_by_metric(df, metric, plot_out, inverse_metric=False)

metric = 'prcnt_lessthan1'
plot_out = f'{output_dir}/{basename}.{in_extension}.{metric}.boxplot'
create_plot_by_metric(df, metric, plot_out, inverse_metric=False)


# %%
