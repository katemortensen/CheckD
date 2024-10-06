#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 8/22/2024
# Purpose: Determine MAG in cohort with least diversity 


# %% 
import argparse
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.stats as stats

# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)


parser.add_argument("-ad_norm_file", "--ad_norm_file", metavar="<mag>.dgr_ad_norm_file", help="specify the path where input files are (<mag>.dgr_ad_norm_file)", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)
parser.add_argument("-assembly_name", "--assembly_name", metavar="assembly_name", help="output assembly_name", required=True)

args = parser.parse_args()
ad_norm_file = args.ad_norm_file
output_dir = args.output_dir
assembly_name = args.assembly_name
mag_desc_stats_out = f'{output_dir}/{assembly_name}.dgr_mag_desc_stats'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


delete this hardcoded comment before running

# %% 
# Test Vars

# humanO1
basename = 'humanO1'
metric = 'mean_ad_norm'
# metric = 'stdev_ad_norm'
# metric = 'prcnt_lessthan1'

##################### 

wkdir = f'/home/kmorten/{basename}_CheckD/'
ad_norm_cohort_stats_dir = f'{wkdir}/ad_norm_cohort_stats'
mag_stats_file = f'{ad_norm_cohort_stats_dir}/{basename}.mag_stats'


#####################

# %% 
df = pd.read_csv(mag_stats_file, skiprows=11, sep='\t')
df.insert(loc=1, column='diversity', value='NA')
inverse_metric = False
if ('mean' in metric) or ('median' in metric):
    inverse_metric = True
    df['diversity'] = 1-df[metric]
    ylabel = f'Diversity (1-{metric})'
else:
    df['diversity'] = df[metric]
    ylabel = f'Diversity ({metric})'

df = df.sort_values(by='diversity', ascending=True)

    
fig = sm.qqplot(df['diversity'])  # 'line="45"' adds a 45-degree reference line
ax = plt.gca()
plt.title(f'{basename} Diversity QQ-Plot')
ax.set_ylabel(f'MAG Quantiles\n{ylabel}')

plot_out = f'{output_dir}/{basename}_qqplot'
plt.savefig(f'{plot_out}.png', format="png", bbox_inches="tight")
plt.savefig(f'{plot_out}.pdf', format="pdf", bbox_inches="tight")
plt.show(block=False)
plt.pause(2)  # Pause for 2 seconds before showing the next plot
plt.close(fig)  # Close the plot window to avoid overlapping
        
plt.show()


# %%
# # QQ-Plot against a t-distribution
# n = len(df)
# fig = sm.qqplot(df['mean_ad_norm'], dist=stats.t, distargs=(n,))  # 'distargs' sets the degrees of freedom for the t-distribution
# plt.title("QQ-Plot vs. t-Distribution")
# plt.show()

# %%
