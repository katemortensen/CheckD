#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 8/25/2024
# Purpose: Determine MAG in cohort with least diversity 


# %% 
import argparse
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn.preprocessing import StandardScaler
import umap
import seaborn as sns 
from matplotlib.colors import Normalize
import numpy as np
import matplotlib.gridspec as gridspec
import os


# %% 
# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ad_norm_cohort_stats_dir", "--ad_norm_cohort_stats_dir", metavar="<mag>.dgr_ad_norm_file", help="specify the path where input files are (<mag>.dgr_ad_norm_file)", required=True)
parser.add_argument("-metric", "--metric", metavar="metric", help="diversity metric basis (mean_ad_norm)", required=True)
parser.add_argument("-basename", "--basename", metavar="basename", help="output basename", required=True)
parser.add_argument("-output_dir", "--output_dir", metavar="output_dir", help="specify the path where the outputs will go", required=True)

args = parser.parse_args()
# ad_norm_file = args.ad_norm_file
ad_norm_cohort_stats_dir = args.ad_norm_cohort_stats_dir
output_dir = args.output_dir
basename = args.basename
mag_stats_file = f'{ad_norm_cohort_stats_dir}/{basename}.mag_stats'

file_out = f'{output_dir}/{basename}.cohort_mag_diversity_baseline'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# delete this hardcoded comment before running

# %% 
# Test Vars

# humanO1
# basename = 'humanO1'
# metric = 'mean_ad_norm'
# metric = 'stdev_ad_norm'
# metric = 'prcnt_lessthan1'

##################### 

# wkdir = f'/home/kmorten/{basename}_CheckD/'
# ad_norm_cohort_stats_dir = f'{wkdir}/ad_norm_cohort_stats'
# mag_stats_file = f'{ad_norm_cohort_stats_dir}/{basename}.mag_stats'
# output_dir = f'{ad_norm_cohort_stats_dir}'
# file_out = f'{output_dir}/{basename}.cohort_mag_diversity_baseline'


#####################

# %% 
df = pd.read_csv(mag_stats_file, skiprows=11, sep='\t')
df.insert(loc=2, column='diversity', value='NA')
df.insert(loc=3, column='distance', value='NA')

metric = 'mean_ad_norm'
df['diversity'] = 1-df[metric]
df = df.sort_values(by='diversity', ascending=True).reset_index(drop=True)
baseline = df['diversity'].min()
df['distance'] = df['diversity'] - baseline

df['diversity'] = df['diversity'].round(4)
df['distance'] = df['distance'].round(4)

with open(file_out, 'w') as fp_out:
    header = '# Allelic depth for all positions.\
        \n# Note that the normalized allelic depth of a reference nucleotide\
        \n# at a position will henceforth be reffered to as "ad_norm".\
        \n# Diversity = 1-mean_ad_norm.\
        \n# The MAG with the lowest diversity is the baseline diversity of a cohort. \
        \n# The distance from the baseline is the diffence between a diversity score and the baseline.\
        \n# (1) MAG name\
        \n# (2) total length of assembly (nucleotides)\
        \n# (3) diversity (1-mean_ad_norm)\
        \n# (4) distance (diversity score - baseline diversity score)\
        \n# (5) mean allelic depth, normalized\
        \n# (6) median allelic depth, normalized\
        \n# (7) standard deviation of allelic depth, normalized values\
        \n# (8) percent of ref positions with < 100% of mapped reads equal to ref nucleotide, aka allelic depth, normalized\
        \n# cols: mag, tot_len, diversity, distance, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1\
        \n#\n'
        # \n# (9) percent CDS of assembly\
        # \n# (10) percent non-CDS of assembly\
        # \n# (11) percent DGR of assembly\
        # \n# (12) percent ncRNA of assembly\
        # \n# (13) percent CRISPR of assembly\
        #\n# cols: mag, tot_len, diversity, distance, mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1, prcnt_cds, prcnt_noncds, prcnt_dgr, prcnt_ncrna, prcnt_crispr\
        # \n# mag\ttot_len\tdiversity\tdistance\tmean_ad_norm\tmedian_ad_norm\tstdev_ad_norm\tprcnt_lessthan1'
    fp_out.write(header)
    # df.rename(columns={'mag': '# mag'}, inplace=True)
    df.to_csv(fp_out, index=False, sep="\t")
    
# %% YOU LEFT OFF HERE - MAKE THE FOLLOWING A SEPARATE FILE 

'''creat_magstats_dict'''

def create_magstats_dict(domain_file_dict): 
    
    ''' magDict '''
    ''' mag > tot_len, prcnt_cds, mead_ad_norm, ... prcnt_dgr, prcnt_cds... '''

    magDict = {}

    count = 0
    for domain in domain_file_dict.keys():
        if count == 0:
            assert domain == 'all'
        count += 1
        print(domain)
        file_in = domain_file_dict.get(domain)
        with open(file_in, 'r') as fp_in:
            for line in fp_in:
                if line.startswith('# cols: '):
                    info = line.strip().split('# cols: ')[1].split(', ')
                if line.startswith('#') is False:
                    linep = line.strip().split('\t')
                    tmp_dict = dict(zip(info, linep))
                    assembly = tmp_dict.get('mag')
                    tmp_dict.pop('mag')
                    assert 'mean_ad_norm' in list(tmp_dict.keys())
                    if magDict.get(f'{assembly}') is None:
                        if domain != 'all':
                            print(domain)
                            print(assembly)
                        assert domain == 'all'
                        magDict[assembly] = tmp_dict
                    elif domain != 'all':
                        feature = f'prcnt_{domain}'
                        feature_value = tmp_dict.get(feature)
                        if magDict[assembly].get(feature) is None:
                            magDict[assembly][feature] = feature_value
                        else:
                            assert feature_value == magDict[assembly].get(feature)

    return magDict

diversity_file= f'{output_dir}/{basename}.cohort_mag_diversity_baseline'
domains = ['all', 'cds','noncds','dgr','ncrna','crispr']
file_list = []
for domain in domains:
    if domain == 'all':
        file_list.append(diversity_file)
    else:
        file_in = f'{ad_norm_cohort_stats_dir}/{basename}.{domain}_mag_stats'
        file_list.append(file_in)
domain2file = dict(zip(domains,file_list))
magstatsDict = create_magstats_dict(domain2file)

df = pd.DataFrame.from_dict(magstatsDict, orient='index')
df.index.name = 'mag'
df = df.reset_index()
assert not df.isna().any().any(), "AssertionError: DataFrame contains NA values"
df.head()                  

                      
# %%
# UMAP Plots

# scaler = StandardScaler()
# cols = ['diversity','distance','stdev_ad_norm','prcnt_lessthan1']
# for domain in domains[1:]:
#     cols.append(f'prcnt_{domain}')

# label = 'distance'

# # df_scaled = scaler.fit_transform(df.iloc[:, [2,5,6]]) 
# df_scaled = scaler.fit_transform(df.loc[:, cols]) 

# umap_model = umap.UMAP(n_neighbors=15, 
#                        n_components=len(cols), 
#                        metric='cosine')
# umap_results = umap_model.fit_transform(df_scaled)


# # Creating the plot
# fig, ax = plt.subplots()

# # print(df[label].isna().sum())
# df[label] = df[label].astype(float)

# norm = Normalize(vmin=np.min(df[label].values), 
#                  vmax=np.max(df[label].values))

# # Scatter plot
# scatter = ax.scatter(x=umap_results[:, 0], y=umap_results[:,1], 
#                      c=df[label].values, cmap='viridis', 
#                      edgecolor='k', norm=norm)

# # Add a title and labels
# ax.set_title('Scatter Plot with Gradient Scale')
# ax.set_xlabel('UMAP 1')
# ax.set_ylabel('UMAP 2')

# # Create a colorbar
# cbar = plt.colorbar(scatter, ax=ax, pad=0.01)
# cbar.set_label(label)
# cbar.ax.tick_params(labelsize=10)  # Adjust tick size if necessary

# # Adjust layout to make space for the colorbar
# plt.tight_layout()

# plt.show()


#########################################
# %% 

# Correlation Plots 
# for i in cols:

#     variable1 = 'diversity'
#     variable2 = i

#     # Perform linear regression
#     slope, intercept, r_value, p_value, std_err = stats.linregress(df[f'{variable1}'], df[f'{variable2}'])

#     # Calculate R-squared
#     r_squared = r_value**2

#     # Print the R-squared value
#     print(f'R-squared: {r_squared:.2f}')


#     # Create the plot again to annotate
#     plt.figure(figsize=(10, 6))
#     sns.regplot(x=f'{variable1}', y=f'{variable2}', data=df, scatter_kws={'s': 50}, line_kws={'color': 'red'})

#     # Add text annotation for R-squared
#     plt.text(x=max(df[f'{variable1}']) * 0.8, y=min(df[f'{variable2}']), s=f'R² = {r_squared:.2f}', 
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=12)

#     plt.title('Correlation of Domain Percentage and Diversity')
#     plt.xlabel(f'{variable1}')
#     plt.ylabel(f'{variable2}')
#     plt.show()
    

# %%