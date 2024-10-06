#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 8/14/2024
# Purpose: Explore correlation between ncRNA, CRISPR, and DGR volume
# (percent of genome, instance count) and MAG diversity 


# %% 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt


# %% 
# Test vars

# human01
basename = 'humanO1'


wkdir = f'/home/kmorten/{basename}_CheckD/'
ad_norm_cohort_stats_dir = f'{wkdir}/ad_norm_cohort_stats/'
output_dir = ad_norm_cohort_stats_dir
domain_stats_out = f'{output_dir}/{basename}.mag_domain_stats'

# %% 
fp = f'{ad_norm_cohort_stats_dir}/{basename}.mag_stats'
df_mag_stats = pd.read_csv(fp, skiprows=11, sep='\t')
df_mag_stats.rename(columns={'# mag': 'mag'}, 
                    inplace = True)

df = pd.DataFrame({'mag':[],
                   'domain':[],
                   'prcnt_domain':[],
                   'domain_count':[]})

domains = ['dgr','crispr','ncrna']

for i in range(len(domains)):
    domain = domains[i]
    fp = f'{ad_norm_cohort_stats_dir}/{basename}.{domain}_mag_stats'

    columns_to_include = ['# mag', f'prcnt_{domain}',f'{domain}_count']
    df_tmp = pd.read_csv(fp, skiprows=17, sep='\t',usecols = columns_to_include)
    df_tmp.rename(columns={'# mag': 'mag',
                        f'prcnt_{domain}': 'prcnt_domain', 
                        f'{domain}_count': 'domain_count'}, 
                        inplace = True)
    df_tmp['domain'] = len(df_tmp)*[f'{domain}'] 
    df = pd.concat([df, df_tmp], ignore_index=True)

assert len(list(set(df['mag'])))*len(domains) == len(df)

df = pd.merge(df_mag_stats, df, on='mag', how='outer')
df.loc[df['domain'] == 'ncrna', 'domain'] = 'ncRNA'
df.loc[df['domain'] == 'crispr', 'domain'] = 'CRISPR'
df.loc[df['domain'] == 'dgr', 'domain'] = 'DGR'

# %% 
'''Plot Domains Together'''

def create_plot_by_metric(df, metric, y, plot_out, inverse_metric=True):

    if y == 'prcnt_domain':
        ylabel = f'Domain Percentage in MAG'
    elif y == 'domain_count':
        ylabel = f'Domain Count in MAG'

    inverse_metric = False
    if ('mean' in metric) or ('median' in metric):
        inverse_metric = True
        df['diversity'] = 1-df[metric]
        xlabel = f'Diversity (1-{metric})'
    else:
        df['diversity'] = df[metric]
        xlabel = f'Diversity ({metric})'
        
    domains = ['DGR','CRISPR','ncRNA']
    domain_colors = ['green', 'red', 'purple']   

    df.sort_values(by='diversity', ascending=inverse_metric, inplace=True)

    fig, ax = plt.subplots(figsize=(14, 7))
    custom_palette = dict(zip(domains, domain_colors))
    sns.lineplot(data=df, y=y, x='diversity', 
                hue='domain', palette = custom_palette, 
                marker='o', alpha = 0.5, linewidth=2.5)

    ax.set_title(f'Domain Presence & Diversity in {basename}',fontsize=20)
    ax.set_xlabel(xlabel, fontsize = 15)
    ax.set_ylabel(ylabel, fontsize = 15)    

    # plt.savefig(f'{plot_out}.png', format="png", bbox_inches="tight")
    # plt.savefig(f'{plot_out}.pdf', format="pdf", bbox_inches="tight")
    plt.show()

# %% 
y = 'domain_count'

metric = 'mean_ad_norm'
plot_out = f'{ad_norm_cohort_stats_dir}/plots/{basename}_diversity_by_prcnt_domain_presence'
# plot_out = f'{output_dir}/{basename}.{in_extension}.{metric}.boxplot'
create_plot_by_metric(df, metric, y, plot_out, inverse_metric=True)

metric = 'stdev_ad_norm'
# plot_out = f'{output_dir}/{basename}.{in_extension}.{metric}.boxplot'
create_plot_by_metric(df, metric, y, plot_out, inverse_metric=False)

metric = 'prcnt_lessthan1'
# plot_out = f'{output_dir}/{basename}.{in_extension}.{metric}.boxplot'
create_plot_by_metric(df, metric, y, plot_out, inverse_metric=False)



#  %% 

def create_plot_by_metric_and_domain(df, domain, metric, plot_out, inverse_metric=True):
    
    df = df[df['domain']==domain]
    
    inverse_metric = False
    if ('mean' in metric) or ('median' in metric):
        inverse_metric = True
        df['diversity'] = 1-df[metric]
        ylabel = f'Diversity (1-{metric})'
    else:
        df['diversity'] = df[metric]
        ylabel = f'Diversity ({metric})'
        
    domains = ['DGR','CRISPR','ncRNA']
    domain_colors = ['green', 'red', 'purple']   

    df.sort_values(by='diversity', ascending=inverse_metric, inplace=True)

    fig, ax = plt.subplots(figsize=(14, 7))
    custom_palette = dict(zip(domains, domain_colors))
    sns.lineplot(data=df, x='prcnt_domain', y='diversity', 
                hue='domain', palette = custom_palette, 
                marker='o', alpha = 0.5, linewidth=2.5)

    ax.set_title(f'{domain} Presence & Diversity in {basename}',fontsize=20)
    ax.set_ylabel(ylabel, fontsize = 15)
    ax.set_xlabel(f'{domain} Percentage in MAG', fontsize = 15)    

    # plt.savefig(f'{plot_out}.png', format="png", bbox_inches="tight")
    # plt.savefig(f'{plot_out}.pdf', format="pdf", bbox_inches="tight")
    plt.show()


# %% 


domain = ['DGR','CRISPR','ncRNA']
metrics = ['mean_ad_norm','stdev_ad_norm','prcnt_lessthan1']

for domain in domains:
    for metric in metrics:
        plot_out = f'{ad_norm_cohort_stats_dir}/plots/{basename}_diversity_by_prcnt_domain_presence'
        # plot_out = f'{output_dir}/{basename}.{in_extension}.{metric}.boxplot'
        create_plot_by_metric_and_domain(df, domain, metric, plot_out, inverse_metric=True)




# %%

