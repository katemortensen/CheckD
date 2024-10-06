#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/17/2024
# Purpose: This script houses commonly used functions in the pipeline.

# %% 
import os
import statistics

# %% 
'''Add executable to PATH'''
# (not used anywhere, but keep in case we need it)
def add_exec_to_path(executable_name, executable_dir_path):
    if shutil.which(f'{executable_name}') == None:
        print(f'Adding {executable_name} dir to PATH.')
        os.environ["PATH"] += os.pathsep + executable_dir_path
        if shutil.which(f'{executable_name}') == None:
            print(f'{executable_name} not in PATH - check installation')
        else:
            print(f'{executable_name} successfully added')
    else:
        print(f'{executable_name} already in PATH.')


'''Check if file exists'''    
def already_done(path_to_check):
    result = False
    dir_to_check = "/".join(path_to_check.split('/')[0:-1])
    if os.path.exists(path_to_check):
        if os.path.getsize(path_to_check) != 0:
            result = True
    return result


'''Collect files'''
def list_files_with_extension(directory, extension):
    file_paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(extension):
                file_paths.append(os.path.join(root, file))
    return file_paths


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
                ad_norm = temp_dict.get('ad_norm')
                if result_dict.get(contig_name) == None:
                    result_dict[contig_name] = {pos:temp_dict}
                else:
                    result_dict[contig_name][pos] = temp_dict
        # The assert statements below cause issues in  ad_norm_checkpoint.py
        # if an entire contig is a coding domain (no nucleotides that are non-CDS)
        # if contig_count >= 1:
        #     assert len(result_dict[contig_name].keys()) > 1
        #     assert len(result_dict.keys()) == contig_count
    result_list = [] # ad_norms list
    if len(result_dict.keys()) > 0:
        if '_' in list(result_dict.keys())[0]: 
            ddupl_dict = {}
            for seq in result_dict.keys():
                contig = seq.split('_')[0]
                for pos in result_dict[seq].keys():
                    ad_norm = result_dict[seq][pos].get('ad_norm')
                    if ddupl_dict.get(contig) == None:
                        ddupl_dict[contig] = {pos:ad_norm}
                        result_list.append(ad_norm)
                    elif ddupl_dict[contig].get(pos) == None:
                        ddupl_dict[contig][pos]=ad_norm
                        result_list.append(ad_norm)
                    assert ad_norm == ddupl_dict[contig].get(pos)
        else:
            for contig in result_dict.keys():
                for pos in result_dict[contig].keys():
                    ad_norm = result_dict[contig][pos].get('ad_norm')
                    result_list.append(ad_norm)
    result_list = list(filter(lambda x: x not in ['N', 'NA', 'NaN', None], result_list))
    return result_dict, result_list

'''Compute mean, median, standard deviation, and percent less than 1'''
def mean_med_stdev(list_of_ad_norms):
    cleaned_list = list_of_ad_norms
    if len(cleaned_list) != 0:
        cleaned_list_num = [float(i) for i in cleaned_list]
        mean_ad_norm = round(statistics.mean(cleaned_list_num),4)
        median_ad_norm = round(statistics.median(cleaned_list_num),4)
        stdev_ad_norm = round(statistics.stdev(cleaned_list_num), 4)
        prcnt_lessthan1 = round((len([i for i in cleaned_list_num if i < 1])/len(list_of_ad_norms))*100,4)
    else:
        mean_ad_norm = 'NA'
        median_ad_norm = 'NA'
        stdev_ad_norm = 'NA'
        prcnt_lessthan1 = 'NA'
    result = [mean_ad_norm, median_ad_norm, stdev_ad_norm, prcnt_lessthan1]
    return(result)

'''Create general stats dictionary'''
'''Note: used in *_cohort_desc_stats.py scripts'''
def create_general_stats_dict(list_of_ad_norms):
    cleaned_list = list(filter(lambda x: x not in ['N','NA', 'NaN', None], list_of_ad_norms))
    assert len(cleaned_list) <= len(list_of_ad_norms)
    if len(cleaned_list) != 0:
        length = len(cleaned_list) 
        cleaned_list_num = [float(i) for i in cleaned_list]
        mean_ad_norm = round(statistics.mean(cleaned_list_num),4)
        median_ad_norm = round(statistics.median(cleaned_list_num),4)
        stdev_ad_norm = round(statistics.stdev(cleaned_list_num), 4)
        prcnt_lessthan1 = round((len([i for i in cleaned_list_num if i < 1])/len(list_of_ad_norms))*100,4)
    else:
        length = 'NA'
        mean_ad_norm = 'NA'
        median_ad_norm = 'NA'
        stdev_ad_norm = 'NA'
        prcnt_lessthan1 = 'NA'
    general_stats_dict = {'length':length, 
                          'mean':mean_ad_norm,
                          'median':median_ad_norm,
                          'stdev':stdev_ad_norm,
                          'prcnt_lessthan1':prcnt_lessthan1}
    return(general_stats_dict)


'''create_desc_ad_norm_dict'''
'''desc > ad_norms list, count'''
'''Note: used in *_cohort_desc_stats.py scripts'''

def create_desc_ad_norm_dict(ad_norm_file):
    result_dict = {}
    contig_count = 0
    with open(ad_norm_file,'r') as fp_in:
        contig_count = 0
        for line in fp_in:
            if line.startswith('# cols:'):
                info = line.strip().split('# cols: ')[-1].split(', ')
            elif line.startswith('>'):
                desc = line.strip().split(';')[-1]
                if result_dict.get(desc) == None:
                    result_dict[desc] = {'ad_norms':[], 'count':1}
                else:
                    new_count = 1 + int(result_dict[desc].get('count'))
                    result_dict[desc]['count'] = new_count
                contig_count = 1
            elif contig_count == 1:
                linep = line.strip('\n').split('\t')                
                temp_dict = dict(zip(info,linep))
                ad_norm = temp_dict.get('ad_norm')      
                old_ad_norms = result_dict[desc].get('ad_norms')
                new_ad_norms = [ad_norm] + old_ad_norms
                result_dict[desc]['ad_norms'] = new_ad_norms                
    result_dict = dict(sorted(result_dict.items(), key=lambda item: item[1]['count'], reverse=False))
    return result_dict


'''Compile desc dictionaries across all mags'''
'''Note: used in *_cohort_desc_stats.py scripts'''

def compile_desc_dicts_across_mags(list_of_ad_norm_files):
    result_dict = {}
    for ad_norm_file in list_of_ad_norm_files:
        temp_dict = create_desc_ad_norm_dict(ad_norm_file)
        if len(temp_dict) > 0:
            for desc in temp_dict.keys():
                if result_dict.get(desc) == None:
                    result_dict[desc] = temp_dict.get(desc)
                else:
                    ad_norms = temp_dict[desc].get('ad_norms')      
                    old_ad_norms = result_dict[desc].get('ad_norms')
                    new_ad_norms = ad_norms + old_ad_norms
                    count = temp_dict[desc].get('count')
                    old_count = result_dict[desc].get('count')
                    new_count = int(count) + int(old_count)
                    result_dict[desc]['ad_norms'] = new_ad_norms
                    result_dict[desc]['count'] = new_count
        else:
            print(f'No positions found in: {ad_norm_file}')             
    if len(result_dict) > 0:
        result_dict = dict(sorted(result_dict.items(), key=lambda item: item[1]['count'], reverse=True))
    return result_dict   


'''Create desc dictionary with general stats and sort by mean_ad_norm'''
'''desc > count, length, ad_norm_list, mean, median, st_dev, prcnt_lessthan1'''
'''Note: used in *_cohort_desc_stats.py scripts'''

def create_desc_stats_dict(list_of_ad_norm_files):
    descDict = compile_desc_dicts_across_mags(list_of_ad_norm_files)        
    for desc in descDict.keys():
        ad_norms = descDict[desc].get('ad_norms')
        general_stats_dict = create_general_stats_dict(ad_norms)
        descDict[desc].update(general_stats_dict)
    if len(descDict) > 0:
        descDict = dict(sorted(descDict.items(), key=lambda item: item[1]['mean'], reverse=False))

    total_nucs = 0
    for desc in descDict.keys():
        new_length = descDict[desc].get('length')
        total_nucs += int(new_length)

    # cross check nucleotide count
    total_nucs_check = 0
    for ad_norm_file in list_of_ad_norm_files:
        with open(ad_norm_file,'r') as fp_in:
            contig_count = 0
            for line in fp_in:
                if line.startswith('>'):
                    desc = line.strip().split(';')[-1]
                    contig_count += 1
                elif (line.startswith('>') == False) and (contig_count >= 1):
                    total_nucs_check += 1
    # assert total_nucs == total_nucs_check
    # NOTE: the above statemnet failes when there are 'N' nucs at positions or
    # when there are 'NA' or 'NaN' for the ad_norms of ATCG nucs
    return(descDict)


'''Sort list of mags'''
def sort_list_by_mag(mag_list):
    def extract_number(s):
        return int(s.split('.')[1])
    # Sorting the list by the numeric part of the strings in ascending order
    sorted_list = sorted(mag_list, key=extract_number)
    return(sorted_list)


