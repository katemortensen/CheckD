#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/24/2024
# Purpose: Assess the diversity across MAGs


# %% 
import sys
import os
import subprocess
import time
from time import gmtime, strftime
import argparse
import statistics
import math
import pandas as pd
import numpy as np
import shutil
import pysam
import matplotlib
import checkm
import concurrent.futures
import multiprocessing
from checkd_functions import already_done


# '''NOTE: mag names cannot have underscores - compromises parsing scripts.'''

# '''NOTE: contig names cannot have underscores - compromises parsing scripts.'''

# '''NOTE: mag fasta files must end in ".fa" '''

# '''NOTE: mag fasta files can be the only files in the mag directory and 
# their fasta extensions should be the same (ex: .fa)'''

     
# %% 

# Initialize the argument parser
parser = argparse.ArgumentParser(description="Run the pipeline or specific steps")

# Add a subparser for different sections of the pipeline
subparsers = parser.add_subparsers(dest="command", help="Pipeline sections")

# Full pipeline argument
parser_all = subparsers.add_parser("all", help="Run the entire pipeline")
parser_all.add_argument("-basename", "--basename", metavar="", help="", required=True)
parser_all.add_argument("-output_dir", "--output_dir", metavar="", help="", required=True)
parser_all.add_argument("-mag_dir", "--mag_dir", metavar="", help="", required=True)
parser_all.add_argument("-hifi_reads", "--hifi_reads", metavar="hifi reads", help="", required=False)
parser_all.add_argument("-illumina_forward_reads", "--illumina_forward_reads", metavar="illumina reads", help="", required=False)
parser_all.add_argument("-illumina_reverse_reads", "--illumina_reverse_reads", metavar="illumina reads", help="", required=False)
parser_all.add_argument("-window", "--window", metavar="", help="", required=False)
parser_all.add_argument("-step", "--step", metavar="", help="", required=False)
parser_all.add_argument("-threads", "--threads", metavar="", help="", required=False)

# Individual pipeline steps
parser1 = subparsers.add_parser("map", help="Map Reads to MAGs")
parser1.add_argument("-basename", "--basename", metavar="", help="", required=True)
parser1.add_argument("-output_dir", "--output_dir", metavar="", help="", required=True)
parser1.add_argument("-mag_dir", "--mag_dir", metavar="", help="", required=True)
parser1.add_argument("-hifi_reads", "--hifi_reads", metavar="hifi reads", help="", required=False)
parser1.add_argument("-illumina_forward_reads", "--illumina_forward_reads", metavar="illumina reads", help="", required=False)
parser1.add_argument("-illumina_reverse_reads", "--illumina_reverse_reads", metavar="illumina reads", help="", required=False)
parser1.add_argument("-window", "--window", metavar="", help="", required=False)
parser1.add_argument("-step", "--step", metavar="", help="", required=False)
parser1.add_argument("-threads", "--threads", metavar="", help="", required=False)

parser2 = subparsers.add_parser("region", help="Region Detection")
parser2.add_argument("-basename", "--basename", metavar="", help="", required=True)
parser2.add_argument("-output_dir", "--output_dir", metavar="", help="", required=True)
parser2.add_argument("-mag_dir", "--mag_dir", metavar="", help="", required=True)
parser2.add_argument("-hifi_reads", "--hifi_reads", metavar="hifi reads", help="", required=False)
parser2.add_argument("-illumina_forward_reads", "--illumina_forward_reads", metavar="illumina reads", help="", required=False)
parser2.add_argument("-illumina_reverse_reads", "--illumina_reverse_reads", metavar="illumina reads", help="", required=False)
parser2.add_argument("-threads", "--threads", metavar="", help="", required=False)

parser3 = subparsers.add_parser("stats", help="Calculate Diversity Stats")
parser3.add_argument("-basename", "--basename", metavar="", help="", required=True)
parser3.add_argument("-output_dir", "--output_dir", metavar="", help="", required=True)
parser3.add_argument("-mag_dir", "--mag_dir", metavar="", help="", required=True)
parser3.add_argument("-hifi_reads", "--hifi_reads", metavar="hifi reads", help="", required=False)
parser3.add_argument("-illumina_forward_reads", "--illumina_forward_reads", metavar="illumina reads", help="", required=False)
parser3.add_argument("-illumina_reverse_reads", "--illumina_reverse_reads", metavar="illumina reads", help="", required=False)
parser3.add_argument("-threads", "--threads", metavar="", help="", required=False)

parser4 = subparsers.add_parser("sliding_window_plot", help="Plot Diversity and Regions")
parser4.add_argument("-basename", "--basename", metavar="", help="", required=True)
parser4.add_argument("-output_dir", "--output_dir", metavar="", help="", required=True)
parser4.add_argument("-mag_dir", "--mag_dir", metavar="", help="", required=True)
parser4.add_argument("-hifi_reads", "--hifi_reads", metavar="hifi reads", help="", required=False)
parser4.add_argument("-illumina_forward_reads", "--illumina_forward_reads", metavar="illumina reads", help="", required=False)
parser4.add_argument("-illumina_reverse_reads", "--illumina_reverse_reads", metavar="illumina reads", help="", required=False)
parser4.add_argument("-threads", "--threads", metavar="", help="", required=False)

args = parser.parse_args()
command = args.command
basename = args.basename
output_dir = args.output_dir
input_mag_dir= args.mag_dir
hifi_reads = args.hifi_reads
forward_reads = args.illumina_forward_reads
reverse_reads = args.illumina_reverse_reads
wkdir = f'{output_dir}/{basename}_checkd'


if not os.path.exists(f'{output_dir}'):
    os.makedirs(f'{output_dir}')
    
if not os.path.exists(f'{wkdir}'):
    os.makedirs(f'{wkdir}')

# Redirect stdout to a log file
log_file = open(f'{wkdir}/checkd_log.txt', 'w')
sys.stdout = log_file

if (forward_reads != None) and (reverse_reads != None):
    read_type = 'illumina'
    assert hifi_reads == None
if hifi_reads != None:
    read_type = 'hifi'
assert read_type != None

if args.threads != None:
    threads = int(args.threads)
else:
    threads = 1

if args.window != None:
    window = args.window
else:
    window = None

if args.step != None:
    step = args.step
else:
    step = None


# COMMENT THIS LINE OUT BEFORE RUNNING FROM COMMAND LINE

# %% 

'''Pre-sets'''
# scripts = os.getcwd()
scripts = os.path.dirname(os.path.abspath(__file__))
checkd_bin = '%s/bin' % (scripts.strip('scripts'))
print(f'bin path: {checkd_bin}')

minimap_mag_dir = f'{wkdir}/reads2mags_minimap'
fgs_out_dir = f'{wkdir}/fgs_out'
results_dir = f'{wkdir}/results'
checkm_out_dir = f'{wkdir}/checkm_out'
rfam_out_dir = f'{wkdir}/rfam_out'
dgrscan_out_dir = f'{wkdir}/dgrscan_out'
crisprone_out_dir = f'{wkdir}/crisprone_out'
virsorter_out_dir = f'{wkdir}/virsorter2_out'
ad_norm_dir = f'{wkdir}/ad_norm'
ad_norm_cohort_stats_dir = f'{wkdir}/ad_norm_cohort_stats'
ad_norm_mag_stats_dir = f'{wkdir}/ad_norm_mag_stats'
sliding_window_dir = f'{wkdir}/sliding_window'


'''Software'''
python = 'python3'
minimap2 = f'{checkd_bin}/minimap2/minimap2'
samtools = f'{checkd_bin}/samtools-1.20/samtools'
bcftools = f'{checkd_bin}/bcftools-1.20/bcftools'
fgs = f'{checkd_bin}/CRISPRone-main/bin/FragGeneScan1.31/run_FragGeneScan.pl'
cmscan = f'{checkd_bin}/infernal-1.1.2/bin/cmscan'
esl_seqstat = f'{checkd_bin}/infernal-1.1.2/bin/esl-seqstat'
Rfam_cm = f'{checkd_bin}/Rfam.cm'
Rfam_clanin = f'{checkd_bin}/Rfam.clanin'
crisprone = f'{checkd_bin}/CRISPRone-main/crisprone-local-20240911_KM.py'
summarize = f'{checkd_bin}/metaCRT/summarize-crispr-new.py'
mydgr = f'{checkd_bin}/myDGR_cp/mydgr-local-acc-RT.py'


# %% 

print(
'\n#####################################'
'\nSetup'
'\n#####################################'
)

print('''\nRename assemblies and contigs. Prepare for run.''')

script = f'{scripts}/rename_assemblies_and_contigs.py'
cmd = f'python3 {script} \
    --mag_dir {input_mag_dir} \
    --output_dir {wkdir}'
try:
    subprocess.check_output(cmd, shell=True)
except subprocess.CalledProcessError as e:
    print(e.output)
    
mag_dir = f'{wkdir}/renamed_mags'
path_to_mags = f'{wkdir}/paths_to_mags.txt'

mags = os.listdir(mag_dir)
with open(path_to_mags, "w") as fp:
    for m in mags:
        path = f'{mag_dir}/{m}\n'
        fp.write(path)

mag_extensions = []
with open(path_to_mags, 'r') as file:
    pmags = [line.strip() for line in file.readlines()]
    for pmag in pmags:
        extension = pmag.strip().split('.')[-1]
        mag_extensions.append(extension)
assert len(list(set(mag_extensions))) == 1

mag_extension = f'.{list(set(mag_extensions))[0]}'
   

print('''\nSymlink MAGs''')

if not os.path.exists(f'{minimap_mag_dir}' ):
    os.makedirs(f'{minimap_mag_dir}' )

with open(path_to_mags, 'r') as file:
    pmags = [line.strip() for line in file.readlines()]
    for pmag in pmags:
        mag_name = pmag.split('/')[-1].split(mag_extension)[0]
        cmd = f'''ln -s {pmag} {minimap_mag_dir}/{mag_name}.fa'''
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            print(e.output)

# %% 

if (command == 'map') or (command == 'all'):

    print(
    '\n#####################################'
    '\nMap Reads to MAGs'
    '\n#####################################'
    )

    print('''\nIndex MAG Fasta Files''')

    with open(path_to_mags, 'r') as file:
        pmags = [line.strip() for line in file.readlines()]
        for pmag in pmags:
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            if already_done(f'{minimap_mag_dir}/{mag_name}.fa.fai') is False:
                cmd = f'''{samtools} faidx \
                    {minimap_mag_dir}/{mag_name}.fa'''
                subprocess.check_output(cmd, shell=True)
                
                cmd = f'''chmod 775 {minimap_mag_dir}/{mag_name}.fa.fai'''
                subprocess.check_output(cmd, shell=True)
            else:
                print(f'{mag_name} fasta already indexed')
            

    '''Map Reads to MAGs - Illumina'''

    # input files: forward & reverse reads, reference MAG
    # output files: sorted bam file, bam index file, and coverage

    # ORIGINAL - works (for Illumina reads) 

    if read_type == 'illumina':

        print('''\nMap Reads to MAGs - Illumina''')

        
        if not os.path.exists(f'{minimap_mag_dir}' ):
            os.makedirs(f'{minimap_mag_dir}' )

        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                
                if already_done(f'{minimap_mag_dir}/{mag_name}_sorted.bam') is False:
                    cmd = f'''{minimap2} -ax sr \
                        {pmag} {forward_reads} {reverse_reads} \
                        | {samtools} sort --threads {threads} \
                        -T {minimap_mag_dir}/{mag_name} \
                        -o {minimap_mag_dir}/{mag_name}_sorted.bam - '''
                    subprocess.check_output(cmd, shell=True)
                else:
                    print(f'{mag_name} already mapped')
                    
                if already_done(f'{minimap_mag_dir}/{mag_name}_sorted.bai') is False:
                    cmd = f'''{samtools} index \
                        {minimap_mag_dir}/{mag_name}_sorted.bam \
                        {minimap_mag_dir}/{mag_name}_sorted.bai'''
                    subprocess.check_output(cmd, shell=True) 
                else:
                    print(f'{mag_name} bam already indexed')

                # '''Coverage'''
                # cmd = f'''{samtools} depth \
                #     {minimap_mag_dir}/{mag_name}_sorted.bam \
                #     > {minimap_mag_dir}/{mag_name}.cov'''
                # subprocess.check_output(cmd, shell=True)    


    '''Map Reads to MAGs - HiFi'''

    # input files: reads fastq file, reference MAG
    # output files: sorted bam file, bam index file, and coverage

    if read_type == 'hifi':

        print('''\nMap Reads to MAGs - HiFi''')

        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]

                if already_done(f'{minimap_mag_dir}/{mag_name}.mmi') is False:
                    # minimap index 
                    cmd = f'''{minimap2} -t {threads} -d \
                        {minimap_mag_dir}/{mag_name}.mmi {pmag} '''
                    subprocess.check_output(cmd, shell=True) 
                    
                if already_done(f'{minimap_mag_dir}/{mag_name}_sorted.bam') is False:
                    # minimap reads to mag (sam/bam)
                    cmd = f'''{minimap2} -t {threads} -ax map-hifi \
                        {pmag} {hifi_reads} \
                        | {samtools} sort --threads {threads} \
                        -T {minimap_mag_dir}/{mag_name} \
                        -o {minimap_mag_dir}/{mag_name}_sorted.bam - '''
                    subprocess.check_output(cmd, shell=True) 
                else:
                    print(f'{mag_name} already mapped')
                                                
                if already_done(f'{minimap_mag_dir}/{mag_name}_sorted.bai') is False:
                    # index bam
                    cmd = f'''{samtools} index \
                        {minimap_mag_dir}/{mag_name}_sorted.bam \
                        {minimap_mag_dir}/{mag_name}_sorted.bai'''
                    subprocess.check_output(cmd, shell=True) 
                else:
                    print(f'{mag_name} bam already indexed')
                    
                # '''Coverage'''
                # cmd = f'''{samtools} depth \
                #     {minimap_mag_dir}/{mag_name}.bam \
                #     > {minimap_mag_dir}/{mag_name}_sorted.cov'''
                # subprocess.check_output(cmd, shell=True) 

        

    print('''\nbcftools mpileup - VCF file''')

    # samtools view bin.41_SRR15275213.sorted.bam | cut -f1 | sort | uniq > /home/kmorten/CheckD/example/hifi_unique_reads_list.txt
    # sort input.txt | uniq > output.txt
    # grep -A 3 -Ff reads_list.txt input.fastq > output_subset.fastq

    # input files:
    # <minimap_mag_dir>/<mag_name>_sorted.bam
    # <mag_name>.fa
    # output file: 
    # <subwkdir>/pos_req/<mag_name>.ref_pos_freq

    with open(path_to_mags, 'r') as file:
        pmags = [line.strip() for line in file.readlines()]
        for pmag in pmags:
            mag_name = pmag.split('/')[-1].split(f'{mag_extension}')[0]
            if already_done(f'{minimap_mag_dir}/{mag_name}.vcf') is False:
                cmd = f'''{bcftools} mpileup --threads {threads} \
                    -Ov -f {minimap_mag_dir}/{mag_name}{mag_extension} \
                    {minimap_mag_dir}/{mag_name}_sorted.bam \
                    -O 'v' -o {minimap_mag_dir}/{mag_name}.vcf \
                    -a 'FORMAT/AD','FORMAT/DP','FORMAT/QS','INFO/AD','INFO/NM' '''
                subprocess.check_output(cmd, shell=True)    
            else:
                print(f'VCF for {mag_name} already exists')

# %% 

if (command == 'region') or (command == 'all'): 

    print(
    '\n#####################################'
    '\nRegion Detection'
    '\n#####################################'
    )
                

    print('''\nFragGeneScan: Codon Positions of Putative Genes ''')

    # input: MAG fasta file 
    # output file: <subwkdir>/fgs_out

    if not os.path.exists(f'{fgs_out_dir}' ):
        os.makedirs(f'{fgs_out_dir}' )
        
    if read_type == 'illumina':
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:   
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(mag_name)            
                if already_done(f'{fgs_out_dir}/{mag_name}-fgs') is False:
                    cmd = f'''{fgs} -genome={pmag} \
                        -out={fgs_out_dir}/{mag_name}-fgs \
                        -complete=0  -train=illumina_5'''
                    subprocess.check_output(cmd, shell=True) 

    if read_type == 'hifi':
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(mag_name)
                if already_done(f'{fgs_out_dir}/{mag_name}-fgs') is False:
                    cmd = f'''{fgs} -genome={pmag} \
                        -out={fgs_out_dir}/{mag_name}-fgs \
                        -complete=1  -train=complete'''
                    subprocess.check_output(cmd, shell=True) 

    print('''\nrfam.py (ncRNA)''')

    # input files:
    # fas_file = path/to/<mag_name>.fa

    # output files:
    # seqstat_file_out = ../rfam_out/<mag_name>.seqstat
    # cmscan_file_out = ../rfam_out/<mag_name>.cmscan
    # tblout_file_out = ../rfam_out/<mag_name>.tblout
    # tblout_deovlp_file_out = ../rfam_out/<mag_name>.deoverlapped.tblout

    with open(path_to_mags, 'r') as file:
        pmags = [line.strip() for line in file.readlines()]
        for pmag in pmags:
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(mag_name)
            if already_done(f'{rfam_out_dir}/{mag_name}.tblout') is False:
                cmd = f'''{python} {scripts}/rfam.py \
                    --fas_file {pmag} \
                    --output_dir {rfam_out_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 


    print('''\nmyDGRscan''')

    # input files:
    # path/to/<mag_name>.fa
    # output dir:
    # ../wkdir/dgrscan_out

    if not os.path.exists(dgrscan_out_dir):
        os.makedirs(dgrscan_out_dir)
        
    with open(path_to_mags, 'r') as file:
        pmags = [line.strip() for line in file.readlines()]
        for pmag in pmags:
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(mag_name)
            cmd = f'''{mydgr} {dgrscan_out_dir} {mag_name} {pmag}'''        
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print(e.output)


    print('''\nCRISPRone''')

    with open(path_to_mags, 'r') as file:
        pmags = [line.strip() for line in file.readlines()]
        for pmag in pmags:
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(mag_name)
            'CRISPRone'
            try:
                subprocess.run([
                    f'{python}', f'{crisprone}',
                    '--outbase', f'{crisprone_out_dir}',
                    '--outprefix', f'{mag_name}',
                    '--fna', f'{pmag}',
                    '--no_tracrRNA', 'True'
                ], check=True, capture_output=True, text=True)
                print(f'{crisprone} - done')
            except subprocess.CalledProcessError as e:
                print(f'Error: {crisprone}')
                print(f"Command failed with error: {e.stderr}")


    print('''\nCRISPRone Summary''')

    with open(path_to_mags, 'r') as file:
        pmags = [line.strip() for line in file.readlines()]
        for pmag in pmags:
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(mag_name)        
            'CRISPRone Summary'
            # cmd = f'''{python} {summarize} \
            #     -f {crisprone_out_dir}/{mag_name}/*.crt \
            #     -repeat {crisprone_out_dir}/{mag_name}/repeat.seq \
            #     -consensus {crisprone_out_dir}/{mag_name}/repeatcon.seq'''
            # subprocess.check_output(cmd, shell=True) 
            
            try:
                subprocess.run([
                    f'{python}', f'{summarize}', 
                    '-f', f'{crisprone_out_dir}/{mag_name}/{mag_name}.crt',
                    '-repeat', f'{crisprone_out_dir}/{mag_name}/repeat.seq',
                    '-consensus', f'{crisprone_out_dir}/{mag_name}/repeatcon.seq'],
                    check=True, capture_output=True, text=True)    
            except subprocess.CalledProcessError as e:
                print(f'Error: {summarize}')
                print(f"Command failed with error: {e.stderr}")
            

    print('''\nVirSorter2''')

    # NOTE: This step takes a long time

    script = f'virsorter2_singularity.py'
    with open(path_to_mags, 'r') as file:
        pmags = [line.strip() for line in file.readlines()]
        for pmag in pmags:
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(mag_name)
            if already_done(f'{virsorter_out_dir}/{mag_name}/final-viral-boundary.tsv') is False:
                try:
                    subprocess.run([
                        f'{python}', f'{scripts}/{script}',
                        '--threads', f'{threads}',
                        '--assembly_name', f'{mag_name}',
                        '--fasta_file', f'{pmag}',
                        '--output_dir', f'{virsorter_out_dir}'
                    ], check=True, capture_output=True, text=True)
                    print(f'{script} - done')
                except subprocess.CalledProcessError as e:
                    print(f'Error: {script}')
                    print(f"Command failed with error: {e.stderr}")
            else:
                print(f'VirSorter2 already done for {mag_name}')


# %%

if (command == 'stats') or (command == 'all'):

    print(
    '\n#####################################'
    '\nStats'
    '\n#####################################'
    )


    print('''\nAllelic Depth - All Positions''')

    # input files:
    # vcf_file = f'/home/kmorten/humanO1_CheckD/reads2mags_minimap/{mag_name}.vcf'
    # output_dir = f'/home/kmorten/humanO1_CheckD/allelic_depth'
    # output_file = f'{output_dir}/{mag_name}.ad_norm'

    if threads > 2: 
        print('\nrunning multiprocessing tool for ad_norm.py')
        # Define the function to be used with multiprocessing
        def ad_norm_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/ad_norm.py \
                --vcf_file {minimap_mag_dir}/{mag_name}.vcf \
                --basename {mag_name} \
                --output_dir {ad_norm_dir}'''
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print(output.decode())
            except subprocess.CalledProcessError:
                print(f'Error: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(ad_norm_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()
    else:
        print('\nrunning without multiprocessing tool')
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(mag_name)
                cmd = f'''{python} {scripts}/ad_norm.py \
                    --vcf_file {minimap_mag_dir}/{mag_name}.vcf \
                    --basename {mag_name} \
                    --output_dir {ad_norm_dir}'''
                subprocess.check_output(cmd, shell=True) 


    print('''\nAllelic Depth - CDS Regions''')

    # input files:
    # fgs_file =  f'/home/kmorten/humanO1_CheckD/fgs_out/{mag_name}-fgs.out'
    # ffn_file = f'/home/kmorten/humanO1_CheckD/fgs_out/{mag_name}-fgs.ffn'
    # ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'

    # output files:
    # cds_output_file = f'{output_dir}/{mag_name}.cds_ad_norm'

    if threads > 2: 
        print('\nrunning multiprocessing tool for cds_ad_norm.py')
        def cds_ad_norm_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/cds_ad_norm.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                --fgs_file {fgs_out_dir}/{mag_name}-fgs.out \
                --ffn_file {fgs_out_dir}/{mag_name}-fgs.ffn \
                --output_dir {ad_norm_dir} \
                --assembly_name {mag_name}'''
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(cds_ad_norm_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                cmd = f'''{python} {scripts}/cds_ad_norm.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                    --fgs_file {fgs_out_dir}/{mag_name}-fgs.out \
                    --ffn_file {fgs_out_dir}/{mag_name}-fgs.ffn \
                    --output_dir {ad_norm_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 


    print('''Allelic Depth - non-CDS Regions''')

    # input files:
    # fgs_file =  f'/home/kmorten/humanO1_CheckD/fgs_out/{mag_name}-fgs.out'
    # ad_norm_file = f'{output_dir}/{mag_name}.ad_norm'
    # cds_ad_norm_file = f'{output_dir}/{mag_name}.cds_ad_norm'

    # output files:
    # noncds_output_file = f'{output_dir}/{mag_name}.noncds_ad_norm'

    if threads > 2: 
        print('running multiprocessing tool for noncds_ad_norm.py')
        def cds_ad_norm_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/noncds_ad_norm.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                --cds_ad_norm_file {ad_norm_dir}/{mag_name}.cds_ad_norm \
                --fgs_file {fgs_out_dir}/{mag_name}-fgs.out \
                --output_dir {ad_norm_dir} \
                --assembly_name {mag_name}'''        
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(cds_ad_norm_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                cmd = f'''{python} {scripts}/noncds_ad_norm.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                    --cds_ad_norm_file {ad_norm_dir}/{mag_name}.cds_ad_norm \
                    --fgs_file {fgs_out_dir}/{mag_name}-fgs.out \
                    --output_dir {ad_norm_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 


    # print('''Allelic Depth -- Checkpoint''')
    # '''ad_norm_checkpoint.py'''

    # # NOTE: reconciles CDS and non-CDS nucleotide counts 
    # # NOTE: fails when VCF file and FragGeneScan identify different INDELS

    # input files:
    # ...ad_norm/(all files in dir)

    # try:
    #     subprocess.run([
    #         f'{python}', f'{scripts}/ad_norm_checkpoint.py',
    #         '--ad_norm_dir', f'{ad_norm_dir}'
    #     ], check=True, capture_output=True, text=True)
    # except subprocess.CalledProcessError as e:
    #     print(f"Command failed with error: {e.stderr}")
        

    print('''\nAllelic Depth - CRISPR Regions''')

    if threads > 2: 
        print('\nrunning multiprocessing tool for crispr_ad_norm.py')
        def crispr_ad_norm_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/crispr_ad_norm.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                --crt_file {crisprone_out_dir}/{mag_name}/{mag_name}.crt \
                --crisprone_gff_file {crisprone_out_dir}/{mag_name}/{mag_name}-sm.gff \
                --output_dir {ad_norm_dir} \
                --assembly_name {mag_name}'''
            subprocess.check_output(cmd, shell=True)          
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(crispr_ad_norm_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()           
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                cmd = f'''{python} {scripts}/crispr_ad_norm.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                    --crt_file {crisprone_out_dir}/{mag_name}/{mag_name}.crt \
                    --crisprone_gff_file {crisprone_out_dir}/{mag_name}/{mag_name}-sm.gff \
                    --output_dir {ad_norm_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True)        


    print('''\nAllelic Depth - ncRNA Regions''')

    if threads > 2: 
        print('\nrunning multiprocessing tool for ncrna_ad_norm.py')
        def ncrna_ad_norm_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/ncrna_ad_norm.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                --tblout_deovlp_file {rfam_out_dir}/{mag_name}.deoverlapped.tblout \
                --output_dir {ad_norm_dir} \
                --assembly_name {mag_name}'''
            subprocess.check_output(cmd, shell=True)                    
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(ncrna_ad_norm_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()    
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                cmd = f'''{python} {scripts}/ncrna_ad_norm.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                    --tblout_deovlp_file {rfam_out_dir}/{mag_name}.deoverlapped.tblout \
                    --output_dir {ad_norm_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 
                


    print('''\nAllelic Depth - DGR ''')

    if threads > 2:  
        print('running multiprocessing tool for dgr_ad_norm.py')
        def dgr_ad_norm_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/dgr_ad_norm.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                --dgr_gff_file {dgrscan_out_dir}/{mag_name}/{mag_name}-DGRscan.gff \
                --output_dir {ad_norm_dir} \
                --assembly_name {mag_name}'''
            subprocess.check_output(cmd, shell=True)                    
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(dgr_ad_norm_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()    
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                cmd = f'''{python} {scripts}/dgr_ad_norm.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.ad_norm \
                    --dgr_gff_file {dgrscan_out_dir}/{mag_name}/{mag_name}-DGRscan.gff \
                    --output_dir {ad_norm_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 

            
    print('''\nAllelic Depth - VIR ''')

    script = 'vir_ad_norm.py'
    if threads > 2: 
        print(f'\nrunning multiprocessing tool for {script}')
        def vir_ad_norm_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            try:
                subprocess.run([
                    f'{python}', f'{scripts}/{script}',
                    '--assembly_name', f'{mag_name}',
                    '--ad_norm_file', f'{ad_norm_dir}/{mag_name}.ad_norm',
                    '--final_viral_boundary_file', f'{virsorter_out_dir}/{mag_name}/final-viral-boundary.tsv',
                    '--final_viral_combined_file', f'{virsorter_out_dir}/{mag_name}/final-viral-combined.fa',
                    '--output_dir', f'{ad_norm_dir}'
                ], check=True, capture_output=True, text=True)
                # print(f'{script} - {mag_name} - done')
            except subprocess.CalledProcessError as e:
                print(f'\nError: {mag_name}; {script}')
                print(f"\nCommand failed with error: {e.stderr}")
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(vir_ad_norm_func, pmags) 
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()    
        print(f'\n{script} - done')
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                try:
                    subprocess.run([
                        f'{python}', f'{scripts}/{script}',
                        '--assembly_name', f'{mag_name}',
                        '--ad_norm_file', f'{ad_norm_dir}/{mag_name}.ad_norm',
                        '--final_viral_boundary_file', f'{virsorter_out_dir}/{mag_name}/final-viral-boundary.tsv',
                        '--final_viral_combined_file', f'{virsorter_out_dir}/{mag_name}/final-viral-combined.fa',
                        '--output_dir', f'{ad_norm_dir}'
                    ], check=True, capture_output=True, text=True)
                    # print(f'{script} - done')
                except subprocess.CalledProcessError as e:
                    print(f'\nError: {mag_name}; {script}')
                    print(f"\nCommand failed with error: {e.stderr}")
        print(f'\n{script} - done')



    print('''\nDomain -- Allelic Depth -- Cohort Stats''')
    '''cohort_domain_stats.py'''

    # input files:
    # ...ad_norm/(all files in dir)

    # output files:
    # ...ad_norm_cohort_stats/<basename>.cohort_domain_stats

    try:
        subprocess.run([
            f'{python}', f'{scripts}/cohort_domain_stats.py',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"\nCommand failed with error: {e.stderr}")
        


    print('''\nAll Positions - Allelic Depth - Cohort MAG Stats ''')
        
    script = 'cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")

    script = 'cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


    print('''\nCDS Positions - Allelic Depth - Cohort MAG Stats ''')

    # input files:
    # ...ad_norm/<mag_name>.ad_norm
    # ...ad_norm/<mag_name>.cds_ad_norm

    # output files:
    # ...ad_norm_cohort_stats/<basename>.cds_cohort_mag_stats

    script = 'cds_cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


    print('''\nnon-CDS positions - Allelic Depth - Cohort MAG Stats''')

    # input files:
    # ...ad_norm/<mag_name>.ad_norm
    # ...ad_norm/<mag_name>.noncds_ad_norm

    # output files:
    # ...ad_norm_cohort_stats/<basename>.noncds_cohort_mag_stats

    script = 'noncds_cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


  
    print('''\nCRISPR - Allelic Depth - Cohort MAG Stats''')

    # input files:
    # ...ad_norm/<mag_name>.ad_norm
    # ...ad_norm/<mag_name>.crispr_ad_norm

    # output files:
    # ...ad_norm_cohort_stats/<basename>.crispr_cohort_mag_stats

    script = 'crispr_cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


    print('''\nncRNA - Allelic Depth - Cohort MAG Stats''')

    # input files:
    # ...ad_norm/<mag_name>.ad_norm
    # ...ad_norm/<mag_name>.ncrna_ad_norm

    # output files:
    # ...ad_norm_cohort_stats/<basename>.ncrna_cohort_mag_stats

    script = 'ncrna_cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")
        
        
    print('''\nDGR - Allelic Depth - Cohort MAG Stats''')

    # input files:
    # ...ad_norm/<mag_name>.ad_norm
    # ...ad_norm/<mag_name>.dgr_ad_norm

    # output files:
    # ...ad_norm_cohort_stats/<basename>.dgr_cohort_mag_stats

    script ='dgr_cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


    print('''\nVIR - Allelic Depth - Cohort MAG Stats''')

    script ='vir_cohort_mag_stats.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


    print('''\nCohort MAG Diversity Baseline''')

    # input files:
    # ...ad_norm_cohort_stats/<basename>.mag_stats
    # ...ad_norm_cohort_stats/<basename>.<domain>_mag_stats

    # output files:
    # ...ad_norm_cohort_stats/<basename>.cohort_mag_diversity_baseline

    # NOTE: This script also contains UMAP and correlation analyses but 
    # does not save these plots (not useful).

    script ='cohort_mag_diversity_baseline.py'
    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_cohort_stats_dir', f'{ad_norm_cohort_stats_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}',
            '--metric', f'mean_ad_norm'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


    print('''\nCohort Desc Stats Section''')

    # Note: This section does not inlcude CDS or non-CDS 

    print('''\nCRISPR - Allelic Depth - Cohort Desc Stats''')
    script = f'crispr_cohort_desc_stats.py'

    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")


    print('''\nncRNA - Allelic Depth - Cohort Desc Stats''')
    script = f'ncrna_cohort_desc_stats.py'

    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}")

    print('''\nDGR - Allelic Depth - Cohort Desc Stats''')
    script = f'dgr_cohort_desc_stats.py'

    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}") 


    print('''\nVIR - Allelic Depth - Cohort Desc Stats''')
    script = 'vir_cohort_desc_stats.py'

    try:
        subprocess.run([
            f'{python}', f'{scripts}/{script}',
            '--ad_norm_dir', f'{ad_norm_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
        print(f'\n{script} - done')
    except subprocess.CalledProcessError as e:
        print(f'\nError: {script}')
        print(f"\nCommand failed with error: {e.stderr}") 


    print('''\nCompiled - Allelic Depth - Cohort Desc Stats''')
    '''cohort_desc_stats.py'''

    # input files:
    # ../ad_norm_cohort_stats/<basename>.<type>_desc_stats
    # ouput files:
    # ../ad_norm_cohort_stats/<basename>.cohort_desc_stats

    try:
        subprocess.run([
            f'{python}', f'{scripts}/cohort_desc_stats.py',
            '--ad_norm_cohort_stats_dir', f'{ad_norm_cohort_stats_dir}',
            '--output_dir', f'{ad_norm_cohort_stats_dir}',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"\nCommand failed with error: {e.stderr}")

    
    print('''\nCohort MAG Stats Section''')

    # Note: This section does not inlcude CDS or non-CDS 

    print('''\nCRISPR - Allelic Depth - MAG Desc Stats''')

    if threads > 2: 
        print('\nrunning multiprocessing tool for crispr_mag_desc_stats.py')
        def crispr_mag_desc_stats_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/crispr_mag_desc_stats.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.crispr_ad_norm \
                --output_dir {ad_norm_mag_stats_dir} \
                --assembly_name {mag_name}'''
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print(output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(crispr_mag_desc_stats_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()    
        print('\ncrispr_mag_desc_stats.py - done')
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\nmag_name')
                cmd = f'''{python} {scripts}/crispr_mag_desc_stats.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.crispr_ad_norm \
                    --output_dir {ad_norm_mag_stats_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 
        print('\ncrispr_mag_desc_stats.py - done')


    print('''\nncRNA - Allelic Depth - MAG Desc Stats''')

    if threads > 2: 
        print('\nrunning multiprocessing tool for ncrna_mag_desc_stats.py')
        def ncrna_mag_desc_stats_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/ncrna_mag_desc_stats.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.ncrna_ad_norm \
                --output_dir {ad_norm_mag_stats_dir} \
                --assembly_name {mag_name}'''
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(ncrna_mag_desc_stats_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()    
        print('\nncrna_mag_desc_stats.py - done')
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                cmd = f'''{python} {scripts}/ncrna_mag_desc_stats.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.ncrna_ad_norm \
                    --output_dir {ad_norm_mag_stats_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 
        print('\nncrna_mag_desc_stats.py - done')


    print('''\nDGR - Allelic Depth - MAG Desc Stats''')

    if threads > 2: 
        print('\nrunning multiprocessing tool for dgr_mag_desc_stats.py')
        def dgr_mag_desc_stats_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/dgr_mag_desc_stats.py \
                --ad_norm_file {ad_norm_dir}/{mag_name}.dgr_ad_norm \
                --output_dir {ad_norm_mag_stats_dir} \
                --assembly_name {mag_name}'''
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(dgr_mag_desc_stats_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()    
        print('\ndgr_mag_desc_stats.py - done')
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\nmag_name')
                cmd = f'''{python} {scripts}/dgr_mag_desc_stats.py \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.dgr_ad_norm \
                    --output_dir {ad_norm_mag_stats_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 
        print('\ndgr_mag_desc_stats.py - done')


    print('''\nVIR - Allelic Depth - MAG Desc Stats''')

    script = 'vir_mag_desc_stats.py'
    if threads > 2: 
        print(f'\nrunning multiprocessing tool for {script}')
        def vir_mag_desc_stats_func(pmag):
            mag_name = pmag.split('/')[-1].split(mag_extension)[0]
            print(f'\n{mag_name}')
            cmd = f'''{python} {scripts}/{script} \
                --ad_norm_file {ad_norm_dir}/{mag_name}.vir_ad_norm \
                --output_dir {ad_norm_mag_stats_dir} \
                --assembly_name {mag_name}'''
            try:
                output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                print('\n',output.decode())
            except subprocess.CalledProcessError:
                print(f'\nError: {mag_name}')
        # Read the file to get the list of pmags
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
        # Create a multiprocessing pool and map the function to the list of pmags
        pool = multiprocessing.Pool(threads)
        pool.map(vir_mag_desc_stats_func, pmags)
        # Close the pool and wait for the work to finish
        pool.close()
        pool.join()    
        print(f'\n{script} - done')
    else:
        with open(path_to_mags, 'r') as file:
            pmags = [line.strip() for line in file.readlines()]
            for pmag in pmags:
                mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                print(f'\n{mag_name}')
                cmd = f'''{python} {scripts}/{script} \
                    --ad_norm_file {ad_norm_dir}/{mag_name}.vir_ad_norm \
                    --output_dir {ad_norm_mag_stats_dir} \
                    --assembly_name {mag_name}'''
                subprocess.check_output(cmd, shell=True) 
        print(f'\n{script} - done')
    
        
        
    print('''\nplot_cohort_desc_stats_boxplot.py''')

    # inputs: 
    # .../ad_norm_cohort_stats/<basename>.cohort_desc_stats
    # outputs:
    # .../ad_norm_cohort_stats/plots/<basename>.desc_stats_boxplot.png/pdf

    try:
        subprocess.run([
            f'{python}', f'{scripts}/plot_cohort_desc_stats_boxplot.py',
            '--cohort_desc_stats_file', f'{ad_norm_cohort_stats_dir}/{basename}.cohort_desc_stats',
            '--output_dir', f'{wkdir}/plots/ad_norm_cohort_stats_plots',
            '--basename', f'{basename}'
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"\nCommand failed with error: {e.stderr}")
        

# %% 

if (command == 'sliding_window_plot') or (command == 'all'):

    print(
    '\n#####################################'
    '\n Sliding Window Plots'
    '\n#####################################'
    )

    print('''\nSliding Window Section''')

    win_step_list = []
    if (window != None) and (step != None):
        win_step_list=[[window, step]]
    elif (window != None) and (step == None):
        win_step_list = [[window, window]]
    elif (window == None):
        win_step_list = [[50000, 40000],[5000, 4000]]
    print(win_step_list)

    '''sliding_window.py'''

    script = 'sliding_window.py'

    def sliding_window_func(pmag):
        mag_name = pmag.split('/')[-1].split(mag_extension)[0]
        print(f'\n{mag_name}')
        try:
            subprocess.run([
                f'{python}', f'{scripts}/{script}',
                '--ad_norm_file', f'{ad_norm_dir}/{mag_name}.ad_norm',
                '--window', f'{window}',
                '--step', f'{step}',
                '--output_dir', f'{sliding_window_dir}',
                '--assembly_name', f'{mag_name}'
            ], check=True, capture_output=True, text=True)
            print(f'\n{script} - done')
        except subprocess.CalledProcessError as e:
            print(f'\nError: {script}')
            print(f"\nCommand failed with error: {e.stderr}")
        
    for win_step in win_step_list:
        window = win_step[0]
        print(window)
        step = win_step[1]
        print(step)
        if threads > 2: 
            print('\nrunning multiprocessing tool for sliding_window.py')
            # Read the file to get the list of pmags
            with open(path_to_mags, 'r') as file:
                pmags = [line.strip() for line in file.readlines()]
            # Create a multiprocessing pool and map the function to the list of pmags
            pool = multiprocessing.Pool(threads)
            pool.map(sliding_window_func, pmags)
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()    
            print('\nsliding_window.py - done')
        else:
            with open(path_to_mags, 'r') as file:
                pmags = [line.strip() for line in file.readlines()]
                for pmag in pmags:
                    mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                    print(mag_name)
                    sliding_window_func(pmag)
            print('\nsliding_window.py - done')


    print('''\nPlot Scripts''')

    print('''\nplot_mag_hotspot.py''')

    # inputs: 
    # ../sliding_window/<assembly>.w50000_s40000
    # ../ad_norm/<assembly>.ncrna_ad_norm,"".crispr_ad_norm,"".dgr_ad_norm

    # outputs:
    # .../plots//sliding_window/{assembly_name}_w{window}_s{step}_{metric}.png/pdf


    script = 'plot_mag_hotspot.py'

    def plot_mag_hotspot_func(pmag):
        mag_name = pmag.split('/')[-1].split(mag_extension)[0]
        print(f'\n{mag_name}')
        try:
            subprocess.run([
                f'{python}', f'{scripts}/{script}',
                '--sliding_window_dir', f'{sliding_window_dir}',
                '--metric', f'mean_ad_norm',
                '--window', f'{window}',
                '--step', f'{step}',
                '--assembly_name', f'{mag_name}',
                '--ad_norm_dir', f'{ad_norm_dir}',
                '--output_dir', f'{wkdir}/plots/sliding_window_plots'
            ], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"\nCommand failed with error: {e.stderr}")

    for win_step in win_step_list:
        window = win_step[0]
        step = win_step[1]
        if threads > 2: 
            print(f'\nrunning multiprocessing tool for {script}')
            with open(path_to_mags, 'r') as file:
                pmags = [line.strip() for line in file.readlines()]
            # Create a multiprocessing pool and map the function to the list of pmags
            pool = multiprocessing.Pool(threads)
            pool.map(plot_mag_hotspot_func, pmags)
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()    
            print(f'\n{script} - done')
        else:
            with open(path_to_mags, 'r') as file:
                pmags = [line.strip() for line in file.readlines()]
                for pmag in pmags:
                    mag_name = pmag.split('/')[-1].split(mag_extension)[0]
                    print(f'\n{mag_name}')
                    plot_mag_hotspot_func(pmag)
                    

# %% 

# '''plot_codon_dist.py'''

# # NOTE: error: image size of pixels too large... non-issue when script is run directly?

# # input files:
# # ../pos_freq/<mag_name>.cds_ref_pos_freq

# # output files:
# # ../results/cohort_plots/<basename>.codon_dist_plot.png/pdf

# cmd = f'''{python} {scripts}/plot_codon_dist.py \
# --ad_norm_dir {ad_norm_dir} \
# --output_dir {wkdir}/plots/ad_norm_cohort_stats_plots \
# --basename {basename}'''
# subprocess.check_output(cmd, shell=True) 

# %% 
log_file.close()