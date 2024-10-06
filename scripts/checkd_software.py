

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
import statistics
import matplotlib
import seaborn
import scipy
import statsmodels.api
import csv

# %% 
def install_and_import(package):
    try:
        __import__(package)
        print(f"{package} is already installed.")
    except ImportError:
        print(f"{package} not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"{package} installed successfully.")
        
# Example usage
required_libraries = ["numpy", "requests", "pandas"]

for library in required_libraries:
    install_and_import(library)



# %% 

'''Pre-sets'''
threads = 8
python = 'python3'
scripts = os.path.dirname(os.path.abspath(__file__))
checkd_bin = '%s/bin' % (scripts.strip('scripts'))
# checkd_bin = checkd_bin.split('/nfs2')[1] # not sure if this is necesary 
'''Software'''
minimap2 = f'{checkd_bin}/minimap2/minimap2'
samtools = f'{checkd_bin}/samtools-1.14/samtools'
bcftools = f'{checkd_bin}/bcftools/bcftools'
fgs = f'{checkd_bin}/FragGeneScan1.31/run_FragGeneScan.pl'
cmscan = f'{checkd_bin}/infernal-1.1.2/bin/cmscan'
esl_seqstat = f'{checkd_bin}/infernal-1.1.2/bin/esl-seqstat'
Rfam_cm = f'{checkd_bin}/Rfam.cm'
Rfam_clanin = f'{checkd_bin}/Rfam.clanin'
mydgr = f'{checkd_bin}/mgDGR'
summarize = f'{checkd_bin}/CRISPRpk/metaCRT/summarize-crispr-new.py'
crisprone = f'{checkd_bin}/CRISPRone/crisprone-local-2022.py'
crisprgraph = f'{checkd_bin}/CRISPRpk/pipelines/crispr_ann_ort.py'


# %%
