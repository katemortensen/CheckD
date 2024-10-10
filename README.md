# ![Logo](./diagrams/checkD_logo.drawio.png) 


## Overview

CheckD is a metagenomics tool used to assess the diversity of metagenome-assembled genomes (MAGs).

[ADD QUICK DETAILS ABOUT OUTPUT]

[MOVE BELOW TO ANOTHER PAGE THAT BREAKS DOWN THE STEPS]

 CheckD can be used for MAGs assembled from short reads or long reads. CheckD uses software to identify various regions including coding domain sequences (CDS), non-coding RNA (ncRNA), Diversity Generating Retroelements (DGR) systems, CRISPR regions and artifacts, and viral DNA. 

## Software Requirements: 
- Python 3.8 
- apptainer 1.3.3 (Singularity)


## Installation and Setup 

Note that CheckD's bin.tar.gz file (below) is nearly 7G zipped and 15G unzipped.

```
git clone https://github.com/katemortensen/CheckD.git
wget -P ./CheckD/. https://omics.informatics.indiana.edu/CheckD/bin.tar.gz
tar -xzvf ./CheckD/bin.tar.gz -C ./CheckD/.
rm ./CheckD/bin.tar.gz
```
## Conda Environment

```
conda env create -f ./CheckD/scripts/checkd_conda_env.yml 
conda activate checkd_env
```

## Test Run

The following test runs stress test the functionality of CheckD on the user's platform, as well as, provide examples of the appropriate setup.

```
cd CheckD/example/

sh ./illumina_example.sh
ls illumina_example_checkd

sh ./hifi_example.sh
ls hifi_example_checkd

cd ../../ 
```

## Usage

Inputs:
- Metagenome-assembled genomes (MAGs)
- Illumina or HiFi reads from which the MAGs were assembled


Illumina:

```
python3 <path_to>/CheckD/scripts/checkd.py \
	all
        -basename healthy_human_gut \
        -output_dir <path_to_output> \
        -mag_dir <path_to_mags> \
        -illumina_forward_reads <path_to_reads1.fastq> \
        -illumina_reverse_reads <path_to_reads2.fastq> \
        -threads 8
```

PacBio HiFi:

```
python3 <path_to>/CheckD/scripts/checkd.py \
	all
        -basename lake_griffy_soil \
        -output_dir <path_to_output> \
        -mag_dir <path_to_mags> \
        -hifi_reads <path_to_hifi_reads.fastq> \
        -threads 8
```

The user can also select window and step sizes for output plots using the "-step" and "-window" arguments. Otherwise, default sizes are applied (50k window, 40k step and 5k window, 4k step).

For more details, see the link below.

[usage tutorial](https://github.com/katemortensen/CheckD/blob/main/USAGE_TUTORIAL.md)

