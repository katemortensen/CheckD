# CheckD

CheckD is a metagenomics tool used to assess the diversity of metagenome-assembled genomes (MAGs). 

CheckD can be used for MAGs assembled from short reads or long reads. 

# Setup 

Software Requirements: 
python3
apptainer (singularity)

Inputs:
Metagenome-assembled genomes (MAGs)
Illumina or HiFi reads

# Overview

CheckD provides insight into the diversity of a MAG based on allelic depth of mapped reads. 


# Usage

Illumina short reads:

```
python3 <path_to>/CheckD/scripts/checkd.py \
        -basename human_gut \
        -output_dir <path_to_output> \
        -mag_dir <path_to>/human_gut_mags \
        -illumina_forward_reads <path_to>/human_gut_reads_1.fastq \
        -illumina_reverse_reads <path_to>/human_gut_reads_2.fastq \
        -threads 8
```

PacBio HiFi:

```
python3 <path_to>/CheckD/scripts/checkd.py \
        -basename soil \
        -output_dir <path_to_output> \
        -mag_dir <path_to>/soil_mags \
        -hifi_reads <path_to>/hifi_reads.fastq \
        -threads 8
```

The user can also select window and step sizes for output plots using the "-step" and "-window" arguments. Otherwise, default sizes are applied (50k window, 40k step and 5k window, 4k step).



[usage tutorial](https://github.com/katemortensen/Hypervariable-region-aware-co-assembly-of-metagenomes/blob/585052fb24fd959ab47ce077d4daa5e0ba7511d7/hypervar-pipeline/usage_tutorial.md)

## Pipeline Details

![alt text](https://github.com/katemortensen/Hypervariable-region-aware-co-assembly-of-metagenomes/blob/9f18217452eb15362cc46c4c5fc05f9e9709498e/images/repeat_guided_spacer_disc_20220510-Page-1.drawio.png)
