#!/bin/bash

python3 $(pwd)/../../CheckD/scripts/checkd.py \
	all \
	-basename illumina_example \
	-output_dir $(pwd) \
	-mag_dir $(pwd)/illumina_mags \
	-illumina_forward_reads $(pwd)/illumina_reads/illumina_reads_1.fastq \
	-illumina_reverse_reads $(pwd)/illumina_reads/illumina_reads_2.fastq \
	-threads 8

