#!/bin/bash

python3 $(pwd)/../../CheckD/scripts/checkd.py \
	all \
	--basename hifi_example \
	--output_dir $(pwd) \
	--mag_dir $(pwd)/hifi_mags \
	--hifi_reads $(pwd)/hifi_reads/hifi_reads.fastq \
	--threads 8

