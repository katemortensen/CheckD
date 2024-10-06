#!/usr/bin/env python3

# Author: Kate Mortensen
# Last Modified: 9/17/2024
# Purpose: Create mini example MAGs and reads

# %% 
from Bio import SeqIO
import pysam
import subprocess
import os
from pbcore.io import BamReader, FastqWriter

# %% 



# def extract_nucleotide_range(fasta_file, contig_id, start, end, output_file):
#     # Parse the FASTA file
#     with open(fasta_file, "r") as handle:
#         for record in SeqIO.parse(handle, "fasta"):
#             if record.id == contig_id:
#                 # Extract the range of nucleotides (start is inclusive, end is exclusive)
#                 sequence_range = record.seq[start-1:end]
#                 # Write the result to an output file
#                 with open(output_file, "w") as output_handle:
#                     output_handle.write(f">{contig_id}_{start}_{end}\n{sequence_range}\n")
#                 print(f"Extracted nucleotides from {start} to {end} for contig {contig_id}:")
#                 print(sequence_range)
#                 return sequence_range

#     print(f"Contig {contig_id} not found in the FASTA file.")
#     return None


# %% 

# # DGR & CRISPR
# mag_name = f'bin.62'
# contig_id = f's0.ctg000338l'
# start_position = 1000
# end_position = 1500

# mag_name = f'bin.4'
# contig_id = f's441.ctg000682c'
# start_position = 250
# end_position = 750

# mag_name = f'bin.4'
# contig_id = f's441.ctg000682c'
# start_position = 1750
# end_position = 2250

# wkdir = f'/home/kmorten/CheckD/'
# fasta_file = f'{wkdir}/hifi_example/mags/{mag_name}.fa'
# fasta_file = f'/home/kmorten/hifi/NM2022/humanO1/MAGs//bin.4.fa'
# output_dir = f'{wkdir}/example/hifi_mags/'
# output_file = f'{output_dir}/example_hifi_{mag_name}.fa'



# %% 
# Illumina MAGs

# # ncRNA & CDS
# mag_name = f'bin.4'
# mag_out_name = f'bin.1'
# contig_id = f'contig.1'
# start_position = 142400
# end_position = 142800

# # virus
# mag_name = f'bin.4'
# mag_out_name = f'bin.2'
# contig_id = f'contig.79'
# start_position = 1
# end_position = 1544

# # CRISPR
# mag_name = f'bin.4'
# mag_out_name = f'bin.3'
# contig_id = f'contig.89'
# start_position = 31
# end_position = 1296


# wkdir = f'/home/kmorten/CheckD/illumina_example/illumina_example_checkd/'
# fasta_file = f'{wkdir}/renamed_mags/{mag_name}.fa'
# output_dir = f'/home/kmorten/CheckD/example/illumina_mags'
# output_file = f'{output_dir}/{mag_out_name}.fa'


# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# extract_nucleotide_range(fasta_file, contig_id, start_position, end_position, output_file)
                
                
# %% 

# Illumina Reads

# wkdir = f'/home/kmorten/CheckD/example'
# bam_path = f'{wkdir}/illumina_example_checkd/reads2mags_minimap'
# mags = ['bin.1','bin.2', 'bin.3']

# fastq1_in = f'/home/kmorten/CheckD/illumina_example/reads/READS_1.fastq'
# fastq2_in = f'/home/kmorten/CheckD/illumina_example/reads/READS_2.fastq'
# fastq1_out = f'{wkdir}/illumina_reads/illumina_reads_1.fastq'
# fastq2_out = f'{wkdir}/illumina_reads/illumina_reads_2.fastq'

# bam_list = []
# for mag in mags:                
#     cmd = f'samtools view -b -q 30 -F 4 -F 256 -f 2 {bam_path}/{mag}_sorted.bam > {wkdir}/{mag}_well_mapped_illumina.bam'
#     bam_list.append(f'{wkdir}/{mag}_well_mapped_illumina.bam')    
#     try:
#         subprocess.check_output(cmd, shell=True)
#     except subprocess.CalledProcessError as e:
#         print(e.output)

# bams = " ".join(bam_list)
# cmd = f'samtools merge {wkdir}/merged_filtered.bam {bams}'
# try:
#     subprocess.check_output(cmd, shell=True)
# except subprocess.CalledProcessError as e:
#     print(e.output)

# cmd = f'samtools fastq -1 {fastq1_out} -2 {fastq2_out} {wkdir}/merged_filtered.bam'
# try:
#     subprocess.check_output(cmd, shell=True)
# except subprocess.CalledProcessError as e:
#     print(e.output)

# def deduplicate_fastq(input_fastq, output_fastq):
#     sequences = set()  # To store unique sequences
#     with open(output_fastq, "w") as output_handle:
#         for record in SeqIO.parse(input_fastq, "fastq"):
#             seq_str = str(record.seq)
#             if seq_str not in sequences:
#                 sequences.add(seq_str)
#                 SeqIO.write(record, output_handle, "fastq")

# # Deduplicate forward and reverse FASTQ files
# deduplicate_fastq(fastq1_out, f'{fastq1_out}.dedup')
# deduplicate_fastq(fastq2_out, f'{fastq2_out}.dedup')   
                
                
# %% 
# HiFi Reads

wkdir = f'/home/kmorten/CheckD/example'
bam_path = f'{wkdir}/hifi_example_checkd/reads2mags_minimap'
mags = ['bin.0','bin.1', 'bin.2']
fastq_in = f'{wkdir}/hifi_reads/SRR15275213_subset.fastq'
fastq_out = f'{wkdir}/hifi_reads/hifi_reads.fastq'
                
bam_list = []
for mag in mags:     
    print(mag)           
    cmd = f'samtools view -b -F 4 -F 256 {bam_path}/{mag}_sorted.bam > {wkdir}/{mag}_well_mapped_hifi.bam'
    bam_list.append(f'{wkdir}/{mag}_well_mapped_hifi.bam')    
    try:
        subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(e.output)


bams = " ".join(bam_list)
cmd = f'samtools merge {wkdir}/merged_filtered.bam {bams}'
try:
    subprocess.check_output(cmd, shell=True)
except subprocess.CalledProcessError as e:
    print(e.output)

# %% 
# cmd = f'samtools fastq -o {fastq_out} {wkdir}/merged_filtered.bam'
# try:
#     subprocess.check_output(cmd, shell=True)
# except subprocess.CalledProcessError as e:
#     print(e.output)


# def bam_to_fastq(bam_file, output_fastq):
#     with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_fastq, "w") as fastq_out:
#         for read in bam:
#             # Filter to only include reads that are HiFi (CCS) based on PacBio tag 'ccs'
#             if read.has_tag("zm"):  # Assuming PacBio tag
#                 seq = read.query_sequence
#                 qual = read.qual
#                 fastq_out.write(f"@{read.query_name}\n{seq}\n+\n{qual}\n")

# # Convert HiFi BAM to FASTQ
# bam_to_fastq(f'{wkdir}/merged_filtered.bam', fastq_out)

# %%
def extract_hifi_reads(bam_file, output_fastq):
    print(bam_file)
    print(output_fastq)
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_fastq, "w") as fastq_out:
        for read in bam:
            # print(read.flag)
            # Check if the read is marked as CCS (typically in HiFi BAM files)
            # and that the read is mapped (not unmapped)
            # if read.has_tag("fp") and not read.is_unmapped:
            print((f"@{read.query_name}\n"))
            fastq_out.write(f"@{read.query_name}\n")
            print((f"{read.query_sequence}\n"))
            fastq_out.write(f"{read.query_sequence}\n")
            fastq_out.write("+\n")
            print(f"{read.qual}\n")
            fastq_out.write(f"{read.qual}\n")

# Run the function for the BAM file
extract_hifi_reads(f'{wkdir}/merged_filtered.bam', fastq_out)


# %% 


                           
def deduplicate_fastq(input_fastq, output_fastq):
    sequences = set()  # To store unique sequences
    with open(output_fastq, "w") as output_handle:
        for record in SeqIO.parse(input_fastq, "fastq"):
            seq_str = str(record.seq)
            if seq_str not in sequences:
                sequences.add(seq_str)
                SeqIO.write(record, output_handle, "fastq")

# Deduplicate forward and reverse FASTQ files
deduplicate_fastq(fastq_out, f'{fastq_out}.dedup')
                
                
                
 
# %%

# samtools view bin.41_SRR15275213.sorted.bam | cut -f1 | sort | uniq > /home/kmorten/CheckD/example/hifi_unique_reads_list.txt
# sort input.txt | uniq > output.txt
# grep -A 3 -Ff reads_list.txt input.fastq > output_subset.fastq

