#!/usr/bin/env python

"""
align.py
"""

import os
from massivetools.pipeline import Pipeline
from massivetools.toolkit import Toolkit
from massivetools.utilities import filename_mangler, which
from massivetools.toolkits.bowtie2 import Bowtie2
from massivetools.toolkits.samtools import Samtools
from massivetools.toolkits.picard import Picard


bowtie2 = Bowtie2(num_threads=24)
samtools = Samtools()
picard = Picard()


def pe_align(fastq_filenames, target_genome_fasta, target_genome_index):
    sam_filename = filename_mangler(fastq_filenames[0], directory="align", replace=[("fastq.gz", "sam"), ("_R1", "")])
    bam_filename = filename_mangler(sam_filename, extension="bam")
    sorted_bam_filename = filename_mangler(bam_filename, suffix="_sorted")
    mark_dup_filename = filename_mangler(sorted_bam_filename, suffix="_remove_dup")
    print sam_filename, bam_filename, sorted_bam_filename, mark_dup_filename
    
    jid1 = bowtie2.pe_align(
        input_genome=target_genome_index, 
        input_fastq_1=fastq_filenames[0], 
        input_fastq_2=fastq_filenames[1], 
        output_filename=sam_filename, force_out_of_date=False)
    
    jid2 = samtools.make_bam(input_genome=target_genome_fasta, 
        input_sam=sam_filename, output_bam=bam_filename, 
        dependencies=[jid1], force_out_of_date=False)
    
    jid3 = samtools.sort(input_filename=bam_filename, 
        out_prefix=os.path.splitext(sorted_bam_filename)[0], 
        dependencies=[jid2], force_out_of_date=False)
    
    jid4 = samtools.index(input_filename=sorted_bam_filename, 
        dependencies=[jid3], force_out_of_date=False)
    
    jid5 = picard.MarkDuplicates(
        input_file=sorted_bam_filename, 
        output_file=mark_dup_filename,
        metrics_file="picard_mark_dup.txt",
        dependencies=[jid4], force_out_of_date=False)
    
    jid6 = samtools.index(input_filename=mark_dup_filename, 
        dependencies=[jid5], force_out_of_date=False)


if __name__=="__main__":
    target_genome_fasta = "/home/users/lab0605/papenfuss/Papenfuss_lab/projects/reference_genomes/human/hg19.fa"
    target_genome_index = "/home/users/lab0605/papenfuss/Papenfuss_lab/projects/reference_genomes/human/hg19.bowtie2"
    fastq_filenames = [
        ["GOT3NC_NoIndex_L008_R1_001.fastq.gz", "GOT3NC_NoIndex_L008_R2_001.fastq.gz"],
    ]
    for _ in fastq_filenames:
        pe_align(_, target_genome_fasta, target_genome_index)
    Pipeline.run()
