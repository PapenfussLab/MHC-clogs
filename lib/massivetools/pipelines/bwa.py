"""
pipelines.bwa
"""

import os
from massivetools.pipeline import Pipeline
from massivetools.toolkit import Toolkit
from massivetools.utilities import filename_mangler, which
from massivetools.toolkits.bwa import BWA
from massivetools.toolkits.samtools import Samtools
from massivetools.toolkits.picard import Picard


bwa = Bwa(num_threads=12)
samtools = Samtools()
picard = Picard()


def pe_align(genome_prefix, genome_fasta, fastq_filenames, 
            sai_filenames=None, sam_filename=None, bam_filename=None, 
            sorted_bam_filename=None, mark_dup_filename=None,
            output_directory="align_bwa", force_out_of_date=False):
    
    if not sai_filenames:
        sai_filenames = [
            filename_mangler(fastq_filenames[0], directory=output_directory, replace=[("fastq.gz", "sam"), ("_R1", "")]),
            filename_mangler(fastq_filenames[1], directory=output_directory, replace=[("fastq.gz", "sam"), ("_R2", "")])
        ]
    
    if not sam_filename:
        sam_filename = filename_mangler(fastq_filenames[0], directory=output_directory, replace=[("fastq.gz", "sam"), ("_R1", "")])
    
    if not bam_filename:
        bam_filename = filename_mangler(sam_filename, extension="bam")
    
    if not sorted_bam_filename:
        sorted_bam_filename = filename_mangler(bam_filename, suffix="_sorted")
    
    if not mark_dup_filename:
        mark_dup_filename = filename_mangler(sorted_bam_filename, suffix="_remove_dup")
    
    print sam_filename, bam_filename, sorted_bam_filename, mark_dup_filename
    
    jid1 = bwa.pe_align(genome_prefix, fastq_filenames[0], fastq_filenames[1], 
        sai_filenames[0], sai_filenames[1], sam_filename, 
        force_out_of_date=False)
    
    jid2 = samtools.make_bam(input_genome=genome_fasta, 
        input_sam=sam_filename, output_bam=bam_filename, 
        dependencies=[jid1], force_out_of_date=force_out_of_date)
    
    jid3 = samtools.sort(input_filename=bam_filename, 
        out_prefix=os.path.splitext(sorted_bam_filename)[0], 
        dependencies=[jid2], force_out_of_date=force_out_of_date)
    
    jid4 = samtools.index(input_filename=sorted_bam_filename, 
        dependencies=[jid3], force_out_of_date=force_out_of_date)
    
    jid5 = picard.MarkDuplicates(
        input_file=sorted_bam_filename, 
        output_file=mark_dup_filename,
        metrics_file="picard_mark_dup.txt",
        dependencies=[jid4], force_out_of_date=force_out_of_date)
    
    jid6 = samtools.index(input_filename=mark_dup_filename, 
        dependencies=[jid5], force_out_of_date=force_out_of_date)


if __name__=="__main__":
    target_genome_fasta = "/home/users/lab0605/papenfuss/Papenfuss_lab/projects/reference_genomes/human/hg19.fa"
    target_genome_index = "/home/users/lab0605/papenfuss/Papenfuss_lab/projects/reference_genomes/human/hg19.bowtie2"
    fastq_filenames = [
        ["GOT3NC_NoIndex_L008_R1_001.fastq.gz", "GOT3NC_NoIndex_L008_R2_001.fastq.gz"],
    ]
    for _ in fastq_filenames:
        pe_align(_, target_genome_fasta, target_genome_index)
    Pipeline.run()
