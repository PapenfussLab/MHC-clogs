import os
from massivetools.pipeline import Pipeline
from massivetools.toolkit import Toolkit
from massivetools.utilities import filename_mangler, which


class Bwa(Toolkit):
    def __init__(self, pipeline=Pipeline.get_default_instance(), aligner_path=None,
                    quality_filter=3, num_threads=12):
        Toolkit.__init__(self)
        
        if not aligner_path:
            self.aligner_path = which("bwa")
            if not self.aligner_path:
                raise MissingExecutableError("Unable to find bwa in path")
        else:
            self.aligner_path = aligner_path
        
        self.register_tool("aln", self.aligner_path \
            + " aln %(options)s -f %(output_sai)s %(genome_prefix)s %(input_fastq)s", 
                defaults={"options": "-q %i -t %i" % (quality_filter, num_threads)})
        
        self.register_tool("samse", self.aligner_path \
            + " samse %(options)s -f %(output_sam)s %(genome_prefix)s %(input_sai)s %(input_fastq)s", 
                defaults={"options": ""})
        
        self.register_tool("sampe", self.aligner_path \
            + " sampe %(options)s -f %(output_sam)s %(genome_prefix)s %(input_sai1)s %(input_sai2)s %(input_fastq1)s %(input_fastq2)s", 
                defaults={"options": ""})
    
    @staticmethod
    def pe_align(genome_prefix, fastq_filename1, fastq_filename2, sai_filename1, sai_filename2, sam_filename,
                pipeline=Pipeline.get_default_instance(), aligner_path=None, 
                quality_filter=3, num_threads=12, force_out_of_date=False):
        
        bwa = BWA(pipeline=pipeline, aligner_path=aligner_path, quality_filter=quality_filter, num_threads=num_threads)
        
        jid1 = bwa.aln(
            genome_prefix=genome_prefix, 
            input_fastq=fastq_filename1, 
            output_sai=sai_filename1, 
            dependencies=[], 
            force_out_of_date=force_out_of_date)
        
        jid2 = bwa.aln(
            genome_prefix=genome_prefix, 
            input_fastq=fastq_filename2, 
            output_sai=sai_filename2, 
            dependencies=[], 
            force_out_of_date=force_out_of_date)
        
        jid3 = bwa.sampe(
            genome_prefix=genome_prefix, 
            input_sai1=sai_filename1, 
            input_sai2=sai_filename2, 
            input_fastq1=fastq_filename1, 
            input_fastq2=fastq_filename2, 
            output_sam=sam_filename, 
            dependencies=[jid1, jid2],
            force_out_of_date=force_out_of_date)
        
        return jid3


def BwaToolkit(*args, **kw):
    import sys
    print sys.stderr, "Use Bwa"
    return Bwa(*args, **kw)
