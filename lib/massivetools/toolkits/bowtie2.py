#!/usr/bin/env python

"""
align.py
"""

import os
from massivetools.pipeline import Pipeline
from massivetools.toolkit import Toolkit
from massivetools.utilities import which
from massivetools.exceptions import MissingExecutableError


class Bowtie2(Toolkit):
    def __init__(self, pipeline=Pipeline.get_default_instance(), aligner_path=None,
                    local=True, qc_filter=True, num_threads=12):
        Toolkit.__init__(self)
        
        if not aligner_path:
            self.aligner_path = which("bowtie2")
            if not self.aligner_path:
                raise MissingExecutableError("Unable to find bowtie2 in path")
        else:
            self.aligner_path = aligner_path
        
        self.register_tool("pe_align", self.aligner_path \
            + " %(options)s -x %(input_genome)s -1 %(input_fastq_1)s -2 %(input_fastq_2)s -S %(output_filename)s", 
                defaults={
                    "options": ' '.join( [
                        "--local" if local else "", 
                        "--qc-filter" if qc_filter else "",
                        "-p %i" % num_threads])})


def Bowtie2Toolkit(*args, **kw):
    import sys
    print sys.stderr, "Use Bowtie2"
    return Bowtie2(*args, **kw)
