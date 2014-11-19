"""
mhc.tk
"""

import os
from mhc.data import HMMER2_DIR, HMMER3_DIR, NUMBER_HMMER_CPUs
from massivetools.toolkit import Toolkit
from massivetools.pipeline import Pipeline


class HMMerTk(Toolkit):
    def __init__(self, pipeline=Pipeline.get_default_instance(), aligner_path=None):
        Toolkit.__init__(self)

        self.register_tool("hmmsearch2",
            "hmmsearch %(options)s %(input_model)s %(input_fasta)s > %(output_filename)s",
            defaults={"options":"--cpu %i" % NUMBER_HMMER_CPUs})

        self.register_tool("hmmsearch3",
            "hmmsearch %(options)s %(max_option)s -o %(output_filename)s --tblout %(output_table)s --domtblout %(output_domtable)s %(input_model)s %(input_fasta)s",
            defaults={"options":"--cpu %i" % NUMBER_HMMER_CPUs, "max_option": ""})

        self.register_tool("translate",
            "translate -a -q -o %(output_6frame_filename)s %(input_genome_filename)s")

        self.register_tool("split",
            "fastaBlockSplit.py %(input_genome_filename)s %(output_block_genome_filename)s")

        self.register_tool("remove", "rm %(input_filename)s")
