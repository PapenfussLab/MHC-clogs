import os
from os.path import basename

from massivetools.pipeline import Pipeline
from massivetools.utilities import which
from massivetools.utilities import replace_extension
from massivetools.exceptions import MissingExecutableError
from massivetools.toolkit import Toolkit


class GATK(Toolkit):
    def __init__(self, pipeline=Pipeline.get_default_instance(), gatk_path=None):
        Toolkit.__init__(self)
        
        if not gatk_path:
            self.gatk_path = which("gatk")
            if not self.gatk_path:
                raise MissingExecutableError("Unable to find 'gatk' in path")
        else:
            self.gatk_path = gatk_path
        
        self.register_tool("count_reads", self.gatk_path + 
            " -T CountReads -R %(input_genome)s -I %(input_bam)s")
        
        self.register_tool("recalibrate", self.gatk_path +
            " -T BaseRecalibrator %(options)s -R %(input_genome)s -I %(input_bam)s %(known_sites)s -o %(recal_report)s")
        
        self.register_tool("print_reads", self.gatk_path + 
            " -T PrintReads %(options)s -R %(input_genome)s -I %(input_bam)s -BQSR %(recal_report)s -o %(output_bam)s")
