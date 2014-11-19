#!/usr/bin/env python

import os
from os.path import basename

from massivetools.pipeline import Pipeline
from massivetools.utilities import which
from massivetools.utilities import replace_extension
from massivetools.exceptions import MissingExecutableError
from massivetools.toolkit import Toolkit


JAVA_PATH = "/usr/local/bioinf/bin/java"
PICARD_PATH = "/usr/local/bioinfsoftware/picard-tools/current/jars/"

class PicardToolkit(Toolkit):
    def __init__(self, pipeline=Pipeline.get_default_instance(), java_path=JAVA_PATH, picard_path=PICARD_PATH):
        Toolkit.__init__(self)
        
        if not java_path:
            self.java_path = which("java")
            if not self.java_path:
                raise MissingExecutableError("Unable to find 'java' in path")
        else:
            self.java_path = java_path
            self.picard_path = picard_path
        
        self.register_tool("mark_duplicates", self.java_path + " -jar -Xmx4g " + self.picard_path + "MarkDuplicates.jar \
                INPUT=%(input_file)s OUTPUT=%(output_file)s METRICS_FILE=%(metrics_file)s REMOVE_DUPLICATES=%(remove_duplicates)s ASSUME_SORTED=%(assume_sorted)s %(others)s", \
                defaults = {
                    "remove_duplicates": "true", 
                    "assume_sorted": "true",
                    "metrics_file": 
                    "others": "VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900"
                }
        )
        
        self.register_tool("sort_sam", self.java_path +" -jar -Xmx4g " + self.picard_path + "SortSam.jar \
                INPUT=%(input_file)s OUTPUT=%(output_file)s SORT_ORDER=%(sort_order)s %(others)s", \
                defaults = {"sort_order":"coordinate", "others":"VALIDATION_STRINGENCY=SILENT"})


