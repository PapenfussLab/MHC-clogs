#!/usr/bin/env python
###############################################################################
#
# samtools.py
#
# Copyright (C) 2011 Christopher Davoren
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

import os
from os.path import basename

from ..pipeline import Pipeline
from ..utilities import which
from ..utilities import replace_extension
from ..exceptions import MissingExecutableError
from ..toolkit import Toolkit


class Samtools(Toolkit):
    def __init__(self, pipeline=Pipeline.get_default_instance(), samtools_path=None):
        Toolkit.__init__(self)
        
        if not samtools_path:
            self.samtools_path = which("samtools")
            if not self.samtools_path:
                raise MissingExecutableError("Unable to find 'samtools' in path")
        else:
            self.samtools_path = samtools_path
        
        self.register_tool("make_bam", self.samtools_path + " view -bt %(input_genome)s -o %(output_bam)s %(input_sam)s")
        self.register_tool("index", self.samtools_path + " index %(input_filename)s")
        self.register_tool("sort", self.samtools_path + " sort %(options)s %(input_filename)s %(out_prefix)s", defaults={"options":""})
        self.register_tool("rmdup", self.samtools_path + " rmdup %(options)s %(input_filename)s %(output_filename)s", defaults={"options":""})
        self.register_tool("merge", self.samtools_path + " merge %(options)s %(output_filename)s %(input_filenames)s", defaults={"options":""})
        self.register_tool("filter", "./scripts/filterBam.py %(options)s %(input_filename)s %(output_filename)s", defaults={"options":""})


def SamtoolsToolkit(*args, **kw):
    import sys
    print sys.stderr, "Use Samtools"
    return Samtools(*args, **kw)
