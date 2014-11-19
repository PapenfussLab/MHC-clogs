#!/usr/bin/env python

import os
from os.path import basename

from massivetools.pipeline import Pipeline
from massivetools.utilities import which
from massivetools.utilities import replace_extension
from massivetools.exceptions import MissingExecutableError
from massivetools.toolkit import Toolkit


class Tophat(Toolkit):
    def __init__(self, pipeline=Pipeline.get_default_instance(), tophat_path=None):
        Toolkit.__init__(self)

        if not tophat_path:
            self.tophat_path = which("tophat")
            if not self.tophat_path:
                raise MissingExecutableError("Unable to find 'tophat' in path")
        else:
            self.tophat_path = tophat_path

        self.register_tool("tophat", self.tophat_path + " %(options)s -o %(output_dir)s %(index_base)s %(reads1)s %(reads2)s",
            defaults={"reads2": "", "output_dir":"./tophat_out"})


def TophatToolkit(*args, **kw):
    import sys
    print sys.stderr, "Use Tophat"
    return Tophat(*args, **kw)
