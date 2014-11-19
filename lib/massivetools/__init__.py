#!/usr/bin/env python
###############################################################################
#
# __init__.py
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

from job import Job
from graph import Graph
from pipeline import Pipeline
from execution_engine import ExecutionEngine
from execution_engine import ThreadPooledExecutionEngine
from execution_engine import QSubTPExecutionEngine
from execution_engine import ThreadedExecutionEngine
from execution_engine import QSubExecutionEngine
from execution_engine import MultiprocessExecutionEngine
from execution_engine import ExecutionEngineFactory
from toolkit import Toolkit
from exceptions import MissingInputFileError, DependencyNotFoundError, MissingExecutableError, MissingArgumentError
from log import get_logger
from log import log_to_stderr
from utilities import filename_mangler
from utilities import get_paired_read_files
from utilities import replace_extension
from utilities import replace_directory
from utilities import strip_directory
from utilities import strip_dir_and_ext
from utilities import prepare_directory
from utilities import which
from threadpool import ThreadPool
from threadpool import WorkerThread
from threadpool import WorkRequest

from toolkits.bwa import BwaToolkit
from toolkits.bowtie2 import Bowtie2Toolkit
from toolkits.tophat import TophatToolkit
from toolkits.samtools import SamtoolsToolkit
from toolkits.picard import PicardToolkit


__all__ = [
    "Job",
    "Graph",
    "Pipeline",
    "ExecutionEngine", 
    "ThreadedExecutionEngine",
    "QSubExecutionEngine",
    "ThreadPooledExecutionEngine",
    "QSubTPExecutionEngine",
    "MultiprocessExecutionEngine",
    "ExecutionEngineFactory",
    "Toolkit",
    "BwaToolkit",
    "SamtoolsToolkit",
    "MissingInputFileError", 
    "DependencyNotFoundError",
    "MissingExecutableError",
    "MissingArgumentError",
    "log_to_stderr",
    "get_logger",
    "get_paired_read_files",
    "replace_extension",
    "replace_directory",
    "strip_directory",
    "strip_dir_and_ext",
    "prepare_directory",
    "which",
    "ThreadPool",
    "WorkerThread",
    "WorkRequest",
    "BwaToolkit",
    "SamtoolsToolkit",
    "PicardToolkit",
    "TophatToolkit",
]

