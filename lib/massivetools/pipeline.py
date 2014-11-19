#!/usr/bin/env python
###############################################################################
#
# pipeline.py
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
"""
Represents a pipeline of jobs - provides a generic interface for making job graphs with dependencies
"""

# from types import *

from job import Job
from graph import Graph
from exceptions import DependencyNotFoundError
from execution_engine import ExecutionEngine
from execution_engine import ExecutionEngineFactory
from config import get_configuration
from log import log_to_stderr

class Pipeline:

    def __init__(self, verbose = 0):
        self.job_counter = 0
        self.jobs = {}
        self.job_graph = Graph()
        self.previous_job_id = -1
        self.verbose = verbose
        
    def add_job(self, command, args = (), name = None, dependencies = None, input_files = [], output_files = [], depend_last_job = True, force_out_of_date = False):
        if name == None:
            name = "Job " + str(self.job_counter)

        new_job = Job(command = command % args, name = name, file_inputs = input_files, file_outputs = output_files, force_out_of_date = force_out_of_date)
        
        # Store new job both by ID and by name
        self.jobs[self.job_counter] = new_job
        self.jobs[name] = new_job

        self.job_graph.add_node(new_job)

        if dependencies:
            for dependency in dependencies:
                if dependency in self.jobs.keys():
                    job_parent = self.jobs[dependency]
                    self.job_graph.add_arc(job_parent, new_job)
                    if self.verbose > 0: print "*** Dependency registered: " + str(job_parent) + " -> " + str(new_job)
                else:
                    raise DependencyNotFoundError("Specified dependency '" + str(dependency) + "' is not in pipeline job list")

        if self.previous_job_id >= 0 and depend_last_job and dependencies is None:
            self.job_graph.add_arc(self.jobs[self.previous_job_id], new_job)

        self.previous_job_id = self.job_counter
        self.job_counter += 1

        return self.previous_job_id

    def execute(self, engine = ExecutionEngineFactory.get_best_engine()):
        config = get_configuration()
        log_setting = config.get("logging", None)
        if log_setting == "stderr":
            log_to_stderr()

        engine.set_graph(self.job_graph)
        return engine.execute()

    @staticmethod
    def run(engine = ExecutionEngineFactory.get_best_engine(), max_parallel=None, max_parallel_buffer=0, terminate_on_nonzero_exit=True, stop_all_on_failure=True):
        if max_parallel is not None:
            engine.set_max_parallel(max_parallel, max_parallel_buffer)
        elif max_parallel_buffer > 0:
            engine.set_max_parallel_buffer(max_parallel_buffer)

        engine.terminate_on_nonzero_exit = terminate_on_nonzero_exit
        engine.stop_all_on_failure = stop_all_on_failure

        engine.set_graph(Pipeline.get_default_instance().job_graph)
        return Pipeline.get_default_instance().execute(engine)

    @staticmethod
    def get_default_instance():
        if not hasattr(Pipeline, "default_pipeline"):
            Pipeline.default_pipeline = Pipeline()

        return Pipeline.default_pipeline
