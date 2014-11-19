#!/usr/bin/env python
###############################################################################
#
# job.py
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
Represents a single job to be run 
"""

import os
import sys
import subprocess
import time
import logging

from exceptions import MissingInputFileError
from log import get_logger

def format_time(t = time.localtime()):
    if isinstance(t, float):
        t = time.localtime(t)

    return time.strftime("%Y/%m/%d %H:%M:%S", t)

class Job:
    # Job type constants
    COMMAND_JOB = 0
    FUNCTION_JOB = 1

    def __init__(self, command = None, func = None, args = (), kwargs = {}, file_inputs = [], file_outputs = [], name = "", force_out_of_date = False, logger = None): 
        self.name = name
        if command is not None:
            self.type = Job.COMMAND_JOB
            self.command = command
        else:
            self.type = Job.FUNCTION_JOB
            self.func = func
            self.fargs = args
            self.fkwargs = kwargs

        self.force_out_of_date = force_out_of_date
        self.file_inputs = file_inputs
        self.file_outputs = file_outputs
        if not logger:
            self.logger = self.get_logger()
        else:
            self.logger = logger
    
    @staticmethod
    def get_logger():
        return get_logger()

    def execute(self, echo_output = True, command_function = None, cf_args = (), cf_kwargs = {}):
        # print >> sys.stderr, "Logging"
        # print >> sys.stderr, self.logger.getEffectiveLevel()
        self.logger.debug("Job '%s' started at %s" % (self.name, format_time()))
        # print >> sys.stderr, "Done"
        self.logger.debug("Job '%s' inputs: %s" % (self.name, repr(self.file_inputs)))
        for input in self.file_inputs:
            self.logger.debug("Job '%s' checking for input file '%s'" % (self.name, input))
            if not os.path.exists(input) or not os.path.isfile(input):
                raise MissingInputFileError("The required input file '" + input + "' for job " + str(self.name) + " does not exist")

        if self.type == Job.COMMAND_JOB:
            self.logger.debug("Job '%s' command: %s'" % (self.name, self.command))
            if command_function is None:
                if echo_output:
                    # print >> sys.stderr, "Opening process"
                    # self.logger.debug("Job '%s' opening process" % str(job))
                    job_process = subprocess.Popen(self.command, shell=True)
                    # print >> sys.stderr, "Waiting for job to finish"
                    # self.logger.debug("Job '%s' waiting for process to finish..." % str(job))
                    job_process.wait()
                    # print >> sys.stderr, "Job is finished"
                    # self.logger.debug("Job '%s' process finished." % str(job))
                else:
                    # Swallow output via IPC pipes
                    job_process = subprocess.Popen(self.command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    output, error = job_process.communicate()

                return job_process.returncode
            else:
                cf_kwargs["command"] = self.command
                result = command_function(*cf_args, **cf_kwargs)
        else:
            result = self.func(*self.fargs, **self.fkwargs)

        self.logger.debug("Job '%s' finished at %s" % (self.name, format_time()))
        for output in self.file_outputs:
            self.logger.debug("Job '%s' output '%s' modification time: %s" % (self.name, output, format_time(os.path.getmtime(output))))
        return result

    def check_inputs(self):
        for input in self.file_inputs:
            if not os.path.exists(input):
                return False

        return True

    def get_missing_inputs(self):
        missing_inputs = []
        for input in self.file_inputs:
            if not os.path.exists(input):
                missing_inputs.append(input)

        return missing_inputs

    def out_of_date(self, comp_in_out = True):
        if self.force_out_of_date:
            self.logger.debug(str(self) + ": Forced out-of-date")
            return True

        if len(self.file_outputs) == 0:
            self.logger.debug(str(self) + ": No outputs, automatically out of date")
            return True

        earliest_output = -1
        for output in self.file_outputs:
            if not os.path.exists(output) and not os.path.isdir(output):
                self.logger.debug(str(self) + ": Output '%s' does not exist - out of date" % (output,))
                return True

            self.logger.debug("Checking output file '%s' for job '%s' : modified %s" % (output, self.name, format_time(os.path.getmtime(output))))
            if earliest_output == -1:
                earliest_output = os.path.getmtime(output)
            else:
                earliest_output = min(earliest_output, os.path.getmtime(output))

        if not len(self.file_inputs):
            self.logger.debug(str(self) + ": All output files present, no input files, not out-of-date")
            return False

        latest_input = 0
        for input in self.file_inputs:
            if os.path.exists(input):
                self.logger.debug("Checking input file '%s' for job '%s' : modified %s" % (input, self.name, format_time(os.path.getmtime(input))))
                latest_input = max(latest_input, os.path.getmtime(input))
            else:
                self.logger.debug(str(self) + ": input '%s' does not exist" % (input,))
                return True

        if not comp_in_out:
            self.logger.debug(str(self) + ": All output files present, file date comparison skipped, not out-of-date")
            return False

        comparison_result = earliest_output < latest_input
        self.logger.debug(str(self) + ": job is %sout of date, file date comparison: (%s < %s)" % (
            ("" if comparison_result else "not "),
            format_time(earliest_output), 
            format_time(latest_input)))

        return comparison_result

    def __str__(self):
        return self.name
