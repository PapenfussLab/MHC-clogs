#!/usr/bin/env python
###############################################################################
#
# toolkit.py
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

import sys
import os
import types
import inspect
import glob

from pipeline import Pipeline
from exceptions import MissingArgumentError
from utilities import extract_fields

class Toolkit:
    def __init__(self, name="<default>", pipeline=Pipeline.get_default_instance()):
        self.name = name
        self.registered_tools = {}
        self.pipeline = pipeline

    def run_tool(self, tool_name, job_name=None, dependencies=None, force_out_of_date=False, extra_inputs=[], extra_outputs=[], **kwargs):
        def process_spec(spec, dictionary):
            spec = spec % dictionary
            if spec.find("*") >= 0 or spec.find("?") >= 0 or spec.find("[") >= 0 or spec.find("]") >= 0:
                spec = glob.glob(spec)
            else:
                spec = [spec]

            return spec

        tool = self.registered_tools[tool_name]
        input_files = [kwargs[input_name] for input_name in tool.get_input_names()] + extra_inputs
        output_files = [kwargs[output_name] for output_name in tool.get_output_names()] + extra_outputs

        for input_spec in tool.inputs:
            input_files += process_spec(input_spec, kwargs)

        for output_spec in tool.outputs:
            output_files += process_spec(output_spec, kwargs)

        # print "Inputs for job", tool_name, ":", input_files
        # print "Outputs for job", tool_name, ":", output_files

        if not job_name:
            job_name = ("%s%s%s%s%s%s" %
                (self.name,
                 " " if self.name and tool.name else "",
                 tool.name,
                 (" " + ", ".join(map(os.path.basename, input_files))) if len(input_files) else "",
                 " -> " if len(map(os.path.basename, input_files)) and len(output_files) else "",
                 ", ".join(map(os.path.basename, output_files)) if len(output_files) else ""))

        tool_args = kwargs

        for expected_arg in tool.args:
            if not expected_arg in tool_args.keys():
                if not expected_arg in tool.defaults.keys():
                    raise MissingArgumentError("Tool '%s' from toolkit '%s' expects argument named '%s' which was not provided" % (tool.name, self.name, expected_arg))
                else:
                    tool_args[expected_arg] = tool.defaults[expected_arg]
        
        return self.pipeline.add_job(tool.cmd_template, args=kwargs, name=job_name, input_files=input_files, output_files=output_files, dependencies=dependencies, force_out_of_date=force_out_of_date)

    def register_tool(self, name, cmd_template, args=None, defaults={}, inputs=[], outputs=[]):
        def dispatcher(self, **kwargs):
            return self.run_tool(name, **kwargs)

        # print extract_fields(cmd_template)

        if not args:
            args = extract_fields(cmd_template)

        self.registered_tools[name] = Tool(name, self, cmd_template, args, defaults, inputs, outputs)

        method_obj = types.MethodType(dispatcher, self, self.__class__)
        setattr(self, name, method_obj)

class Tool:
    def __init__(self, name, owner, cmd_template, args={}, defaults={}, inputs=[], outputs=[]):
        self.name = name
        self.owner = owner
        self.cmd_template = cmd_template
        self.args = args
        self.defaults = defaults
        self.inputs = inputs
        self.outputs = outputs

    def get_input_names(self):
        input_names = []
        for arg in self.args:
            if arg.startswith("input"):
                input_names.append(arg)

        return input_names

    def get_output_names(self):
        output_names = []
        for arg in self.args:
            if arg.startswith("output"):
                output_names.append(arg)

        return output_names
