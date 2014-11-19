#!/usr/bin/env python
###############################################################################
#
# log.py
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
import logging

_log_set = False

def get_logger():
    return logging.getLogger("massivetools.debug")

def get_output_logger():
    return logging.getLogger("massivetools.output")

def get_multiprocess_logger():
    return logging.getLogger("massivetools.multiprocess")

def log_to_stderr(log_level = logging.DEBUG):
    global _log_set

    if _log_set:
        return

    log_format = "%(asctime)s.%(msecs)03d  %(levelname)-5s : %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"

    logger = get_logger()
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(log_format, date_format))
    logger.addHandler(handler)
    logger.setLevel(log_level)

    mlogger = get_multiprocess_logger()
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(log_format, date_format))
    mlogger.addHandler(handler)
    mlogger.setLevel(log_level)

    _log_set = True

class NullHandler(logging.Handler):
    def emit(self, record):
        pass

get_logger().addHandler(NullHandler())
get_multiprocess_logger().addHandler(NullHandler())

output_logger = get_output_logger()
output_logger_handler = logging.StreamHandler()
output_logger_handler.setFormatter(logging.Formatter("\t%(message)s", "%Y-%m-%d %H:%M:%S"))
# output_logger_handler.setFormatter(logging.Formatter("[%(asctime)s.%(msecs)03d] %(message)s", "%Y-%m-%d %H:%M:%S"))
output_logger.addHandler(output_logger_handler)
output_logger.setLevel(logging.INFO)


