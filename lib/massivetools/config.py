#!/usr/bin/env python
###############################################################################
#
# config.py
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
import sys
import traceback

config = None

DEFAULT_CONFIG_DIRECTORY = os.path.join(os.environ["HOME"], ".massivetools")
DEFAULT_CONFIG_FILENAME = "config.py"

def set_defaults(settings):
    settings["max_cluster_queue_size"] = 31
    settings["logging"] = None
    settings["qsub_file_tries"] = 3
    settings["qsub_file_try_delay"] = 1.0
    settings["qsub_use_sync"] = True
    settings["qsub_poll_interval"] = 2.0

def load_configuration(filename=None):
    if filename is None:
        filename = os.path.join(DEFAULT_CONFIG_DIRECTORY, DEFAULT_CONFIG_FILENAME)

    settings = {}
    set_defaults(settings)
    if os.path.isfile(filename):
        try:
            config_file = open(filename, "r")
            exec config_file in settings
            config_file.close()
        except Exception, e:
            print >> sys.stderr, "Error loading configuration file: " + str(e)
            (type, value, tb) = sys.exc_info()
            traceback.print_tb(tb)

            settings = {}
            set_defaults(settings)
            return settings

    return settings


def get_configuration(filename=None):
    global config

    if filename is not None:
        return load_configuration(filename)

    if config is None:
        config = load_configuration()

    return config
