#!/usr/bin/env python
###############################################################################
#
# utilities.py
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
import logging
import glob
import re
import string
import copy


def which(program):
    """
    Emulates the behaviour of the \*nix 'which' command - if the given program exists in the path, 
    the full executable path is returned.

    Credit for this function go to an answer on Stack Overflow:
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    
    return None


def filename_mangler(filename, directory="", replace=None, extension="", suffix="", mkdir=True):
    new_filename = copy.copy(filename)
    if directory: new_filename = replace_directory(new_filename, directory)
    if replace:
        for _match,_replace in replace:
            new_filename = new_filename.replace(_match, _replace)
    if extension: new_filename = replace_extension(new_filename, extension)
    if suffix: new_filename = add_suffix(new_filename, suffix)
    
    if directory and mkdir:
        try:
            os.makedirs(directory)
        except:
            pass
    
    return new_filename


# def filename_mangler(filename, directory="", extension="", append="", mkdir=False):
#     _directory,_filename = os.path.split(filename)
#     _name,_extension = os.path.splitext(_filename)
#     if not directory: directory = _directory
#     if not extension: extension = _extension
#     
#     if directory and mkdir:
#         try:
#             os.makedirs(directory)
#         except:
#             pass
#     
#     return os.path.join(directory, _name) + append + "." + extension


def add_suffix(filename, suffix):
    basename, extension = os.path.splitext(filename)
    return basename + suffix + extension


def replace_extension(filepath, new_extension, new_path=None):
    if len(new_extension) > 0:
        new_filename = os.path.splitext(filepath)[0] + "." + new_extension
    else:
        new_filename = os.path.splitext(filepath)[0]
    
    return (new_filename if new_path is None else os.path.join(new_path, os.path.basename(new_filename)))


def strip_directory(filepath):
    return os.path.split(filepath)[1]


def strip_dir_and_ext(filepath):
    return os.path.splitext(strip_directory(filepath))[0]


def replace_directory(filepath, new_directory):
    return os.path.join(new_directory, strip_directory(filepath))


def prepare_directory(directory, filename=""):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return os.path.join(directory, os.path.split(filename)[1])


def common_prefix(filepath1, filepath2, remove_trailing_underscores=True):
    filename1 = os.path.split(filepath1)[1]
    filename2 = os.path.split(filepath2)[1]
    
    prefix = ""
    for i in xrange(0, min(len(filename1), len(filename2))):
        if not filename1[i] == filename2[i]:
            break
        else:
            prefix += filename1[i]
    
    if remove_trailing_underscores:
        while prefix[-1] == "_":
            prefix = prefix[:-1]
    
    return prefix


def extract_fields(astring):
    field_list = []
    
    formatter = string.Formatter()
    new_fields = formatter.parse(astring)
    
    # print astring
    
    for (literal, field_name, format_spec, conversion) in new_fields:
        if field_name:
            # print "Found standard field name: " + field_name
            field_list.append(field_name)
    
    # print field_list
    
    old_field_matcher = re.compile("(^|[^%])%\((.*?)\)[diouxXeEfFgGcrs]")
    pos = 0
    match = old_field_matcher.search(astring, pos)
    while match:
        # print match.group(0)
        pos = match.end()
        # print "Found old field name: " + match.group(2) + " with whole match of " + match.group(0) + " and pos is now " + str(pos) + " (" + astring[pos:] + ")"
        field_list.append(match.group(2))
        match = old_field_matcher.search(astring, pos)
    
    # print "whole field list: " + str(field_list)
    return field_list


def get_paired_read_files(path, prefix = "", suffix = ""):
    if not path:
        return None
    
    # Get a list of all matching files (including optionally specified prefix and suffix)
    # If path does not point to an existing directory, it is assumed to be a glob pattern
    # (prefix and suffix are ignored in this case)
    if os.path.isdir(path):
        file_list = glob.glob(os.path.join(path, prefix + "*" + suffix))
    else:
        file_list = glob.glob(path)
    
    # Must match at least one file
    if not len(file_list):
        return None
    
    # Must match an even number of files
    if len(file_list) % 2:
        return None
    
    # Sort the list - paired files will now be adjacent assuming numerical file naming
    file_list.sort()
    
    # Group pairs into a list of tuples
    pairs = zip(file_list[::2], file_list[1::2])
    
    return [(f1, f2, common_prefix(f1, f2)) for f1, f2 in pairs] 
