#!/usr/bin/env python
###############################################################################
#
# exceptions.py
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

class MassiveToolsError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class MissingInputFileError(MassiveToolsError):
    pass

class DependencyNotFoundError(MassiveToolsError):
    pass

class MissingExecutableError(MassiveToolsError):
    pass

class UnsynchronizedUnavailableError(MassiveToolsError):
    pass

class CyclicDependencyError(MassiveToolsError):
    pass

class NoHeadNodesError(MassiveToolsError):
    pass
    
class NonZeroExitCodeError(MassiveToolsError):
    pass

class NoJobGraphError(MassiveToolsError):
    pass

class MissingArgumentError(MassiveToolsError):
    pass

