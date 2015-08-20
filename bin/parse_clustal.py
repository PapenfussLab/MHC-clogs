#!/usr/bin/env python

from mungolite.align import *

aln = Alignment.load("output/all.aln", format="clustal")
print len(aln), aln.numberOfSeqs()
