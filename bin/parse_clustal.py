#!/usr/bin/env python

from mungo.align import *

aln = Alignment.load("output/all.aln", format="clustal")
print len(aln), aln.numberOfSeqs()
