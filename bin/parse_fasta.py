#!/usr/bin/env python

from mungolite.fasta import FastaReader


input_file = FastaReader("output/all.fa")

seq = {}
for h,s in input_file:
    tokens = h.split()
    acc = tokens[0]
    try:
        seq[acc] += 1
    except KeyError:
        seq[acc] = 1

print seq

for acc in seq:
    print acc, seq[acc]
    if seq[acc]>1:
        print "***"

