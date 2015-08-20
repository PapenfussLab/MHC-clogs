#!/usr/bin/env python

"""
count_proteins_per_species.py <input-fasta-file>
"""

import os
import glob
from argparse import ArgumentParser
from mungolite.fasta import FastaFile


parser = ArgumentParser()
parser.add_argument('input_fasta', type=str, help='Input fasta file')
args = parser.parse_args()

data = {}
for h,s in FastaFile(args.input_fasta):
    spp = h.split("_")[0]
    try:
        data[spp] += 1
    except KeyError:
        data[spp] = 1

for spp in data:
    print "%s\t%i" % (spp, data[spp])
