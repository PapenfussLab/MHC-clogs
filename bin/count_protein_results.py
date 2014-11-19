#!/usr/bin/env python

"""
count_protein_results.py <input-dir>
"""

import os
import glob
from argparse import ArgumentParser
from mungo.fasta import FastaFile


parser = ArgumentParser()
parser.add_argument('input_dir', type=str, help='Input directory (protein_search)')
args = parser.parse_args()

input_filenames = glob.glob(os.path.join(args.input_dir, "results/*_classI_*.fa"))

combined = set()
data = []
for i,input_filename in enumerate(input_filenames):
    data.append(set())
    for h,s in FastaFile(input_filename):
        data[i].add(h)
    print input_filename, len(data[i])
    combined = combined.union(data[i])
print len(combined)
print
print combined.symmetric_difference(data[3])
print combined.symmetric_difference(data[4])
