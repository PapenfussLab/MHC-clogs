#!/usr/bin/env python

"""
rename_fasta.py
"""

from mungolite.fasta import FastaFile
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument("input_filename", type=str, help="Input filename")
parser.add_argument("output_filename", type=str, help="Output filename")
parser.add_argument("mapping_filename", type=str, help="Mapping filename")
args = parser.parse_args()


data = {}
for line in open(args.mapping_filename):
    line = line.strip()
    if line=="" or line[0]=="#": continue
    name, new_name = line.split("\t")
    data[name] = new_name


output_file = FastaFile(args.output_filename, "w")
for h,s in FastaFile(args.input_filename, "r"):
    header_tokens = h.split()
    name = header_tokens[0]
    try:
        new_name = data[name]
        print "fasta_raname.py:", name, new_name
        h = "%s %s" % (new_name, " ".join(header_tokens[0:]))
    except KeyError:
        pass
    output_file.write(h, s)

output_file.close()
