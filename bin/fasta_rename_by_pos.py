#!/usr/bin/env python

"""
fasta_rename_by_pos.py
"""

from argparse import ArgumentParser
from mungolite.fasta import FastaFile
from srt.intervals import GenomeIntersector


def get_position(tokens):
    for token in tokens:
        if "Position" in token:
            for delim in ["=", ":", "(", ")"]:
                token = token.replace(delim, " ")
            tokens2 = token.split()
            chrom = tokens2[1]
            pos = tokens2[2]
            strand = tokens2[3]
            start, end = [int(x) for x in pos.split("-")]
            return chrom, start, end, strand
    raise Exception("Could not find position")



parser = ArgumentParser()
parser.add_argument("data_filename", type=str, help="Data filename")
parser.add_argument("input_filename", type=str, help="Input filename")
parser.add_argument("output_filename", type=str, help="Output filename")
args = parser.parse_args()


intersector = GenomeIntersector()

# >mdUT1 Chain=chain39 Position=1:705246778-705258088(+) GeneID=None ProteinID=None Score=82.8 E-value=6.8e-24 Length=11311 Comment=No overlapping annotations
for h,s in FastaFile(args.data_filename):
    tokens = h.split()
    name = tokens[0]
    chrom, start, end, strand = get_position(tokens)
    intersector.add((chrom, strand), start, end, name)


output_filename = FastaFile(args.output_filename, "w")
for h,s in FastaFile(args.input_filename):
    tokens = h.split()
    name = tokens[0]

    chrom, start, end, strand = get_position(tokens)
    rs = intersector.find((chrom, strand), start, end)
    if rs:
        new_name = rs[0].value
        print "fasta_rename_UTs_by_pos.py:", name, new_name
        h = "%s %s" % (new_name, " ".join(tokens[1:]))
    output_filename.write(h, s)
