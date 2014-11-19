#!/usr/bin/env python

"""
rename_seqs.py
"""

from argparse import ArgumentParser
from mungo.fasta import FastaFile


parser = ArgumentParser()
parser.add_argument("input_filename", type=str, help="Input filename")
parser.add_argument("output_filename", type=str, help="Output filename")
args = parser.parse_args()


order = []
seqs = {}
for h, s in FastaFile(args.input_filename):
    h = h.replace("dr_si:", "dr_")
    h = h.replace("dr_zgc:", "dr_")

    tokens = h.split()
    gene_symbol = tokens[0]

    try:
        seqs[gene_symbol].append((h, s))
    except KeyError:
        seqs[gene_symbol] = [(h, s)]
        order.append(gene_symbol)


writer = FastaFile(args.output_filename, 'w')
count = 0
for gene_symbol in order:
    count += 1
    i = 0
    if len(seqs[gene_symbol]) == 1:
        h, s = seqs[gene_symbol][0]
        writer.write(h, s)
    elif len(seqs[gene_symbol]) > 1:
        for h, s in seqs[gene_symbol]:
            i += 1
            tokens = h.split()
            new_gene_symbol = "%s.%i" % (gene_symbol, i)
            # print count, i, gene_symbol, new_gene_symbol, tokens[1]
            tokens = [new_gene_symbol] + tokens[1:]
            h = " ".join(tokens)
            writer.write(h, s)
writer.close()
