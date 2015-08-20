#!/usr/bin/env python

from mungolite.fasta import FastaFile


filenames = [
    "output/species/opossum.fa",
    "output/species/wallaby.fa",
    "output/species/devil.fa",
    "output/species/platypus.fa"
]

save = []
for filename in filenames:
    for h,s in FastaFile(filename):
        if "UT" in h:
            save.append((h, s))

save.sort()

f = FastaFile("output/phylo/UTs.fa", "w")
for h,s in save:
    f.write(h, s)
f.close()
