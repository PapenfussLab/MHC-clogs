#!/usr/bin/env python

"""
count_hits_per_species.py
"""

import os
from mungo.fasta import FastaFile
from mhc.data import *


search_species = [
    'human',
    'mouse',
    'dog',
    'cow',
    'opossum',
    'wallaby',
    'devil',
    'platypus',
    'chicken',
    'zebrafinch',
    'turkey',    
    'lizard',
    'x_tropicalis',
    'zebrafish',
    'tetraodon',
    'lamprey',
    'c_intestinalis',
    'fruitfly',
    's_cerevisiae']


count = {}

# protein search
for species in search_species:
    input_filename = os.path.join("output/proteins", species, data_species_subdir, "custom_hmmer2_classI_proteins.fa")
    count[species] = [0]
    for h,s in FastaFile(input_filename):
        count[species][0] += 1

# genome search
for species in search_species:
    input_filename = os.path.join("output/genomes", species, data_species_subdir, "custom_hmmer2_chains.fa")
    count[species].append(0)
    for h,s in FastaFile(input_filename):
        count[species][1] += 1

# merged search
for species in search_species:
    input_filename = os.path.join("output/genomes", species, data_species_subdir, "custom_hmmer2_merged_classI_proteins.fa")
    count[species].append(0)
    for h,s in FastaFile(input_filename):
        count[species][2] += 1

# Output
print "Species\tProtein\tGenome\tMerged"
for species in search_species:
    print "%s\t%i\t%i\t%i" % (species, count[species][0], count[species][1], count[species][2])
