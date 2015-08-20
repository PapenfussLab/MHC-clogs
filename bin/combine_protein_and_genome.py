#!/usr/bin/env python

"""
combine_protein_and_genome.py
"""

from mhc.data import *
from mungolite.fasta import FastaFile


def extract_geneid(header):
    tokens = header.split()
    for token in tokens:
        if "GeneID" in token:
            return token.split("=")[1]
    return None


models = "custom_hmmer2"
input_protein_dir = "output_proteins"
input_genome_dir = "output_genomes"
output_dir = "output"

search_species = sorted(species_short.keys())


for species in search_species:
    print species

    if models=="hmmer2":
        input_protein_filename = os.path.join(input_protein_dir, species, data_species_subdir, "hmmer2_classI_proteins.fa")
        input_genome_filename = os.path.join(input_genome_dir, species, data_species_subdir, "hmmer2_chains_pep.fa")
        output_filename = os.path.join(output_dir, "hmmer2_%s.fa") % species
    elif models=="custom_hmmer2":
        input_protein_filename = os.path.join(input_protein_dir, species, data_species_subdir, "custom_hmmer2_classI_proteins.fa")
        input_genome_filename = os.path.join(input_genome_dir, species, data_species_subdir, "custom_hmmer2_chains_pep.fa")
        output_filename = os.path.join(output_dir, "custom_hmmer2_%s.fa") % species
    else:
        print "No such model"
        sys.exit(-1)

    sequences = {}
    for h,s in FastaFile(input_genome_filename):
        geneid = extract_geneid(h)
        sequences[geneid] = (h,s)

    # Overwrite genomic sequences if proteins are available
    for h,s in FastaFile(input_protein_filename):
        geneid = extract_geneid(h)
        sequences[geneid] = (h,s)

    writer = FastaFile(output_filename, "w")
    for geneid in sequences:
        h,s = sequences[geneid]
        writer.write(h,s)
    writer.close()
