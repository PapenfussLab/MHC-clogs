#!/usr/bin/env python

"""
rename_seqs.py <work_dir>
"""

from argparse import ArgumentParser
from mhc.data import *


parser = ArgumentParser()
parser.add_argument("work_dir", type=str, help="Work directory")
args = parser.parse_args()


search_species = sorted(species_short.keys())
for species in search_species:
    print species
    input_merged_filename = os.path.join(args.work_dir, species, data_species_subdir, "custom_hmmer2_merged_classI_proteins.fa")
    output_merged_filename = os.path.join(args.work_dir, species, data_species_subdir, "custom_hmmer2_merged_classI_proteins_renamed.fa")
    os.system("fasta_unique_labels.py %s %s" % (input_merged_filename, output_merged_filename))


# Rename opossum UT chains by position
# data_filename = "data/UTs_ensembl_R68/all_UTs.fa"
data_filename = "data/md_chains_Ensembl_R68_named.fa"
species = "opossum"
input_merged_filename = os.path.join(args.work_dir, species, data_species_subdir, "custom_hmmer2_merged_classI_proteins_renamed.fa")
output_merged_filename = os.path.join(args.work_dir, species, data_species_subdir, "custom_hmmer2_merged_classI_proteins_renamed2.fa")
os.system("fasta_rename_by_pos.py %s %s %s" % (data_filename, input_merged_filename, output_merged_filename))
os.system("mv %s %s" % (output_merged_filename, input_merged_filename))


# Rename other species and maybe other opossum genes by mapping
mapping_filename = "data/Fasta_accession_mapping.txt"
for species in ["opossum", "wallaby", "devil", "platypus"]:
    print "###", species
    input_merged_filename = os.path.join(args.work_dir, species, data_species_subdir, "custom_hmmer2_merged_classI_proteins_renamed.fa")
    output_merged_filename = os.path.join(args.work_dir, species, data_species_subdir, "custom_hmmer2_merged_classI_proteins_renamed2.fa")
    os.system("fasta_rename.py %s %s %s" % (input_merged_filename, output_merged_filename, mapping_filename))
    os.system("mv %s %s" % (output_merged_filename, input_merged_filename))
