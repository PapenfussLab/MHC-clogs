"""
ensembl

ensembl_gene_id  ensembl_transcript_id  associated_gene_name  description  chrom  gene_start  gene_end  strand
"""

import sys
from useful import *


def ensembl_protein_source(description_str):
    desc_fields = ensembl_protein_description(description_str)
    tests = ["chromosome", "scaffold", "genescaffold", "supercontig", "ultracontig"]
    for test in tests:
        if test in desc_fields.keys():
            return ensembl_field_value(desc_fields[test])
    # if no obvious identifier found, just use whole description
    return description_str


def ensembl_field_value(field_str):
    return field_str.split(":")[1]


def get_ensembl_description(description):
    result = {}
    tokens = description.split()
    result["gene"] = tokens[0]
    for attributes in tokens[1:]:
        attribute_tokens = attributes.split(':')
        if len(attribute_tokens) > 1:
            result[attribute_tokens[0]] = ':'.join(attribute_tokens[1:])
    return result


def read_ensembl_sequences(filename):
    sequences = {}
    fasta_file = open(filename, 'r')
    for line in fasta_file:
        line = line.strip()
        if len(line) == 0: continue
        if line.strip()[0] == '>':
            ensembl_attributes = ensembl_protein_description(line[1:])
            sequences[ensembl_attributes["protein_id"]] = ensembl_attributes
    return sequences
