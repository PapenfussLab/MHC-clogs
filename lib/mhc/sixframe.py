#!/usr/bin/env python

import os
import sys
import optparse

codon_table = {
    "TTT" :	"F",
    "TTC" :	"F",
    "TTA" :	"L",
    "TTG" :	"L",
    "TCT" :	"S",
    "TCC" :	"S",
    "TCA" :	"S",
    "TCG" :	"S",
    "TAT" :	"Y",
    "TAC" :	"Y",
    "TAA" :	"*", # Stop codon (ochre)
    "TAG" :	"*", # Stop codon (amber)
    "TGT" :	"C",
    "TGC" :	"C",
    "TGA" :	"*", # Stop codon (opal/umber)
    "TGG" :	"W",
    "CTT" :	"L",
    "CTC" :	"L",
    "CTA" :	"L",
    "CTG" :	"L",
    "CCT" :	"P",
    "CCC" :	"P",
    "CCA" :	"P",
    "CCG" :	"P",
    "CAT" :	"H",
    "CAC" :	"H",
    "CAA" :	"Q",
    "CAG" :	"Q",
    "CGT" :	"R",
    "CGC" :	"R",
    "CGA" :	"R",
    "CGG" :	"R",
    "ATT" :	"I",
    "ATC" :	"I",
    "ATA" :	"I",
    "ATG" :	"M",
    "ACT" :	"T",
    "ACC" :	"T",
    "ACA" :	"T",
    "ACG" :	"T",
    "AAT" :	"N",
    "AAC" :	"N",
    "AAA" :	"K",
    "AAG" :	"K",
    "AGT" :	"S",
    "AGC" :	"S",
    "AGA" :	"R",
    "AGG" :	"R",
    "GTT" :	"V",
    "GTC" :	"V",
    "GTA" :	"V",
    "GTG" :	"V",
    "GCT" :	"A",
    "GCC" :	"A",
    "GCA" :	"A",
    "GCG" :	"A",
    "GAT" :	"D",
    "GAC" :	"D",
    "GAA" :	"E",
    "GAG" :	"E",
    "GGT" :	"G",
    "GGC" :	"G",
    "GGA" :	"G",
    "GGG" :	"G",
}

# Deals with all IUPAC DNA codes
#
# Code    Description    Complement
#
# A       Adenine        T
# C       Cytosine       G
# G       Guanine        C
# T       Thymine        A
# R       A or G         Y
# Y       C or T         R
# M       C or A         K
# K       T or G         M
# W       T or A         W
# S       C or G         S
# B       C, T or G      V
# V       A, C or G      B
# D       A, T or G      H
# H       A, T or C      D
# N       Any base       N
#
# Others:
# X       Any base       X
# .       Gap            .
# -       Gap            -
complement_table = {
    'A' : 'T',
    'C' : 'G',
    'G' : 'C',
    'T' : 'A',
    'R' : 'Y',
    'Y' : 'R',
    'M' : 'K',
    'K' : 'M',
    'W' : 'W',
    'S' : 'S',
    'B' : 'V',
    'V' : 'B',
    'D' : 'H',
    'H' : 'D',
    'N' : 'N',
    'X' : 'X',
    '.' : '.',
    '-' : '-'
}

def reverse_complement(base_list):
    rcomp = base_list[::-1]
    return ''.join([complement_table[base] for base in rcomp])

    
def base_translate(base):
    """
    if base in codon_table.keys():
        return codon_table[base]
    else:
        return "X"
    """

    try:
        return codon_table[base]
    except KeyError:
        return "X"


def print_width(string, width):
    pos = 0
    while pos < len(string):
        print string[pos:pos+width]
        pos += width


def print_frame_translation(base_list, offset, forward=True, line_length=50):
    print_width(sixframe_translate(base_list, offset, forward), line_length)


def frame_translate(base_list, offset, forward=True):
    if not forward:
        base_list = reverse_complement(base_list)

    aa_list = []
    pos = offset
    while pos <= len(base_list) - 3:
        base = base_list[pos:pos+3]
        aa_list.append(base_translate(base))

        pos += 3

    return "".join(aa_list)


def print_sixframe_translation(base_list, old_header, header_format="translation:%s"):
    forward_translations = [
        ("+1", 0),
        ("+2", 1),
        ("+3", 2)]
    reverse_translations = [
        ("-1", 0),
        ("-2", 1),
        ("-3", 2)]
    base_list = base_list.upper()
    for translation in forward_translations:
        print ">%s %s" % (old_header, header_format % translation[0])
        print_frame_translation(base_list, translation[1], True)
    base_list = reverse_complement(base_list)
    for translation in reverse_translations:
        print ">%s %s" % (old_header, header_format % translation[0])
        print_frame_translation(base_list, translation[1], True)


def main():
    parser = optparse.OptionParser("%prog <input.fasta>")
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)
    input_fasta = open(args[0], 'r')
    bases = ""
    description = None
    for line in input_fasta:
        if len(line.strip()) > 0 and line.strip()[0] == ">":
            if description is not None:
                print_sixframe_translation(bases, description)
            bases = ""
            description = line.strip()[1:]
            continue

        bases += line.strip()

    if len(bases) > 0:
        print_sixframe_translation(bases, description)


if __name__ == "__main__":
    main()
