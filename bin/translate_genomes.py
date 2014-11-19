#!/usr/bin/env python

"""
translate_genomes.py <input-dir>
"""

import os
from argparse import ArgumentParser
from mhc.data import *
from mhc.tk import HMMerTk
from massivetools.pipeline import Pipeline


parser = ArgumentParser()
parser.add_argument("input_dir", type=str, help="Input directory (data)")
args = parser.parse_args()


hmmer_tk = HMMerTk()

for species in sorted(genome_filenames.keys()):
    genome_filename = os.path.join(args.input_dir, genome_filenames[species])
    dirname = os.path.split(genome_filename)[0]
    block_genome_filename = genome_filename.replace(".dna.toplevel.fa", ".block.fa")
    translated_filename = genome_filename.replace(".dna.toplevel.fa", ".6_frames.fa")
    
    jid1 = hmmer_tk.split(
        input_genome_filename=genome_filename,
        output_block_genome_filename=block_genome_filename,
        dependencies=[])
    
    jid2 = hmmer_tk.translate(
        input_genome_filename=block_genome_filename,
        output_6frame_filename=translated_filename,
        dependencies=[jid1])
    
    jid3 = hmmer_tk.remove(
        input_filename=block_genome_filename,
        dependencies=[jid2])

Pipeline.run(max_parallel=4)
