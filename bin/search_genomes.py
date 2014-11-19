#!/usr/bin/env python

"""
genome_search.py <input-dir> <output-dir>
"""

import os
import sys
from argparse import ArgumentParser
from mhc.data import *
from mhc.tk import HMMerTk
from massivetools.pipeline import Pipeline


parser = ArgumentParser()
parser.add_argument("-s", "--species", dest="species", type=str, default="all", help="all (default), human, dog, ...")
parser.add_argument("-m", "--models", dest="models", type=str, default="hmmer3", 
     help="Pre-defined set of hmmer models [hmmer3 (default), hmmer2, custom_hmmer3, custom_hmmer2, exon_hmmer2]")
parser.add_argument("-t", "--threads", dest="threads", type=int, default=8,
    help="Number of parallel threads (default 8)")
parser.add_argument("-c", "--clobber", action="store_true", dest="clobber")
parser.add_argument("-x", "--max", action="store_true", dest="max")
parser.add_argument("input_dir", type=str, help="Input directory (data)")
parser.add_argument("output_dir", type=str, help="Output directory (genome_search)")
args = parser.parse_args()


if args.models=="hmmer2":
    models = hmmer2_models
elif args.models=="hmmer3":
    models = hmmer3_models
elif args.models=="custom_hmmer2":
    models = custom_hmmer2_models
elif args.models=="custom_hmmer3":
    models = custom_hmmer3_models
elif args.models=="exon_hmmer2":
    models = exon_hmmer2_models
else:
    print "No such model"
    sys.exit(-1)

if "hmmer3" in args.models and args.max:
    max_option = "--max"
else:
    max_option = ""

if not args.species or args.species=="all":
    search_species = sorted(species_short.keys())
else:
    search_species = [args.species]


hmmer_tk = HMMerTk()
for species in search_species:
    genome_filename = os.path.join(args.input_dir, genome_filenames[species])
    translated_filename = genome_filename.replace(".dna.toplevel.fa", ".6_frames.fa")
    print
    print species, genome_filename, translated_filename
    for model_name in sorted(models.keys()):
        model_filename = models[model_name]

        results_dirname = os.path.join(args.output_dir, species, data_species_subdir)
        if not os.path.exists(results_dirname):
            os.makedirs(results_dirname)

        model_basename = os.path.splitext(os.path.basename(model_filename))[0]
        fasta_basename = os.path.splitext(os.path.basename(genome_filename))[0]
        results_basename = "_".join([args.models, model_basename])
        results_filename = os.path.join(results_dirname, results_basename) + ".txt"
        print " "*3, model_name, results_filename

        if "hmmer2" in args.models:
            jid = hmmer_tk.hmmsearch2(
                input_model=model_filename,
                input_fasta=translated_filename,
                output_filename=results_filename,
                force_out_of_date=args.clobber,
                dependencies=[])
        elif "hmmer3" in args.models:
            table_filename = os.path.join(results_dirname, results_basename) + "_table.txt"
            domtable_filename = os.path.join(results_dirname, results_basename) + "_domtable.txt"
            jid = hmmer_tk.hmmsearch3(
                input_model=model_filename,
                input_fasta=translated_filename,
                output_filename=results_filename,
                output_table=table_filename,
                output_domtable=domtable_filename,
                max_option=max_option,
                force_out_of_date=args.clobber,
                dependencies=[])
        else:
            print "No such model"
            sys.exit(-1)

Pipeline.run(max_parallel=args.threads)
