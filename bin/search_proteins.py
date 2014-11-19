#!/usr/bin/env python

"""
protein_search.py [--clobber|--models [hmmer3|hmmer2|custom]] <input-dir> <output-dir>
"""

import os
import sys
from argparse import ArgumentParser
from mhc.tk import HMMerTk
from massivetools.pipeline import Pipeline
from mhc.data import *


parser = ArgumentParser()
parser.add_argument(
    "-m", "--models", dest="models", default="hmmer3", 
    help="Pre-defined set of hmmer models [hmmer3 (default), hmmer2, custom_hmmer3, custom_hmmer2]")
parser.add_argument("-t", "--threads", dest="threads", type=int, default=NUMBER_THREADS)
parser.add_argument("-x", "--max", action="store_true", dest="max")
parser.add_argument("-c", "--clobber", action="store_true", dest="clobber")
parser.add_argument("input_dir", type=str, help="Input directory (data)")
parser.add_argument("output_dir", type=str, help="Output directory (protein_search)")
args = parser.parse_args()


if args.models=="hmmer3":
    models = hmmer3_models
elif args.models=="custom_hmmer3":
    models = custom_hmmer3_models
elif args.models=="hmmer2":
    models = hmmer2_protein_models
elif args.models=="custom_hmmer2":
    models = custom_hmmer2_protein_models
else:
    print "No such model"
    sys.exit(-1)

if "hmmer3" in args.models and args.max:
    max_option = "--max"
else:
    max_option = ""


hmmer_tk = HMMerTk()

for species in sorted(protein_filenames.keys()):
    protein_filename = os.path.join(args.input_dir, protein_filenames[species])
    print
    print species, protein_filename
    for model_name in sorted(models.keys()):
        model_filename = models[model_name]

        results_dirname = os.path.join(args.output_dir, species, data_species_subdir)
        if not os.path.exists(results_dirname):
            os.makedirs(results_dirname)

        model_basename = os.path.splitext(os.path.basename(model_filename))[0]
        results_basename = "_".join([args.models, model_basename])
#        if max_option:
#            results_basename = "max_" + results_basename
        results_filename = os.path.join(results_dirname, results_basename) + ".txt"
        print " "*3, model_name, results_filename

        if "hmmer2" in args.models:
            jid = hmmer_tk.hmmsearch2(
                input_model=model_filename,
                input_fasta=protein_filename,
                output_filename=results_filename,
                force_out_of_date=args.clobber,
                dependencies=[])
        elif "hmmer3" in args.models:
            table_filename = os.path.join(results_dirname, results_basename) + "_table.txt"
            domtable_filename = os.path.join(results_dirname, results_basename) + "_domtable.txt"
            jid = hmmer_tk.hmmsearch3(
                input_model=model_filename,
                input_fasta=protein_filename,
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
