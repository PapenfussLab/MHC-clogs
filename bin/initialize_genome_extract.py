#!/usr/bin/env python

"""
initialize_genome_extract.py
"""

import os
import sys
import glob
import copy
from argparse import ArgumentParser
from srt.intervals import GenomeIntersector
from mhc.data import *
from mhc.hmmer import *


def is_not_class_II(classI, classII_intersector):
    rs = classII_intersector.find(classI.accession, classI.sStart, classI.sEnd)
    if rs:
        rs.sort(key=lambda x: -x.value.score)
        if classI.score<rs[0].value.score:
            return False
    return True


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("-m", "--method", dest="method", default="hmmer2", help="hmmer2 (default), custom_hmmer2")
    parser.add_argument("-s", "--species", dest="species", default="", help="all (default), human, dog, ...")
    parser.add_argument("work_dir", type=str, help="Working directory (e.g. output_genomes)")
    args = parser.parse_args()

    work_dir = args.work_dir

    if not args.species or args.species=="all":
        search_species = sorted(species_short.keys())
    else:
        search_species = [args.species]

    for species in search_species:
        # Setup input filenames
        if args.method=="hmmer2":
            MHC_I_filename = os.path.join(work_dir, species, data_species_subdir, "hmmer2_supfam_mhc_fs.txt")
            MHC_II_beta_filename = os.path.join(work_dir, species, data_species_subdir, "hmmer2_MHC_II_beta_fs.txt")
            Ig_filename = os.path.join(work_dir, species, data_species_subdir, "hmmer2_C1-set_fs.txt")
            C_term_filename = os.path.join(work_dir, species, data_species_subdir, "hmmer2_MHC_I_C_fs.txt")
            output_filename = os.path.join(work_dir, species, data_species_subdir, "hmmer2_genomic.txt")
        elif args.method=="custom_hmmer2":
            MHC_I_filename = os.path.join(work_dir, species, data_species_subdir, "custom_hmmer2_supfam_mhc_fs.txt")
            MHC_II_beta_filename = os.path.join(work_dir, species, data_species_subdir, "hmmer2_MHC_II_beta_fs.txt")
            Ig_filename = os.path.join(work_dir, species, data_species_subdir, "custom_hmmer2_C1-set_fs.txt")
            C_term_filename = os.path.join(work_dir, species, data_species_subdir, "hmmer2_MHC_I_C_fs.txt")
            output_filename = os.path.join(work_dir, species, data_species_subdir, "custom_hmmer2_genomic.txt")

        # Load domains
        if "hmmer2" in args.method:
            MHC_I_hits = DomainHit.parse(MHC_I_filename, domain="supfam_mhc_fs")
            MHC_II_beta_hits = DomainHit.parse(MHC_II_beta_filename, domain="MHC_II_beta_fs")
            Ig_hits = DomainHit.parse(Ig_filename, domain="C1-set_fs")
            C_term_hits = DomainHit.parse(C_term_filename, domain="MHC_I_C_fs")

        MHC_I_hits = [d.to_genomic() for d in MHC_I_hits]
        MHC_II_beta_hits = [d.to_genomic() for d in MHC_II_beta_hits]
        Ig_hits = [d.to_genomic() for d in Ig_hits]
        C_term_hits = [d.to_genomic() for d in C_term_hits]

        print "%s\t%i\t%i\t%i\t%i" % (species, len(MHC_I_hits), len(Ig_hits), len(C_term_hits), len(MHC_II_beta_hits))

        # Filter class I/II overlaps hits
        classII_intersector = GenomeIntersector()
        for d in MHC_II_beta_hits:
            classII_intersector.add(d.accession, d.sStart, d.sEnd, d)
        MHC_I_hits = [d for d in MHC_I_hits if is_not_class_II(d, classII_intersector)]

        output_file = open(output_filename, 'w')
        for d in MHC_I_hits+Ig_hits+C_term_hits:
            output_file.write(str(d)+"\n")
        output_file.close()
