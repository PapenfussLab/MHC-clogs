#!/usr/bin/env python

import os
import sys
from argparse import ArgumentParser
from mhc import biomart
from mhc.data import *


ensembl_ftp_address = 'ftp://ftp.ensembl.org/pub/'

biomart_gene_datasets = {
    'c_intestinalis' : 'cintestinalis_gene_ensembl',
    'chicken'        : 'ggallus_gene_ensembl',
    'cow'            : 'btaurus_gene_ensembl',
    'devil'          : 'sharrisii_gene_ensembl',
    'dog'            : 'cfamiliaris_gene_ensembl',
    'fruitfly'       : 'dmelanogaster_gene_ensembl',
    'human'          : 'hsapiens_gene_ensembl',
    'lamprey'        : 'pmarinus_gene_ensembl',
    'lizard'         : 'acarolinensis_gene_ensembl', 
    'mouse'          : 'mmusculus_gene_ensembl',
    'opossum'        : 'mdomestica_gene_ensembl',
    'platypus'       : 'oanatinus_gene_ensembl',
    's_cerevisiae'   : 'scerevisiae_gene_ensembl',
    'tetraodon'      : 'tnigroviridis_gene_ensembl',
    'turkey'         : 'mgallopavo_gene_ensembl',
    'wallaby'        : 'meugenii_gene_ensembl',
    'x_tropicalis'   : 'xtropicalis_gene_ensembl', 
    'zebrafinch'     : 'tguttata_gene_ensembl',
    'zebrafish'      : 'drerio_gene_ensembl',
}



def main(output_dir, force=False):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # Retrieve genomes
    print >> sys.stderr, "Retrieving reference genomes using wget (large files; may take a long time)..."
    for species in sorted(ensembl_genome_filenames.keys()):
        if not os.path.exists(os.path.join(output_dir, species)):
            os.mkdir(os.path.join(output_dir, species))
        genome_filename = os.path.join(output_dir, genome_filenames[species])
        genome_uncompressed_filename = genome_filename[:-3]
        print >> sys.stderr, "\t%s..." % species,
        sys.stderr.flush()
        if (os.path.exists(genome_filename) and os.path.getsize(genome_filename) > 0) or (
            os.path.exists(genome_uncompressed_filename) and os.path.getsize(genome_uncompressed_filename) > 0):
            print >> sys.stderr, "already exists."
            continue
        print >> sys.stderr, "downloading..."
        command = "wget --directory-prefix=%s %s" % (os.path.dirname(genome_filename), ensembl_ftp_address+ensembl_genome_filenames[species])
        print >> sys.stderr, "\t\t%s" % command
        es = os.system(command)
        if es != 0:
            print >> sys.stderr, "Errors encountered (see above), exiting."
            sys.exit(1)
    
    # Retrieve annotated protein data
    print >> sys.stderr, "Retrieving annotated protein data using wget..."
    for species in sorted(ensembl_genome_filenames.keys()):
        print >> sys.stderr, "\t%s..." % species,
        protein_filename = ensembl_protein_filenames[species]
        protein_full_path = os.path.join(output_dir, protein_filenames[species])
        protein_uncompressed = protein_full_path[:-3]
        if (os.path.exists(protein_full_path) and os.path.getsize(protein_full_path) > 0) or (
            os.path.exists(protein_uncompressed) and os.path.getsize(protein_uncompressed) > 0):
            print >> sys.stderr, "already exists."
            continue
        print >> sys.stderr, "downloading..."
        command = "wget --directory-prefix=%s %s" % (os.path.dirname(protein_full_path), ensembl_ftp_address+protein_filename)
        print >> sys.stderr, "\t\t%s" % command
        es = os.system(command)
        if es != 0:
            print >> sys.stderr, "Errors encountered (see above), exiting."
            sys.exit(1)

    # Retrieve annotated transcript data
    print >> sys.stderr, "Retrieving annotated transcript data using wget..."
    for species in sorted(ensembl_genome_filenames.keys()):
        print >> sys.stderr, "\t%s..." % species,
        transcript_filename = ensembl_transcript_filenames[species]
        transcript_full_path = os.path.join(output_dir, transcript_filenames[species])
        transcript_uncompressed = transcript_full_path[:-3]
        if (os.path.exists(transcript_full_path) and os.path.getsize(transcript_full_path) > 0) or (
            os.path.exists(transcript_uncompressed) and os.path.getsize(transcript_uncompressed) > 0):
            print >> sys.stderr, "already exists."
            continue
        print >> sys.stderr, "downloading..."
        command = "wget --directory-prefix=%s %s" % (os.path.dirname(transcript_full_path), ensembl_ftp_address+transcript_filename)
        print >> sys.stderr, "\t\t%s" % command
        es = os.system(command)
        if es != 0:
            print >> sys.stderr, "Errors encountered (see above), exiting."
            sys.exit(1)
    
    # Decompress all downloaded FASTA files
    print >> sys.stderr, "Decompressing downloaded genomes..."
    try:
        os.system("find . -name '*.fa.gz' -print -execdir gunzip {} ';'")
    except:
        pass
    
    # Build blast databases
    for species in sorted(ensembl_genome_filenames.keys()):
        genome_uncompressed_filename = os.path.join(output_dir, genome_filenames[species])
        print "!!!", genome_uncompressed_filename
        if not os.path.exists(genome_uncompressed_filename + ".nin"):
            os.system("makeblastdb -max_file_sz 10GB -in %s -dbtype nucl -parse_seqids" % genome_uncompressed_filename)
    
    # Retrieve Biomart gene datasets
    print >> sys.stderr, "Retrieving Ensembl gene datasets..."
    for species in sorted(biomart_gene_datasets.keys()):
        print >> sys.stderr, "\t%s..." % species,
        sys.stderr.flush()
        biomart_filename = os.path.join(output_dir, biomart_filenames[species])
        if os.path.exists(biomart_filename) and os.path.getsize(biomart_filename) > 0 and not force:
            print >> sys.stderr, "already exists."
            continue
        if not os.path.exists(os.path.join(output_dir, species)):
            os.mkdir(os.path.join(output_dir, species))
        if not os.path.exists(os.path.join(output_dir, species, data_species_subdir)):
            os.mkdir(os.path.join(output_dir, species, data_species_subdir))
        fp = open(biomart_filename, 'w')
        print >> fp, biomart.BiomartGene.retrieve_dataset(biomart_gene_datasets[species])
        fp.close()
        print >> sys.stderr, "done."
    print >> sys.stderr, "Finished."


if __name__ == "__main__":
    parser = ArgumentParser("%prog [-f] <output-dir>")
    parser.add_argument("-f", "--force", dest="force", default=False, action="store_true", help="Retrieve datasets even if the file already exists.")
    parser.add_argument("output_dir", type=str, help="Output directory (data)")
    args = parser.parse_args()
    main(args.output_dir, args.force)
