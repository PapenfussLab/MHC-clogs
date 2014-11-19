#!/usr/bin/env python

"""
extract_hmmer2_protein_match.py [--models hmmer2|custom_hmmer2] <input dir> <output dir>

Approach:
1. Search with Pfam only & require significant MHC_I or MHC_I and C1, and Score(MHC_I)>Score(MHC_II_beta)
2. Extract 1 representative match per gene
3. Update models and re-search

NB: This version supports hmmer2 only

To do:
- Need to deal with 2 domains (split or otherwise).
"""

import os
from argparse import ArgumentParser
from mhc.data import *
from mhc.proteintools import *
from mungo.fasta import FastaFile


def is_weird(protein):
    try:
        return species=="human" \
            and "HSCHR6_MHC" in protein.biomart_gene.chromosome \
            and not protein.gene_symbol=="MICB" \
            or protein.biomart_gene.chromosome=="HG1322_PATCH"
    except:
        print protein
        sys.exit(-1)


def main(input_dir, output_dir, species, args, output_all_classI_text_file, output_all_classI_fasta_file, output_all_classII_fasta_file):
    # Input data for each species
    fasta_filename = os.path.join(input_dir, protein_filenames[species])
    biomart_filename = os.path.join(input_dir, species, 
        data_species_subdir, "biomart_export_tsv.txt")

    # Setup input filenames
    if args.models=="custom_hmmer2":
        MHC_I_filename = os.path.join(output_dir, species, 
            data_species_subdir, "custom_hmmer2_supfam_mhc_fs.txt")
        MHC_II_beta_filename = os.path.join(output_dir, species, 
            data_species_subdir, "hmmer2_MHC_II_beta_fs.txt")
        Ig_filename = os.path.join(output_dir, species, 
            data_species_subdir, "custom_hmmer2_C1-set_fs.txt")
    else:
        MHC_I_filename = os.path.join(output_dir, species, 
            data_species_subdir, "hmmer2_supfam_mhc_fs.txt")
        MHC_II_beta_filename = os.path.join(output_dir, species, 
            data_species_subdir, "hmmer2_MHC_II_beta_fs.txt")
        Ig_filename = os.path.join(output_dir, species, 
            data_species_subdir, "hmmer2_C1-set_fs.txt")

    # Setup species-specific output filenames
    if args.models=="custom_hmmer2":
        output_spp_text_filename = os.path.join(output_dir, species, 
            data_species_subdir, "custom_hmmer2_protein_search.txt")
        output_spp_classII_fasta_filename = os.path.join(output_dir, species, 
            data_species_subdir, "custom_hmmer2_classII_proteins.fa")
        output_spp_classI_fasta_filename = os.path.join(output_dir, species, 
            data_species_subdir, "custom_hmmer2_classI_proteins.fa")
    else:
        output_spp_text_filename = os.path.join(output_dir, species, 
            data_species_subdir, "hmmer2_protein_search.txt")
        output_spp_classII_fasta_filename = os.path.join(output_dir, species, 
            data_species_subdir, "hmmer2_classII_proteins.fa")
        output_spp_classI_fasta_filename = os.path.join(output_dir, species, 
            data_species_subdir, "hmmer2_classI_proteins.fa")

    # Open output files
    output_spp_text_file = open(output_spp_text_filename, "w")
    output_spp_classII_fasta_file = FastaFile(output_spp_classII_fasta_filename, "w")
    output_spp_classI_fasta_file = FastaFile(output_spp_classI_fasta_filename, "w")

    # Load protein/gene data
    protein_database = ProteinDatabase(species, fasta_filename, biomart_filename)

    # Load domain data
    combiner = MHCClassIDomainCombiner()
    combiner.add_domains("MHC_I", MHC_I_filename)
    combiner.add_domains("C1_set", Ig_filename)
    combiner.add_domains("MHC_II_beta", MHC_II_beta_filename)

    # Collapse isoforms
    results = {}
    for protein_id in combiner:
        protein = protein_database[protein_id]
        protein.combiner = combiner[protein_id]
        key = protein.ensembl_gene_id

        # Chuck out alternative human MHC haplotypes, except MICB
        if is_weird(protein):
            continue
        elif species=="human" and protein.gene_symbol=="MICB":
            key = "MICB"

        try:
            results[key].append(protein)
        except KeyError:
            results[key] = [protein]

    # Sort isoforms by length
    for gene_id in results:
        results[gene_id].sort(key=lambda x: -len(x.seq))

    # Convert to list of lists
    results = results.values()

    # Sort the matches by score
    results.sort(key=lambda x: -x[0].combiner.get_score())

    # Output all class I/II plus a representative (longest) class I for each gene
    classI_gene_count = 0
    format = "\t".join(["%s","%i","%s","%s","%i","%i","%i","%s","%s","%s"])
    fasta_header_format = "%s Position=%s GeneID=%s ProteinID=%s Score=%0.1f E-value=%0.2g"

    # Uniquify gene symbols (corrects naming of bt_BOLA genes)
    symbols = {}
    for result in results:
        protein = result[0]
        try:
            symbols[protein.gene_symbol].append(protein)
            symbols[protein.gene_symbol][-1].gene_symbol += ".%02i" % len(symbols[protein.gene_symbol])
        except KeyError:
            symbols[protein.gene_symbol] = [protein]

    for result in results:
        protein = result[0] ## Keep the longest
        if protein.combiner.is_class_I():
            classI_gene_count += 1
            try:
                position = "%s:%i-%i(%s)" % (protein.biomart_gene.chromosome, protein.biomart_gene.chrom_start, protein.biomart_gene.chrom_end, protein.biomart_gene.strand)
            except:
                print protein.biomart_gene.chromosome, protein.biomart_gene.chrom_start, protein.biomart_gene.chrom_end, protein.biomart_gene.strand
                sys.exit(-1)
            score = protein.combiner.get_score()
            evalue = protein.combiner.get_evalue()

            print >> output_spp_text_file, "###", classI_gene_count
            print >> output_spp_text_file, "Symbol:", protein.gene_symbol
            print >> output_spp_text_file, "Gene id:", protein.ensembl_gene_id
            print >> output_spp_text_file, "Protein id:", protein.ensembl_protein_id
            print >> output_spp_text_file, "# isoforms", len(result)
            print >> output_spp_text_file, "Longest/shortest", len(result[0].seq), len(result[-1].seq)
            print >> output_spp_text_file, "Genomic location", position
            if not protein.combiner.has_domain("MHC_I"):
                protein.combiner.MHC_I = None
            print >> output_spp_text_file, protein.combiner.MHC_I
            if not protein.combiner.has_domain("C1_set"):
                protein.combiner.C1_set = None
            print >> output_spp_text_file, protein.combiner.C1_set
            print >> output_spp_text_file
            print >> output_spp_text_file

            output_all_classI_text_file.write(format % (
                species, 
                classI_gene_count, 
                protein.gene_symbol, 
                protein.ensembl_gene_id, 
                len(result), 
                len(result[0].seq), 
                len(result[-1].seq), 
                position,
                str(protein.combiner.MHC_I).replace("\t", "; "),
                str(protein.combiner.C1_set).replace("\t", "; ")
            )+"\n")

            h = fasta_header_format % (protein.gene_symbol, position, \
                protein.ensembl_gene_id, protein.ensembl_protein_id,
                score, evalue)
            output_all_classI_fasta_file.write(h, protein.seq + "\n")
            output_spp_classI_fasta_file.write(h, protein.seq + "\n")
        elif protein.combiner.is_class_II():
            h = "%s %s %s" % (protein.gene_symbol, protein.ensembl_gene_id, protein.ensembl_protein_id)
            output_all_classII_fasta_file.write(h, protein.seq + "\n")
            output_spp_classII_fasta_file.write(h, protein.seq + "\n")

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("-m", "--models", dest="models", default="custom_hmmer2", 
        help="custom_hmmer2 (default), hmmer2")
    parser.add_argument("-s", "--species", dest="species", default="",
        help="all (default), human, dog, ...")
    parser.add_argument("input_dir", type=str, help='Input directory (data)')
    parser.add_argument("output_dir", type=str, help='Output directory (protein_search)')
    args = parser.parse_args()

    # Setup output file to collect a representative class I from each species
    if args.models=="custom_hmmer2":
        output_all_classI_text_filename = os.path.join(args.output_dir, "results", "custom_hmmer2_classI_proteins.txt")
        output_all_classI_fasta_filename = os.path.join(args.output_dir, "results", "custom_hmmer2_classI_proteins.fa")
        output_all_classII_fasta_filename = os.path.join(args.output_dir, "results", "custom_hmmer2_classII_proteins.fa")
    else:
        output_all_classI_text_filename = os.path.join(args.output_dir, "results", "hmmer2_classI_proteins.txt")
        output_all_classI_fasta_filename = os.path.join(args.output_dir, "results", "hmmer2_classI_proteins.fa")
        output_all_classII_fasta_filename = os.path.join(args.output_dir, "results", "hmmer2_classII_proteins.fa")
    output_all_classI_fasta_file = FastaFile(output_all_classI_fasta_filename, "w")
    output_all_classII_fasta_file = FastaFile(output_all_classII_fasta_filename, "w")

    output_all_classI_text_file = open(output_all_classI_text_filename, "w")
    header = "\t".join(["species", "#", "gene_symbol", "ensembl_gene_id", "# isoforms", "Longest isoform", "Shortest isoform", "Location", "MHC_I", "C1-set"])
    output_all_classI_text_file.write(header + "\n")

    # Can target process to one species
    if not args.species or args.species=="all":
        search_species = sorted(species_short.keys())
    else:
        search_species = [args.species]

    for species in search_species:
        print species
        main(args.input_dir, args.output_dir, species, args, 
            output_all_classI_text_file, output_all_classI_fasta_file, output_all_classII_fasta_file)
