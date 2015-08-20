#!/usr/bin/env python

"""
extract_hmmer2_genome_match.py [options] <data_dir> <input_protein_dir> <output_genome_dir>
Previously DP.py

Author: Tony Papenfuss
Date: Fri Nov 10 21:45:25 EST 2006
"""

import os
import sys
from argparse import ArgumentParser
import numpy
from mungolite.sequence import translate
from mungolite.fasta import FastaFile
from mhc.data import *
from mhc.hmmer import *
from mhc.biomart import *
from mhc.intervals import GenomeIntersector
from mhc.proteintools import Protein
from mhc import blastplus
from mhc.useful import argmax


model = ["S", "alpha1", "alpha2", "Ig", "C-term"]
rescaling = [0, 1, 1, 2, 20]
score_rescale_dict = {
    "S": 0,
    "alpha1": 1,
    "alpha2": 1,
    "Ig": 2,
    "C-term": 20}
GAP_FN_FREE = 5000
GAP_FN_SCALE = 300.0/(20000-GAP_FN_FREE)
MODEL_GAP_PENALTY = 10


def is_weird(accession):
    weird = "HSCHR6_MHC" in accession or "HG1322_PATCH" in accession
    return weird


def load_genomic(filename, reverse_neg_strand=True):
    data = []
    for line in open(filename):
        d = GenomeDomainHit.parse_line(line)
        if is_weird(d.accession): continue
        d.domain = infer_domain_type(d)
        data.append(d)
    data = split_and_sort(data, reverse_neg_strand=reverse_neg_strand)
    return data


def infer_domain_type(d):
    domain_lengths = {"MHC_I_fs": 180, "supfam_MHC_I_fs": 178, "supfam_mhc_fs": 178}
    if d.domain in ["MHC_I", "MHC_I_fs", "supfam_MHC_I_fs", "supfam_mhc_fs"]:
        L = domain_lengths[d.domain]
        bmp = numpy.zeros(L)
        bmp[d.qStart-1:d.qEnd] = 1
        f_alpha1 = float(sum(bmp[0:L/2]))/L
        f_alpha2 = float(sum(bmp[L/2:]))/L
        if f_alpha2 >= f_alpha1:
            name = "alpha2"
        else:
            name = "alpha1"
    elif d.domain in ["Ig", "C1-set_fs", "C1-set_ls", "supfam_Ig_fs", "supfam_ig_fs"]:
        name = "Ig"
    elif d.domain in ["C_term", "MHC_I_C_fs", "MHC_I_C_ls"]:
        name = "C-term"
    else:
        name = ""
    return name


def split_and_sort(data, reverse_neg_strand=True):
    split_data = {}
    for d in data:
        key = (d.accession, d.strand)
        try:
            split_data[key].append(d)
        except KeyError:
            split_data[key] = [d]
    for key in split_data:
        if reverse_neg_strand and key[1] == "-":
            split_data[key].sort(key=lambda d: -d.sStart)
        else:
            split_data[key].sort(key=lambda d: d.sStart)
    return split_data


def match_score(domain, state):
    if domain.domain == state:
        # Match
        score = score_rescale_dict[state]*domain.score
    else:
        # Mis-match
        score = -10000
    return score


def gap_penalty_fn(d):
    free = GAP_FN_FREE
    scale = GAP_FN_SCALE
    if d < free:
        penalty = 0
    else:
        penalty = scale*(d-free)
    return penalty


def build_dp_matrix(data):
    n = len(data) # number of matches in subject
    m = len(model) # number of states in model

    mat = numpy.zeros((n, m), dtype=float)
    gap = numpy.zeros((n, m), dtype=int)
    ptr = numpy.zeros((n, m, 2), dtype=int)

    ## First match

    # Loop over states (domains)
    for j in xrange(1, m):
        if data[0].domain == model[j]:
            mat[0, j] = rescaling[j]*data[0].score
            ptr[0, j, :] = [-1, -1]
        else:
            mat[0, j] = 0
            ptr[0, j, :] = [-1, 0]

    ## Subsequent matches
    for i in xrange(1, n):
        # Loop over states
        for j in xrange(1, m):
            # Distance between features
            dx = abs(data[i].sStart - data[i-1].sEnd)+1

            # Cumulative penalty
            match_gap_penalty = gap_penalty_fn(gap[i-1, j-1]+dx)
            mismatch_gap_penalty = gap_penalty_fn(gap[i-1, j]+dx)

            # Test states
            if j == 1:
                mat[i, j], choice = argmax(
                    mat[i-1, j-1] + match_score(data[i], model[j]),
                    mat[i-1, j] - mismatch_gap_penalty)
            elif j in [2, 3]:
                mat[i, j], choice = argmax(
                    mat[i-1, j-1] + match_score(data[i], model[j]) - match_gap_penalty, # match
                    mat[i-1, j] - mismatch_gap_penalty) # mismatch
            elif j == 4:
                mat[i, j], choice = argmax(
                    mat[i-1, j-1] + match_score(data[i], model[j]) - match_gap_penalty,
                    mat[i-1, j] - mismatch_gap_penalty,
                    mat[i, j-1] - MODEL_GAP_PENALTY)

            if choice == 0:
                ptr[i, j, :] = [-1, -1]
                gap[i, j] = 0
            elif choice == 1:
                ptr[i, j, :] = [-1, 0]
                gap[i, j] = gap[i-1, j]+dx
            elif choice == 2:
                ptr[i, j, :] = [0, -1]
                gap[i, j] = gap[i, j-1]+dx

    return mat, ptr, gap


def traceback(data, mat, ptr, threshold=-150):
    n = len(data)
    m = len(model)

    chains = []
    i = n-1
    j = m-1
    while i > -1:
        if j == m-1 and mat[i, j] < threshold:
            # print i,j,mat[i,j], (ptr[i,j,0],ptr[i,j,1]), data[i]
            i -= 1
            continue
        elif j == m-1 and mat[i, j] >= threshold:
            if ptr[i, j, 0] == -1 and ptr[i, j, 1] == -1:
                chain = [data[i]]
                # print "****",
            elif ptr[i, j, 0] == 0 and ptr[i, j, 1] == -1:
                chain = []
                # print "***",
        elif j < m-1:
            if ptr[i, j, 0] == -1 and ptr[i, j, 1] == -1:
                chain.append(data[i])
                # print "+",

        # print i,j,mat[i,j], (ptr[i,j,0],ptr[i,j,1]), data[i]

        i, j = i+ptr[i, j, 0], j+ptr[i, j, 1]
        if j == 0:
            chain.reverse()
            chains.append(chain)
            j = m-1

    return chains


def output_traceback(data, mat, ptr, gap, stream=sys.stdout, print_header=False):
    n = len(data)
    m = len(model)

    if print_header:
        print >> stream, "\t".join(
            ["Domain", "Score", "Chrom", "Strand", "sStart", "sEnd", "Gap"]
            + model + [""] + model + [""] + model)
    i = 0
    print >> stream, "\t".join([
        data[i].domain, "%0.1f" % data[i].score, data[i].accession,
        data[i].strand, "%i" % data[i].sStart, "%i" % data[i].sEnd, "0"]
        + ["%0.2f" % F for F in mat[i, :]] + [""]
        + ["%i,%i" % tuple(tb) for tb in ptr[i, :]] + [""]
        + ["%i" % x for x in gap[i, :]])

    for i in xrange(1, n):
        print >> stream, "\t".join([
            data[i].domain, "%0.1f" % data[i].score, data[i].accession,
            data[i].strand, "%i" % data[i].sStart, "%i" % data[i].sEnd,
            "%i" % (data[i].sStart-data[i-1].sEnd)]
            + ["%0.2f" % F for F in mat[i, :]] + [""]
            + ["%i,%i" % tuple(tb) for tb in ptr[i, :]] + [""]
            + ["0"] + ["%i" % x for x in gap[i, 1:]])
    print >> stream


def chain_position(chain):
    extrema = []
    for c in chain:
        extrema.append(c.sStart)
        extrema.append(c.sEnd)
    return chain[0].accession, min(extrema), max(extrema), chain[0].strand


def chain_score(chain):
    # Neglects mismatches/gaps !!!
    score = 0.0
    for c in chain:
        score += c.score
    return score


def chain_evalue(chain):
    evalue = 1.0
    for c in chain:
        evalue *= c.eValue
    return evalue


CHAIN_FORMAT = ">chain%0.2i Position=%s Name=%s GeneID=%s ProteinID=%s Score=%0.1f E-value=%0.2g Length=%i Comment=%s"

def output_chains(output_chain_filename, chains):
    output_chain_file = open(output_chain_filename, "w")
    i = 0
    for chain in chains:
        i += 1
        chrom, start, end, strand = chain_position(chain)
        position = "%s:%i-%i(%s)" % (chrom, start, end, strand)
        output_chain_file.write(CHAIN_FORMAT % (i, position, \
            chain[0].name, chain[0].ensembl_gene_id, \
            chain[0].ensembl_protein_id, chain_score(chain), \
            chain_evalue(chain), end-start+1, chain[0].comment)  + "\n")

        for feature in chain:
            output_chain_file.write(str(feature)+"\n")
        output_chain_file.write("\n")


def output_sequences(output_chain_filename, chains, translated=False):
    output_chain_file = open(output_chain_filename, "w")
    i = 0
    for chain in chains:
        i += 1
        chrom, start, end, strand = chain_position(chain)
        position = "%s:%i-%i(%s)" % (chrom, start, end, strand)
        output_chain_file.write(CHAIN_FORMAT % (i, position, \
            chain[0].name, chain[0].ensembl_gene_id, \
            chain[0].ensembl_protein_id, chain_score(chain), \
            chain_evalue(chain), end-start+1, chain[0].comment) + "\n")

        if translated:
            for feature in chain:
                output_chain_file.write(feature.pep + "\n")
        else:
            for feature in chain:
                output_chain_file.write(feature.seq + "\n")

        output_chain_file.write("\n")
    output_chain_file.close()


def output_merged_sequences(output_merged_filename, chains, proteins, genes, translated=True):
    if translated:
        # Output merged peptide sequences
        # Merge sequences
        merged_data = []
        for chain in chains:
            ensembl_gene_id = chain[0].ensembl_gene_id
            try:
                p = copy.copy(proteins[ensembl_gene_id])
                merged_data.append((chain, p))
                del proteins[ensembl_gene_id]
            except:
                merged_data.append((chain, None))

        # Output
        format = "%s Chain=%s Position=%s GeneID=%s ProteinID=%s Score=%0.1f E-value=%0.2g Length=%i Comment=%s"
        output_merged_file = FastaFile(output_merged_filename, "w")
        i = 0
        for chain, p in merged_data:
            i += 1
            chrom, start, end, strand = chain_position(chain)
            position = "%s:%i-%i(%s)" % (chrom, start, end, strand)
            L = end-start+1

            spp = species_short[species]
            if p:
                make_pretty = True
                s = p.seq
            else:
                make_pretty = False
                s = "\n".join([feature.pep for feature in chain])

            chain_name = "chain%0.2i" % i
            h = format % (chain[0].name, chain_name, position, \
                chain[0].ensembl_gene_id, chain[0].ensembl_protein_id, \
                chain_score(chain), chain_evalue(chain), L, chain[0].comment)
            output_merged_file.write(h, s + "\n", make_pretty=make_pretty)

        # format = "%s %s %s"
        for p in proteins.values():
            gene = genes[p.ensembl_gene_id]
            position = "%s:%i-%i(%s)" % (gene.chromosome, gene.chrom_start, gene.chrom_end, gene.strand)
            score = float(p.header.split()[4].split("=")[-1])
            evalue = float(p.header.split()[5].split("=")[-1])
            h = format % (p.gene_symbol, "None", position, p.ensembl_gene_id, \
                p.ensembl_protein_id, score, evalue, 0, "Not detected in genome search")
            output_merged_file.write(h, p.seq + "\n")
        output_merged_file.close()
    else: 
        # Output merged nucleotide sequences
        # Merge sequences
        merged_data = []
        for chain in chains:
            ensembl_gene_id = chain[0].ensembl_gene_id
            try:
                p = copy.copy(proteins[ensembl_gene_id])
                merged_data.append((chain, p))
                del proteins[ensembl_gene_id]
            except:
                merged_data.append((chain, None))

        # Output
        format = "%s Chain=%s Position=%s GeneID=%s ProteinID=%s Score=%0.1f E-value=%0.2g Length=%i Comment=%s"
        output_merged_file = FastaFile(output_merged_filename, "w")
        i = 0
        for chain, p in merged_data:
            i += 1
            chrom, start, end, strand = chain_position(chain)
            position = "%s:%i-%i(%s)" % (chrom, start, end, strand)
            L = end-start+1

            spp = species_short[species]
            if p:
                make_pretty = True
                s = p.seq
            else:
                make_pretty = False
                s = "\n".join([feature.pep for feature in chain])

            chain_name = "chain%0.2i" % i
            h = format % (chain[0].name, chain_name, position, \
                chain[0].ensembl_gene_id, chain[0].ensembl_protein_id, \
                chain_score(chain), chain_evalue(chain), L, chain[0].comment)
            output_merged_file.write(h, s + "\n", make_pretty=make_pretty)

        # format = "%s %s %s"
        for p in proteins.values():
            gene = genes[p.ensembl_gene_id]
            position = "%s:%i-%i(%s)" % (gene.chromosome, gene.chrom_start, gene.chrom_end, gene.strand)
            score = float(p.header.split()[4].split("=")[-1])
            evalue = float(p.header.split()[5].split("=")[-1])
            h = format % (p.gene_symbol, "None", position, p.ensembl_gene_id, \
                p.ensembl_protein_id, score, evalue, 0, "Not detected in genome search")
            output_merged_file.write(h, p.seq + "\n")
        output_merged_file.close()


def add_sequence(blast_filename, chain):
    seq = []
    for feature in chain:
        h, s = blastplus.get_sequence(
            blast_filename, "lcl|" + feature.accession, \
            start=feature.sStart, end=feature.sEnd, strand=feature.strand)
        feature.seq = s
        feature.pep = translate(s)
    return chain


def annotate_chains(species, chains, biomart_filename, blast_filename,
                    input_protein_filename, add_seq=True):
    spp = species_short[species]

    # Load biomart data & construct a genome intersector from biomart data
    genes = BiomartGene.parse(biomart_filename)
    gi = GenomeIntersector()
    for gene in genes.values():
        gi.add(gene.chromosome, gene.strand, gene.chrom_start, gene.chrom_end, gene)

    # Load protein search results
    proteins = {}
    for h, s in FastaFile(input_protein_filename):
        tokens = h.split()
        ensembl_gene_id = tokens[2].split("=")[-1]
        proteins[ensembl_gene_id] = Protein(
            gene_symbol=tokens[0],
            ensembl_gene_id=ensembl_gene_id,
            ensembl_protein_id=tokens[3].split("=")[-1],
            header=h,
            seq=s)

    # Annotate each chain
    i = 0
    for chain in chains:
        i += 1
        name = "%s_chain%0.2i" % (spp, i)

        # Add genomic/peptide sequences
        if add_seq:
            chain = add_sequence(blast_filename, chain)
        else:
            for feature in chain:
                feature.seq = ""
                feature.pep = ""

        # Annotate using biomart and protein search results
        chrom, start, end, strand = chain_position(chain)
        rs = gi.find(chrom, strand, start, end)
        if len(rs) == 0:
            # Unannotated region
            ensembl_gene_id = None
            ensembl_protein_id = None
            comment = "No overlapping annotations"
        elif len(rs) == 1:
            # Single overlapping annotation
            r = rs[0]
            if not r.value.gene_symbol is None:
                name = "%s_%s" % (spp, r.value.gene_symbol)

            ensembl_gene_id = r.value.ensembl_gene_id
            if r.value.ensembl_gene_id in proteins:
                ensembl_protein_id = proteins[ensembl_gene_id].ensembl_protein_id
                comment = "Overlaps single annotation & in proteins"
            else:
                ensembl_protein_id = None
                comment = "Overlaps single annotation; not in proteins"
        elif len(rs) > 1:
            # Multiple overlapping annotations
            # Rank hits based on proportion of overlap
            for r in rs:
                position = [start, end, r.value.chrom_start, r.value.chrom_end]
                position.sort()
                r.value.overlap = float(position[2]-position[1])/float(position[3]-position[0])
            rs.sort(key=lambda x: -x.value.overlap)

            # Look for overlap with corresponding protein search match
            for r in rs:
                if r.value.ensembl_gene_id in proteins:
                    match = True
                    if not r.value.gene_symbol is None:
                        name = "%s_%s" % (spp, r.value.gene_symbol)
                    ensembl_gene_id = r.value.ensembl_gene_id
                    ensembl_protein_id = proteins[ensembl_gene_id].ensembl_protein_id
                    comment = "Overlaps multiple annotations; one in proteins"
                    break
                else:
                    match = False

            if not match:
                # Overlaps multiple; none in proteins; use annotation with highest proportion of overlap
                if not r.value.gene_symbol is None:
                    name = "%s_%s" % (spp, rs[0].value.gene_symbol)
                ensembl_gene_id = rs[0].value.ensembl_gene_id
                ensembl_protein_id = None
                comment = "Overlaps multiple annotations; none in proteins: %s" \
                    % (";".join(["(%s, %s)" % (r.value.gene_symbol, r.value.ensembl_gene_id) for r in rs]))

        for feature in chain:
            feature.name = name
            feature.ensembl_gene_id = ensembl_gene_id
            feature.ensembl_protein_id = ensembl_protein_id
            feature.comment = comment

    return chains, proteins, genes


def main(species, input_domain_filename, input_protein_filename,
         biomart_filename, blast_filename, output_chain_filename,
         output_fasta_filename, output_peptide_filename,
         output_merged_filename, output_merged_peptide_filename,
         traceback_filename, debug=False):

    data = load_genomic(input_domain_filename)

    # Align model
    if debug:
        tb_file = open(traceback_filename, "w")
    else:
        tb_file = None
    first_time = True

    chains = []
    for key in data:
        # if debug: print "   ", key
        mat, ptr, gap = build_dp_matrix(data[key])
        if debug:
            output_traceback(data[key], mat, ptr, gap, stream=tb_file, print_header=first_time)
            first_time = False

        for chain in traceback(data[key], mat, ptr):
            if chain_score(chain) > 10 and chain_evalue(chain) < 0.1:
                chains.append(chain)

    if debug: tb_file.close()
    chains.sort(key=lambda chain: -chain_score(chain))


# To do:
# 2. Look at annotate chains (TAP1 issue)
# 1. Look at bug in merging


    # Annotate chains
    chains, proteins, genes = annotate_chains(
        species, chains, biomart_filename, blast_filename, input_protein_filename)

    # Output
    output_chains(output_chain_filename, chains)
    output_sequences(output_peptide_filename, chains, translated=True)
    output_sequences(output_fasta_filename, chains, translated=False)
    output_merged_sequences(output_merged_peptide_filename, chains, proteins, genes, translated=True)
    output_merged_sequences(output_merged_peptide_filename, chains, proteins, genes, translated=False)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-m", "--models", dest="models", default="custom_hmmer2",
        help="Pre-defined set of hmmer models [custom_hmmer2 (default), hmmer2]")
    parser.add_argument(
        "-s", "--species", dest="species", default="",
        help="Species to search [all (default), human, dog, ...]")
    parser.add_argument(
        "-d", "--debug", dest="debug", action="store_true", default=False, 
        help="Output traceback information")
    parser.add_argument("input_data_dir", type=str, help="Input directory (data)")
    parser.add_argument("input_protein_dir", type=str, help="Input directory (data)")
    parser.add_argument("output_genome_dir", type=str, help="Output directory (genome_search)")
    args = parser.parse_args()

    if not args.species or args.species == "all":
        search_species = sorted(species_short.keys())
    else:
        search_species = [args.species]

    for species in search_species:
        print species
        biomart_filename = os.path.join(args.input_data_dir, biomart_filenames[species])
        blast_filename = os.path.join(args.input_data_dir, genome_filenames[species])
        if args.models == "hmmer2":
            input_domain_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "hmmer2_genomic.txt")
            input_protein_filename = os.path.join(args.input_protein_dir, species, data_species_subdir, "hmmer2_classI_proteins.fa")
            output_chain_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "hmmer2_chains.txt")
            output_fasta_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "hmmer2_chains.fa")
            output_peptide_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "hmmer2_chains_pep.fa")
            output_merged_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "hmmer2_merged_classI.fa")
            output_merged_peptide_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "hmmer2_merged_classI_proteins.fa")
            traceback_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "hmmer2_traceback.txt")
        elif args.models == "custom_hmmer2":
            input_domain_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "custom_hmmer2_genomic.txt")
            input_protein_filename = os.path.join(args.input_protein_dir, species, data_species_subdir, "custom_hmmer2_classI_proteins.fa")
            output_chain_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "custom_hmmer2_chains.txt")
            output_fasta_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "custom_hmmer2_chains.fa")
            output_peptide_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "custom_hmmer2_chains_pep.fa")
            output_merged_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "custom_hmmer2_merged_classI.fa")
            output_merged_peptide_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "custom_hmmer2_merged_classI_proteins.fa")
            traceback_filename = os.path.join(args.output_genome_dir, species, data_species_subdir, "custom_hmmer2_traceback.txt")
        else:
            print "No such model"
            sys.exit(-1)

        main(species, input_domain_filename, input_protein_filename, 
            biomart_filename, blast_filename, output_chain_filename, 
            output_fasta_filename, output_peptide_filename, 
            output_merged_filename, output_merged_peptide_filename, 
            traceback_filename, debug=args.debug)
