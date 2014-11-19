"""
proteintools.py
"""

import re
import string
from mhc.data import *
from mhc.biomart import *
from mhc.hmmer3 import *
from mungo.fasta import FastaFile


class Protein:
    """
    Container for protein data.
    """

    def __init__(self, ensembl_protein_id=None, ensembl_transcript_id=None, 
                 ensembl_gene_id=None, gene_symbol=None, description=None, 
                 header=None, seq=None):
        self.ensembl_protein_id = ensembl_protein_id
        self.ensembl_transcript_id = ensembl_transcript_id
        self.ensembl_gene_id = ensembl_gene_id
        self.gene_symbol = gene_symbol
        self.description = description
        self.header = header
        self.seq = seq
        self.biomart_gene = None

    def __repr__(self):
        format = "%(ensembl_protein_id)s\t%(ensembl_transcript_id)s\t%(ensembl_gene_id)s\t%(gene_symbol)s\t%(description)s" 
        return format % self.__dict__


class ProteinDatabase:
    """
    Database of sequences and protein data.
    Lookup protein/gene data by protein id.
    """

    def __init__(self, species, fasta_filename, biomart_filename):
        self.species = species
        self.proteins = {}
        self.proteins_seq = FastaFile(fasta_filename, indexed=True)
        self.biomart_genes = BiomartGene.parse(biomart_filename)

    def __getitem__(self, ensembl_protein_id):
        # print ensembl_protein_id
        h,s = self.proteins_seq.search(ensembl_protein_id)
        protein = Protein()
        tokens = h.split()
        protein.ensembl_protein_id = tokens[0]
        protein.description = " ".join(tokens[1:])
        protein.ensembl_gene_id = tokens[3].split(":")[1]
        protein.ensembl_transcript_id = tokens[4].split(":")[1]
        try:
            protein.biomart_gene = self.biomart_genes[protein.ensembl_gene_id]
            if not protein.biomart_gene.gene_symbol is None:
                rs = re.search(" \((?P<num>[0-9]+) of [0-9]+\)", protein.biomart_gene.gene_symbol)
                if rs:
                    i,j = rs.start(), rs.end()
                    num = int(rs.groups()[0])
                    protein.biomart_gene.gene_symbol = protein.biomart_gene.gene_symbol[0:i] + string.letters[num-1]
                protein.gene_symbol = "%s_%s" % (
                    species_short[self.species], 
                    protein.biomart_gene.gene_symbol)
            else:
                protein.gene_symbol = protein.ensembl_gene_id
        except KeyError:
            protein.biomart_gene = BiomartGene()
            protein.gene_symbol = protein.ensembl_gene_id
        protein.seq = s
        return protein


class DomainCombination:
    """
    Base class for domain combinations in a single protein
    """

    def __init__(self):
        self.domain_names = []

    def has_domain(self, domain_name):
        return domain_name in self.domain_names

    def sort(self, key=lambda x: -x.overall_score):
        for domain_name in self.domain_names:
            self.__dict__[domain_name].sort(key=key)

    def add_domain(self, domain_name, domain):
        self.domain_names.append(domain_name)
        try:
            self.__dict__[domain_name].append(domain)
        except KeyError:
            self.__dict__[domain_name] = [domain]

    def __repr__(self):
        output = []
        for domain_name in self.domain_names:
            output.append(str(self.__dict__[domain_name]))
        return "\n".join(output) + "\n"


class MHCClassIDomainCombination(DomainCombination):
    """
    The MHCClassIDomainCombination contains domain matches to a single protein.
    Makes testing class I-ness or class II-ness easy.
    """

    def get_score(self):
        try:
            self.sort()
            return self.MHC_I[0].overall_score
        except:
            return 0

    def is_class_I(self, cutoff=1e-5):
        self.sort()
        
        # No class I domain
        if not self.has_domain("MHC_I"): return False

        # Class II domain
        if self.has_domain("MHC_II_beta") and \
            self.MHC_I[0].overall_score<=self.MHC_II_beta[0].overall_score:
            return False

        # Strong class I hit
        if self.MHC_I[0].independent_evalue<=cutoff: return True

        # Weak class I, but no Ig
        if not self.has_domain("C1_set"): return False

        # Weak class I and strong Ig
        weak_hit_plus_ig = self.MHC_I[0].overall_score>0 and \
            self.C1_set[0].independent_evalue<=cutoff
        MHC_I_pos = 0.5*(self.MHC_I[0].sequence_start+self.MHC_I[0].sequence_end)
        C1_set_pos = 0.5*(self.C1_set[0].sequence_start+self.C1_set[0].sequence_end)
        good_position = MHC_I_pos<C1_set_pos
        if weak_hit_plus_ig and good_position:
            return True
        else:
            return False

    def is_class_II(self, cutoff=1e-5):
        # No class II domain
        if not self.has_domain("MHC_II_beta"): return False

        # Class I domain
        if self.has_domain("MHC_I") and \
            self.MHC_II_beta[0].overall_score<=self.MHC_I[0].overall_score:
            return False

        # Strong class II hit
        if self.MHC_II_beta[0].independent_evalue<=cutoff: return True

        # Weak class II, but no Ig
        if not self.has_domain("C1_set"): return False

        # Weak class II and strong Ig
        weak_hit_plus_ig = self.MHC_II_beta[0].overall_score>0 and \
            self.C1_set[0].independent_evalue<=cutoff
        MHC_II_pos = 0.5*(self.MHC_II_beta[0].sequence_start+self.MHC_II_beta[0].sequence_end)
        C1_set_pos = 0.5*(self.C1_set[0].sequence_start+self.C1_set[0].sequence_end)
        good_position = MHC_II_pos<C1_set_pos
        if weak_hit_plus_ig and good_position:
            return True
        else:
            return False


class MHCClassIDomainCombiner:
    """
    Combiner lets you add domains from from multiple searches and collects
    them by protein_id (target_name). Provides a convenient iterator over
    the protein hits.
    """

    def __init__(self):
        self.combinations = {}

    def __getitem__(self, protein_id):
        return self.combinations[protein_id]

    def __iter__(self):
        self.index = -1
        return self

    def next(self):
        self.index += 1
        if self.index>=len(self.combinations): raise StopIteration
        k = self.combinations.keys()[self.index]
        return k

    def add_domains(self, domain_name, domain_filename):
        for domain in DomainHit.parse(domain_filename):
            try:
                self.combinations[domain.target_name].add_domain(domain_name, domain)
            except KeyError:
                self.combinations[domain.target_name] = MHCClassIDomainCombination()
                self.combinations[domain.target_name].add_domain(domain_name, domain)


