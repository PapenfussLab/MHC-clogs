#!/usr/bin/env python

"""
Contains classes and functions useful for parsing gene information exports from biomart.org.

If run, parsed output of the given filename will be dumped to standard output.

Author: Chris Davoren
Date: 13/07/2011
"""

import urllib
import urllib2
import sys
import optparse


def zerolen_to_none(string):
    if string is None or len(string) == 0:
        return None
    else:
        return string


def get_with_call(dictionary, key, default, func):
    value = dictionary.get(key, default)
    if value is not None:
        value = func(value)
    return value


def strand_num_to_sign(strand):
    return {"1":"+", "-1":"-", 1:"+", -1:"-"}[strand]


class BiomartGene:
    def __init__(self):
        self.ensembl_gene_id = None
        self.ensembl_transcript_id = None
        self.gene_symbol = None
        self.description = None
        self.chromosome = None
        self.chrom_start = None
        self.chrom_end = None
        self.strand = None

    def __repr__(self):
        format =  "%(ensembl_gene_id)s\t%(ensembl_transcript_id)s\t%(gene_symbol)s\t%(description)s\t%(chromosome)s\t%(chrom_start)s\t%(chrom_end)s\t%(strand)s"
        return format % self.__dict__

    def from_field_dict(self, field_dict):
        self.ensembl_gene_id = field_dict.get("Ensembl Gene ID", None)
        self.ensembl_transcript_id = get_with_call(field_dict, "Ensembl Transcript ID", None, zerolen_to_none)
        self.gene_symbol = get_with_call(field_dict, "Associated Gene Name", None, zerolen_to_none)
        self.description = get_with_call(field_dict, "Description", None, zerolen_to_none)
        self.chromosome = get_with_call(field_dict, "Chromosome Name", None, zerolen_to_none)
        self.chrom_start = get_with_call(field_dict, "Gene Start (bp)", None, int)
        self.chrom_end = get_with_call(field_dict, "Gene End (bp)", None, int)
        self.strand = get_with_call(field_dict, "Strand", None, int)

    def strand_str(self):
        if self.strand == 1:
            return '+'
        elif self.strand == -1:
            return  '-'
        else:
            return '.'

    def __str__(self):
        if self.ensembl_gene_id==None:
            return "No biomart data"
        else:
            return "Ensembl Gene ID: %s, Ensembl Transcript ID: %s, Name: %s, Desc: %s, Location: %s:%d-%d(%s)" % (
                self.ensembl_gene_id,
                self.ensembl_transcript_id,
                self.gene_symbol,
                self.description,
                self.chromosome,
                self.chrom_start,
                self.chrom_end,
                self.strand)

    @staticmethod
    def parse(filename):
        biomart_file = open(filename, 'r')

        # First line contains field headers
        field_names = biomart_file.readline().strip().split("\t")
        genes = {}
        for line in biomart_file:
            # print >> sys.stderr, line
            line = line.strip()
            if len(line) == 0:
                continue
            fields = line.split("\t")
            gene = BiomartGene()
            gene.from_field_dict(dict(zip(field_names, fields)))
            gene.strand = strand_num_to_sign(gene.strand)
            genes[gene.ensembl_gene_id] = gene
        biomart_file.close()
        return genes

    @staticmethod
    def retrieve_dataset(ds_name):
        biomart_xml_request = """
            <?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                <Dataset name = "%s" interface = "default" >
                    <Attribute name = "ensembl_gene_id" />
                    <Attribute name = "ensembl_transcript_id" />
                    <Attribute name = "external_gene_id" />
                    <Attribute name = "description" />
                    <Attribute name = "chromosome_name" />
                    <Attribute name = "start_position" />
                    <Attribute name = "end_position" />
                    <Attribute name = "strand" />
                </Dataset>
            </Query>"""
        req_str = biomart_xml_request.replace("\n", "").strip() % ds_name
#        r = urllib2.urlopen(
#            url="http://www.biomart.org/biomart/martservice",
#            data=urllib.urlencode({'query':req_str}))
        r = urllib2.urlopen(
            url="http://feb2014.archive.ensembl.org/biomart/martservice",
            data=urllib.urlencode({'query':req_str}))
        return r.read()


def main():
    parser = optparse.OptionParser("%prog <biomart_tsv_export.txt>")
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    for ensembl_id, biomart_gene in BiomartGene.load_tsv(args[0]).items():
        print biomart_gene


if __name__ == "__main__":
    main()
