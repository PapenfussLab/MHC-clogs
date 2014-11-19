"""
hmmer3
"""

from useful import *


class DomainHit:
    def __init__(self, **kw):
        self.target_name = kw.get("target_name", None)
        self.target_accession = kw.get("target_accession", None)
        self.target_length = kw.get("target_length", None)
        self.query_name = kw.get("query_name", None)
        self.query_accession = kw.get("query_accession", None)
        self.query_length = kw.get("query_length", None)
        self.overall_evalue = kw.get("overall_evalue", None)
        self.overall_score = kw.get("overall_score", None)
        self.overall_bias = kw.get("overall_bias", None)
        self.domain_number = kw.get("domain_number", None)
        self.total_sequence_domains = kw.get("total_sequence_domains", None)
        self.conditional_evalue = kw.get("conditional_evalue", None)
        self.independent_evalue = kw.get("independent_evalue", None)
        self.domain_score = kw.get("domain_score", None)
        self.domain_bias = kw.get("domain_bias", None)

        # HMM model alignment positions
        self.profile_start = kw.get("profile_start", None)
        self.profile_end = kw.get("profile_end", None)

        # Protein sequence alignment positions
        self.sequence_start = kw.get("sequence_start", None)
        self.sequence_end = kw.get("sequence_end", None)

        self.envelope_start = kw.get("envelope_start", None)
        self.envelope_end = kw.get("envelope_end", None)
        self.aligned_residue_probability = kw.get("aligned_residue_probability", None)
        self.description = kw.get("description", None)

    def __repr__(self):
        format = "%-20s %-10s %5d %-20s %-10s %5d %9g %6g %5g %3d %3d %9g %9g %6g %5g %5d %5d %5d %5d %5d %5d %4f %s"
        values = [self.target_name, self.target_accession, self.target_length,
            self.query_name, self.query_accession, self.query_length,
            self.overall_evalue, self.overall_score, self.overall_bias,
            self.domain_number, self.total_sequence_domains, self.conditional_evalue,
            self.independent_evalue, self.domain_score, self.domain_bias,
            self.profile_start, self.profile_end, self.sequence_start,
            self.sequence_end, self.envelope_start, self.envelope_end,
            self.aligned_residue_probability, self.description]

        output = []
        for f,x in zip(format.split(" "), values):
            try:
                output.append(f % x)
            except:
                output.append(".")
        return "\t".join(output)

    @staticmethod
    def parse(filename):
        table_file = open(filename, 'r')
        domains = []
        for line in table_file:
            if line.strip()[0]=="#" or not line.strip():
                continue

            fields = remove_blank_fields(line.strip().split(' '), 23)
            if len(fields) < 23:
                continue

            domain = DomainHit()
            domain.target_name = fields[0]
            domain.target_accession = fields[1]
            domain.target_length = int(fields[2])
            domain.query_name = fields[3]
            domain.query_accession = fields[4]
            domain.query_length = int(fields[5])
            domain.overall_evalue = float(fields[6])
            domain.overall_score = float(fields[7])
            domain.overall_bias = float(fields[8])
            domain.domain_number = int(fields[9])
            domain.total_sequence_domains = int(fields[10])
            domain.conditional_evalue = float(fields[11])
            domain.independent_evalue = float(fields[12])
            domain.domain_score = float(fields[13])
            domain.domain_bias = float(fields[14])
            domain.profile_start = int(fields[15])
            domain.profile_end = int(fields[16])
            domain.sequence_start = int(fields[17])
            domain.sequence_end = int(fields[18])
            domain.envelope_start = int(fields[19])
            domain.envelope_end = int(fields[20])
            domain.aligned_residue_probability = float(fields[21])

            if len(fields) > 22:
                domain.description = " ".join(fields[22:])

            domains.append(domain)

        table_file.close()
        return domains


def sixframe_to_genome(start, end, frame, length):
    """
    Converts 6-frame amino acid sequence coordinates to genomic coordinates.  All coordinates are 0-based.
    @param start: amino acid start coordinate
    @param end: amino acid end coordinator (non-inclusive)
    @param frame: reading frame (one of +1, +2, +3, -1, -2, -3)
    @param length: reference genome nucleotide sequence length - this parameter is required for negative frame calculations only
    @return: (start, end) tuple of corresponding genome coordinates
    """
    if frame >= 0:
        start = 3*start+(frame-1)
        end = 3*end+(frame-1)
    else:
        start = length-(3*start)-(frame-1)
        end = length-(3*end)-(frame-1)
    return (start, end)
