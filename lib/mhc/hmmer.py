"""
hmmer module
"""

import re


hmmer2frame = {0: 1, 1: 2, 2: 3, 3: -1, 4: -2, 5: -3}
frame2hmmer = dict([(v,k) for k,v in hmmer2frame.iteritems()])


class DomainHit:
    attributes = ['domain', 'accession', 'count', 'sStart', 'sEnd', 
        'sCode', 'qStart', 'qEnd', 'qCode', 'score', 'eValue']

    def __init__(self, **kw):
        self.domain = kw.get("domain", None)
        self.accession = kw.get("accession", None)
        self.count = kw.get("count", None)
        self.sStart = kw.get("sStart", None)
        self.sEnd = kw.get("sEnd", None)
        self.sCode = kw.get("sCode", None)
        self.qStart = kw.get("qStart", None)
        self.qEnd = kw.get("qEnd", None)
        self.qCode = kw.get("qCode", None)
        self.score = kw.get("score", None)
        self.eValue = kw.get("eValue", None)

    def __repr__(self):
        format = "%(domain)s\t%(accession)s\t%(count)s\t%(sStart)i\t%(sEnd)i\t%(sCode)s\t%(qStart)i\t%(qEnd)i\t%(qCode)s\t%(score)f\t%(eValue)g"
        return format % self.__dict__

    def __str__(self):
        format = "%(domain)s\t%(accession)s\t%(count)s\t%(sStart)i\t%(sEnd)i\t%(sCode)s\t%(qStart)i\t%(qEnd)i\t%(qCode)s\t%(score)f\t%(eValue)g"
        return format % self.__dict__

    @staticmethod
    def parse_line(line, domain=None):
        tokens = line.split()
        hit = DomainHit()
        hit.domain = domain
        hit.accession = tokens[0]
        hit.count = tokens[1]
        hit.sStart = int(tokens[2])
        hit.sEnd = int(tokens[3])
        hit.sCode = tokens[4]
        hit.qStart = int(tokens[5])
        hit.qEnd = int(tokens[6])
        hit.qCode = tokens[7]
        hit.score = float(tokens[8])
        hit.eValue = float(tokens[9])
        return hit

    @staticmethod
    def parse(input_filename, domain=None):
        """Return an iterator to a HMMer file."""

        input_file = open(input_filename)

        start_regex = re.compile('^Parsed for domains')
        end_regex = re.compile('^Alignments of top-scoring domains')
        abort_regex = re.compile('\[no hits above thresholds\]')

        if not jump_to_match(input_file, start_regex):
            raise Exception('No match found. File may be empty.')

        # 3. Parse domain details
        line = input_file.next()
        line = input_file.next()
        data = []
        for line in input_file:
            line = line.strip()
            if end_regex.match(line) or abort_regex.match(line):
                break
            elif not line:
                continue
            d = DomainHit.parse_line(line, domain=domain)
            data.append(d)
        return data

    def to_genomic(self):
        """Convert block 6 frame coords to genomic, e.g.
        chrom.blockStart-blockEnd:frame aaStart aaEnd or
        chrom:blockStart-blockEnd:frame aaStart aaEnd
        --> chrom,blockStart,blockEnd,gStart,gEnd,strand
        """
        return GenomeDomainHit.from_domain_hit(self)


def jump_to_match(input_file, regex):
    """Jump to regex match in file.
    @param input_file: File object
    @param regex: Compiled regex object
    @return: True if successful, False otherwise
    """
    for line in input_file:
        if regex.match(line):
            return True
    return False


class GenomeDomainHit:
    attributes = ['domain', 'accession', 'count', 'sStart', 'sEnd', 
                  'sCode', 'qStart', 'qEnd', 'qCode', 'score', 'eValue', 'strand']
    def __init__(self, **kw):
        self.domain = kw.get("domain", None)
        self.accession = kw.get("accession", None)
        self.count = kw.get("count", None)
        self.sStart = kw.get("sStart", None)
        self.sEnd = kw.get("sEnd", None)
        self.sCode = kw.get("sCode", None)
        self.qStart = kw.get("qStart", None)
        self.qEnd = kw.get("qEnd", None)
        self.qCode = kw.get("qCode", None)
        self.score = kw.get("score", None)
        self.eValue = kw.get("eValue", None)
        self.strand = kw.get("strand", None)

    def __repr__(self):
        format = "%(domain)s\t%(accession)s\t%(count)s\t%(sStart)i\t%(sEnd)i\t%(sCode)s\t%(qStart)i\t%(qEnd)i\t%(qCode)s\t%(score)f\t%(eValue)g\t%(strand)s"
        return format % self.__dict__

    def __str__(self):
        format = "%(domain)s\t%(accession)s\t%(count)s\t%(sStart)i\t%(sEnd)i\t%(sCode)s\t%(qStart)i\t%(qEnd)i\t%(qCode)s\t%(score)f\t%(eValue)g\t%(strand)s"
        return format % self.__dict__

    @staticmethod
    def parse_line(line, domain=None):
        tokens = line.strip().split("\t")
        hit = GenomeDomainHit()
        hit.domain = tokens[0]
        hit.accession = tokens[1]
        hit.count = tokens[2]
        hit.sStart = int(tokens[3])
        hit.sEnd = int(tokens[4])
        hit.sCode = tokens[5]
        hit.qStart = int(tokens[6])
        hit.qEnd = int(tokens[7])
        hit.qCode = tokens[8]
        hit.score = float(tokens[9])
        hit.eValue = float(tokens[10])
        hit.strand = tokens[11]
        return hit

    @staticmethod
    def parse(filename):
        data = []
        for line in open(filename):
            data.append(GenomeDomainHit.parse_line(line))
        return data

    @staticmethod
    def from_domain_hit(domain_hit):
        genome_domain_hit = GenomeDomainHit(**domain_hit.__dict__)

        regex = re.compile(
            '(?P<chrom>[\w\.]+)([:](?P<block_start>\d+)[-|,](?P<block_end>\d+))?:(?P<frame>[0-5])')
        rs = regex.search(domain_hit.accession)
        d = rs.groupdict()

        block_start = int(d["block_start"])
        block_end = int(d["block_end"])
        L = block_end-block_start+1

        hmmer_frame = int(int(d["frame"]))
        frame = hmmer2frame[hmmer_frame]
        if frame>=0:
            genome_domain_hit.strand = '+'
            g_start = 3*(domain_hit.sStart-1)+(frame-1)+1
            g_end = 3*(domain_hit.sEnd-1)+(frame-1)+3
        else:
            genome_domain_hit.strand = '-'
            g_start = L-(3*(domain_hit.sStart-1)+abs(frame)-1)
            g_end = L-(3*(domain_hit.sEnd-1)+abs(frame)+1)

        genome_domain_hit.accession = d["chrom"]
        genome_domain_hit.sStart = block_start + g_start - 1
        genome_domain_hit.sEnd = block_start + g_end - 1
        if genome_domain_hit.sStart>genome_domain_hit.sEnd:
            genome_domain_hit.sStart,genome_domain_hit.sEnd = genome_domain_hit.sEnd,genome_domain_hit.sStart

        return genome_domain_hit
