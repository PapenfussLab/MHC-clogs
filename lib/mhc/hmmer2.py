#!/usr/bin/env python

import os
import re
import sys


class SequenceHit:
    def __init__(self):
        self.sequence = None
        self.description = None
        self.score = None
        self.evalue = None
        self.num_domains = None

    @staticmethod
    def parse(filename):
        table_file = open(filename, 'r')
        
        sequence_hits = []
        parse_state = 'NOTFOUND'
        sequence_field_length = None
        description_field_length = None
        for line in table_file:
            line = line.strip()
            if parse_state == 'NOTFOUND':
                match = re.match('(Sequence.*) (Description.*) Score.*E-value.*N', line)
                if match is not None:
                    sequence_field_length = len(match.group(1))
                    description_field_length = len(match.group(2))
                    parse_state = 'HEADERFOUND'
                continue
            elif parse_state == 'HEADERFOUND':
                parse_state = 'FOUND'
                continue
            elif parse_state == 'FOUND':
                if len(line) == 0:
                    break
                sequence_hits.append(SequenceHit.parse_line(line, sequence_field_length, description_field_length))
                continue
        return sequence_hits

    @staticmethod
    def parse_line(line_string, name_field_width=None, desc_field_width=52):
        line_string = line_string.strip()
        if name_field_width is None:
            match = re.match('(\S*)(\s*)')
            name_field_width = len(match.group(1)) + len(match.group(2)) - 1
        sequence_hit = SequenceHit()
        sequence_hit.sequence = line_string[:name_field_width]
        sequence_hit.description = line_string[name_field_width+1:name_field_width+1+desc_field_width]
        fields = line_string[name_field_width+desc_field_width+2:].split()
        sequence_hit.score = float(fields[0])
        sequence_hit.evalue = float(fields[1])
        sequence_hit.num_domains = int(fields[2])
        return sequence_hit

    def __str__(self):
        return "%s %g %d %s" % (self.sequence, self.evalue, self.num_domains, self.description)


class DomainHit:
    def __init__(self, **kw):
        self.sequence = kw.get("sequence", None)
        self.domain_num = kw.get("domain_num", None)
        self.domain_total = kw.get("domain_total", None)
        self.sequence_from = kw.get("sequence_from", None)
        self.sequence_to = kw.get("sequence_to", None)
        self.sequence_from_start = kw.get("sequence_from_start", None)
        self.sequence_to_end = kw.get("sequence_to_end", None)
        self.hmm_from = kw.get("hmm_from", None)
        self.hmm_to = kw.get("hmm_to", None)
        self.hmm_from_start = kw.get("hmm_from_start", None)
        self.hmm_to_end = kw.get("hmm_to_end", None)
        self.score = kw.get("score", None)
        self.evalue = kw.get("evalue", None)

    def __setattr__(self, k, v):
        try:
            k = {
                "target_name": "sequence",
                "sequence_start": "sequence_from",
                "sequence_end": "sequence_to",
            }[k]
        except:
            pass
        self.__dict__[k] = v

    def __getattr__(self, k):
        if k=="target_name":
            return self.sequence
        elif k=="sequence_start":
            return self.sequence_from
        elif k=="sequence_end":
            return self.sequence_to
        elif k=="overall_score":
            return self.score
        elif k=="independent_evalue":
            return self.evalue
        else:
            return self.__dict__[k]

    def __repr__(self):
        format = "%(sequence)s\t%(sequence_from)s\t%(sequence_to)s\t%(score)s\t%(evalue)s"
        return format % self.__dict__

    def sequence_str(self):
        seq_str = ""
        seq_str += ['.','['][self.sequence_from_start]
        seq_str += ['.',']'][self.sequence_to_end]
        return seq_str

    def hmm_str(self):
        h_str = ""
        h_str += ['.','['][self.hmm_from_start]
        h_str += ['.',']'][self.hmm_to_end]
        return h_str

    @staticmethod
    def parse(filename):
        table_file = open(filename, 'r')
        
        domain_hits = []
        parse_state = 'NOTFOUND'
        sequence_field_length = None
        
        for line in table_file:
            line = line.strip()
            if parse_state == 'NOTFOUND':
                match = re.match('(Sequence.*) Domain.*seq-f.*seq-t.*hmm-f.*hmm-t.*score.*E-value', line)
                if match is not None:
                    # Found domain hit header, but still have to skip extra line containing 'underline' characters
                    parse_state = 'HEADERFOUND'
                    sequence_field_length = len(match.group(1))
                continue
            elif parse_state == 'HEADERFOUND':
                # Skip header 'underlines' and start parsing data from next line
                parse_state = 'FOUND'
                continue
            elif parse_state == 'FOUND':
                if len(line) == 0:
                    # No more domain hit data, stop
                    break
                if line!="[no hits above thresholds]":
                    domain_hits.append(DomainHit.parse_line(line, sequence_field_length))
        
        return domain_hits

    @staticmethod
    def parse_line(line_string, name_field_width=None):
        domain_hit = DomainHit()
        line_string = line_string.strip()
        fields = line_string.split()
        domain_hit.sequence = fields[0]
        fields = fields[1:]
        domain_hit.domain_num = int(fields[0].split('/')[0])
        domain_hit.domain_total = int(fields[0].split('/')[1])
        domain_hit.sequence_from = int(fields[1])
        domain_hit.sequence_to = int(fields[2])
        domain_hit.sequence_from_start = fields[3][0] == '['
        domain_hit.sequence_to_end = fields[3][1] == ']'
        domain_hit.hmm_from = int(fields[4])
        domain_hit.hmm_to = int(fields[5])
        domain_hit.hmm_from_start = fields[6][0] == '['
        domain_hit.hmm_to_end = fields[6][1] == ']'
        domain_hit.score = float(fields[7])
        domain_hit.evalue = float(fields[8])
        return domain_hit

    def __str__(self):
        return "%s %d %d %d %d %g" % (self.sequence, self.sequence_from, self.sequence_to, self.hmm_from, self.hmm_to, self.evalue)

    @staticmethod
    def sort_list_on(domainhit_list, fieldname, reverse=False, in_place=False):
        if not inplace:
            return sorted(domainhit_list, key=attrgetter(fieldname), reverse=reverse)
        else:
            domainhit_list.sort(key=attrgetter(fieldname), reverse=reverse)
            return domainhit_list


def main():
    print "Whole-sequence hits:"
    sequence_hits = SequenceHit.parse(sys.argv[1])

    for sequence_hit in sequence_hits:
        print sequence_hit

    print
    print "Domain hits:"
    domain_hits = DomainHit.parse(sys.argv[1])

    for domain_hit in domain_hits:
        print domain_hit


if __name__ == "__main__":
    main()
