#!/usr/bin/env python

import os, sys, optparse, re, hmmer2

class Component:
    def __init__(self):
        self.name = None
        self.chromosome = None
        self.domain_count = None
        self.sequence_start = None
        self.sequence_end = None
        self.sequence_code = None
        self.hmm_start = None
        self.hmm_end = None
        self.hmm_code = None
        self.score = None
        self.evalue = None
        self.strand = None

    def __str__(self):
        return "%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%g\t%g\t%s" % (
            self.name,
            self.chromosome,
            self.domain_count,
            self.sequence_start,
            self.sequence_end,
            self.sequence_code,
            self.hmm_start,
            self.hmm_end,
            self.hmm_code,
            self.score,
            self.evalue,
            self.strand)

    @staticmethod
    def parse(line):
        line = line.split("\t")
        component = Component()
        component.name = line[0]
        component.chromosome = line[1]
        component.domain_count = line[2]
        component.sequence_start = int(line[3])
        component.sequence_end = int(line[4])
        component.sequence_code = line[5]
        component.hmm_start = int(line[6])
        component.hmm_end = int(line[7])
        component.hmm_code = line[8]
        component.score = float(line[9])
        component.evalue = float(line[10])
        component.strand = line[11]
        return component

class Chain:
    def __init__(self):
        self.alpha1 = None
        self.alpha2 = None
        self.ig = None
        self.cterm = None
        self.score = None
        self.length = None
        self.start = None
        self.end = None
        self.strand = None
        self.chromosome = None
        self.num = None

    @staticmethod
    def parse(filename):
        chains = []
        mode = 'CHAIN'
        fp = open(filename, 'r')
        chain = None
        for line in fp:
            line = line.strip()
            if len(line) > 0:
                if mode == 'COMPONENT':
                    component = Component.parse(line)
                    # chain[component.name] = component
                    if component.name == 'alpha1':
                        chain.alpha1 = component
                    elif component.name == 'alpha2':
                        chain.alpha2 = component
                    elif component.name == 'Ig':
                        chain.ig = component
                    elif component.name == 'C-term':
                        chain.cterm = component
                elif line[0] == '>':
                    chain = Chain()
                    lm = re.search('Length=([0-9]+)', line)
                    if lm is not None:
                        chain.length = int(lm.group(1))
                    sm = re.search('Score=([0-9.]+)', line)
                    if sm is not None:
                        chain.score = float(sm.group(1))
                    mode = 'COMPONENT'
            else:
                mode = 'CHAIN'
                if chain is not None:
                    chains.append(chain)
                    chain = None
        if chain is not None:
            chains.append(chain)
        for i, chain in enumerate(chains):
            chain.num = i+1
            (chain.start, chain.end) = chain.get_bounds()
            chain.strand = chain.get_strand()
            chain.chromosome = chain.get_chromosome()
        return chains

    def get_chromosome(self):
        l = [self.alpha1, self.alpha2, self.ig, self.cterm]
        for d in l:
            if d is not None:
                return d.chromosome
        return None

    def get_strand(self):
        l = [self.alpha1, self.alpha2, self.ig, self.cterm]
        for d in l:
            if d is not None:
                return d.strand
        return None

    def get_bounds(self):
        bounds = []
        l = [self.alpha1, self.alpha2, self.ig, self.cterm]
        for d in l:
            if d is not None:
                bounds += [d.sequence_start, d.sequence_end]
        return (min(bounds), max(bounds))


def main():
    parser = optparse.OptionParser("%prog <input_chain_file>")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit()

    chains = Chain.parse(args[0])

    print >> sys.stdout, "Total chains read: %d" % len(chains)
    for i, chain in enumerate(chains):
        print >> sys.stdout, "\tchain %d: length=%d, score=%f" % (i+1, chain.length, chain.score)
        print >> sys.stdout, "\t\t%s" % (chain.alpha1)
        print >> sys.stdout, "\t\t%s" % (chain.alpha2)
        print >> sys.stdout, "\t\t%s" % (chain.ig)
        print >> sys.stdout, "\t\t%s" % (chain.cterm)

if __name__ == "__main__":
    main()
