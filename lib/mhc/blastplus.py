"""mhc.blastplus"""

import os


def get_sequence(blast_bb, accession, start=0, end=0, strand='+', padding=0, debug=False):
    """Load a sequence from a BLAST database.
    
    @param blast_bb: BLAST database
    @param accession: Accession name
    @param start: Start coordinate (Default: 0, extract from start of sequence)
    @param end: End coordinate (Default: 0, extract to the end of sequence)
    @param strand: Strand +/- (Default: +)
    @param padding: Sequence padding (Default: 0)
    @returns: (header,seq)
    """
    
    strand = {
        "plus": "plus", "minus": "minus", 
        "+": "plus", "-": "minus",
        "+1": "plus", "-1": "minus"
    }[strand]
        
    if start>end: start,end = end,start
    
    cmd = 'blastdbcmd -db %s -entry "%s" -range %i-%i -strand %s -outfmt %%f' \
        % (blast_bb, accession, start, end, strand)
    
    p = os.popen(cmd)
    header = p.readline()[1:].strip()
    if not header:
        raise Exception('BLAST failure')
    
    seq = []
    for line in p:
        seq.append(line.strip())
    seq = ''.join(seq)
    
    if not seq:
        print blast_db, accession
        raise NotFoundException()
    
    return header,seq
