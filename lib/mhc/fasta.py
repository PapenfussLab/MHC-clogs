#!/usr/bin/env python
"""
A collection of useful FASTA file functions, including parsing, splitting and
6-frame translation.  

When invoked as main, this script can perform some useful operations on FASTA
files.  See the option help output for details (-h)

TODO: Verify that the six frame translation is valid (especially remember case switching).
TODO: Verify that the six frame coordinate 'untranslation' is valid.

Author: Chris Davoren
"""

import os
import re
import sys
import optparse
import sixframe

split_keyword = 's'
postsplit_keyword = 'ps'
translation_keyword = 't'


def ensembl_match(string):
    """
    Returns decoded information from an Ensembl FASTA sequence header.  See the
    notes on Ensembl FASTA format headers here:
    http://asia.ensembl.org/info/data/ftp/index.html 

    The identifier, data type, coordinate type, version, name, start, end and
    strand coordinates are returned.  For example, given the header
    'HSCHR5_1_CTG5 dna:chromosome
    chromosome:GRCh37:HSCHR5_1_CTG5:161800673:161974131:1', this method would
    return:

    ('HSCHR5_1_CTG5', 'dna', 'chromosome', 'GRCh37', 'HSCHR5_1_CTG5',
    161800673, 161974131, 1)

    If the 'dna:chromosome' field were absent, the second tuple element would
    be None.

    @param string The FASTA header to decode
    @return A tuple containing the decoded field data
    """
    ensembl_regex = "^(\\S*).*?(\\w+):([^:]*):([^:]+):([0-9]*):([0-9]*):([0-9])"
    ensembl_regex_full = "^(\\S*).*?(\\w+):(\\w+).*?(\\w+):([^:]*):([^:]+):([0-9]*):([0-9]*):([0-9])"

    match_full = re.match(ensembl_regex_full, string)
    if match_full is not None:
        (identifier, data_type, coord_type, version, name, start, end, strand) = (
            match_full.group(1),
            match_full.group(2),
            match_full.group(4),
            match_full.group(5),
            match_full.group(6),
            int(match_full.group(7)),
            int(match_full.group(8)),
            int(match_full.group(9)))
        return (identifier, data_type, coord_type, version, name, start, end, strand)

    match = re.match(ensembl_regex, string)
    if match is not None:
        (identifier, coord_type, version, name, start, end, strand) = (
            match.group(1),
            match.group(2),
            match.group(3),
            match.group(4),
            int(match.group(5)),
            int(match.group(6)),
            int(match.group(7)))
        return (identifier, None, coord_type, version, name, start, end, strand)

    return None


def translation_match(string):
    """
    Returns the 6-frame translation offset decoded from a FASTA sequence header.

    @param string The FASTA header to decode
    @return The 6-frame translation offset, or None if no valid translation
            field is found.  The offset is an integer value from the set
            (+1,+2,+3,-1,-2,-3) 
    """
    translation_regex = "%s:([+-][1-3]):([0-9]+)" % translation_keyword
    match = re.search(translation_regex, string)
    if match is not None:
        return (int(match.group(1)), int(match.group(2))) 
    translation_regex = "%s:([+-][1-3])" % translation_keyword
    match = re.search(translation_regex, string)
    if match is not None:
        return (int(match.group(1)), None)
    translation_regex = "translation frame ([+-][1-3])"
    match = re.search(translation_regex, string)
    if match is not None:
        return (int(match.group(1)), None)
    return None


def split_match(string):
    """
    Returns the PRE 6-frame translation split range decoded from a FASTA
    sequence header.

    @param string the FASTA header to decode
    @return A 3-tuple containing the start and end points of the split and the original length of the sequence.  These
            coordinates are usually 1-based and non-inclusive.
    """
    split_regex = '%s:([0-9]+):([0-9]+):([0-9]+)' % split_keyword
    match = re.search(split_regex, string)
    if match is None:
        return None
    return (int(match.group(1)), int(match.group(2)), int(match.group(3)))


def sixframesplit_match(string):
    """
    Returns the post 6-frame translation split range decoded from a FASTA
    sequence header.

    @param string The FASTA header to decode
    @return A 2-tuple containing the start and end points of the split.  These
            coordinates are usually 1-based and non-inclusive.
    """
    sixframesplit_regex = "%s:([0-9]+):([0-9]+):([0-9]+)" % postsplit_keyword
    match = re.search(sixframesplit_regex, string)
    if match is None:
        return None
    return (int(match.group(1)), int(match.group(2)), int(match.group(3)))


def untranslate_coords(sequence_start, sequence_end, real_length, translation):
    """
    Takes start and end point on a 6-frame translated sequence and returns the
    corresponding coordinates on the original, untranslated base sequence.

    The coordinates given and returned are expected to be 1-based and
    non-inclusive.

    @param sequence_start The start point on the translated sequence
    @param sequence_end The end point on the translated sequence
                        (non-inclusive)
    @param real_length The length of the original, untranslated sequence
    @param translation The offset of the 6-frame translated sequence; expected
                       to be an integer in the set (+1,+2,+3,-1,-2,-3) 
    @return A 2-tuple containing the corresponding start and end coordinates on
            the untranslated sequence
    """
    if translation > 0:
        real_sequence_start = ((sequence_start-1)*3)+(translation-1)+1
        real_sequence_end = ((sequence_end-1)*3)+(translation-1)+3
    else:
        translation *= -1
        # Order swapped because these are from the -ve translations
        real_sequence_end = (real_length-((sequence_start-1)*3))-(translation-1)
        real_sequence_start = ((real_length-(sequence_end*3))-(translation-1))+1
    return (real_sequence_start, real_sequence_end)


def unsplit_coords(sequence_start, sequence_end, split_start):
    """
    Takes start and end points on a split segment of a sequence and returns the
    corresponding coordinates on the original, non-split sequence.

    Coordinates given and returned are expected to be 1-based and
    non-inclusive.

    @param sequence_start The start point on the split sequence
    @param sequence_end The end point on the split sequence (non-inclusve)
    @param split_start The point on the original, non-split sequence where this
           split segment starts (1-based)
    @return A 2-tuple containing the corresponding start and points on the
            non-split sequence
    """
    real_sequence_start = (sequence_start + split_start) - 1
    real_sequence_end = (sequence_end + split_start) - 1
    return (real_sequence_start, real_sequence_end)


def ensembl_str(ensembl_data):
    """
    Reconstructs Ensembl FASTA header fields from a tuple containing decoded
    data.  This is the reverse of the function ensembl_match() except that it
    does not include the original identifier (in other words, it includes the
    extra fields only).

    @param ensembl_data A tuple containing decoded field data
    @returns A string representing the given data in the original Ensembl
             format
    """
    if ensembl_data[1] is not None:
        return "%s:%s %s:%s:%s:%d:%d:%d" % (
            ensembl_data[1],
            ensembl_data[2],
            ensembl_data[2],
            ensembl_data[3],
            ensembl_data[4],
            ensembl_data[5],
            ensembl_data[6],
            ensembl_data[7])
    else:
        return "%s:%s:%s:%d:%d:%d" % (
            ensembl_data[2],
            ensembl_data[3],
            ensembl_data[4],
            ensembl_data[5],
            ensembl_data[6],
            ensembl_data[7])


class Fasta:
    def __init__(self):
        self.sequences = {}
        self.ensembl = {}

    def __str__(self):
        return "FASTA: %d sequences" % len(self.sequences)

    def print_summary(self, dest=sys.stdout, sort=False):
        print >> dest, "Contains %d sequences:" % len(self.sequences)
        sequence_names = self.sequences.keys()
        if sort:
            sequence_names.sort()
        for sequence_name in sequence_names:
            sequence = self.sequences[sequence_name]
            print >> dest, "\t%s\t%d" % (sequence_name, len(sequence)),

            if sequence_name in self.ensembl:
                print >> dest, "\t%s" % ensembl_str(self.ensembl[sequence_name])
            else:
                print >> dest

    def split(self, size=10000000, header_format=split_keyword+":%d:%d:%d", append_to_id=False):
        """
        Returns a new Fasta object where each sequence has been split into
        smaller sequences of equal size.  A field is added to each sequence
        header that notes the start and end points of the split.

        Note that split header information is expected to be 1-based and
        non-inclusive.

        @param size The size of each subsequence (for a given sequence, the
                    last subsequence may be smaller than this)
        @param header_format If specified, the format string of the header to
                             add.  The default is 'split:%d:%d'
        @param append_to_id If true, the segment info is appended to the first
                            token in the sequence header instead of concatenated
                            with the existing whole sequence header (added for 
                            HMMER2 compatibility)
        @return A new Fasta object containing the smaller, split sequences
        """
        split_fasta = Fasta()
        for sequence_name, sequence in self.sequences.items():
            pos = 0
            while pos < len(sequence):
                # Last segment is a special case
                if pos+size > len(sequence):
                    segment_sequence = sequence[pos:]
                    segment_end = len(sequence)
                else:
                    segment_sequence = sequence[pos:pos+size]
                    segment_end = pos+size
                header_addition = header_format % (pos+1, segment_end, len(sequence))
                if not append_to_id:
                    split_sequence_name = "%s %s" % (sequence_name, header_addition)
                else:
                    sequence_id = sequence_name.split()[0]
                    sequence_rest = sequence_name[len(sequence_id)+1:]
                    split_sequence_name = "%s|%s %s" % (sequence_id, header_addition, sequence_rest)
                if sequence_name in self.ensembl:
                    split_sequence_name += " " + ensembl_str(self.ensembl[sequence_name])
                split_fasta.sequences[split_sequence_name] = segment_sequence
                pos += size
        return split_fasta

    def translate(self, header_format=translation_keyword+":%s:%d", append_to_id=False):
        """
        Returns a new Fasta object representing the 6-frame translation of this
        Fasta.  Note that this assumes that this Fasta contains genomic (rather
        than protein) data.

        A field is added to each sequence header that contains the translation
        offset; the value of this field will be one of the set
        (+1,+2,+3,-1,-2,-3)

        @param header_format If specified, the custom format string of the
                             header to add.  The default is 'translation:%s'
        @param append_to_id If true, the translation header is appended to the
                            first token in the sequence header instead of 
                            concatenated with the existing whole sequence header
                            (added for HMMER2 compatibility)
        @return A new, 6-frame translated Fasta object
        """
        def create_header(sequence_name, sequence_length, header_format, translation, append_to_id):
            if not append_to_id:
                return "%s %s" % (sequence_name, header_format % (translation, sequence_len))
            else:
                sequence_id = sequence_name.split()[0]
                sequence_rest = sequence_name[len(sequence_id)+1:]
                return "%s|%s %s" % (sequence_id, header_format % (translation, sequence_length), sequence_rest)

        translated_fasta = Fasta()
        for sequence_name, sequence in self.sequences.items():
            # print >> sys.stderr, "Translating: %s" % sequence_name
            # Switch to upper case to match codon table
            sequence = sequence.upper()
            header = create_header(sequence_name, len(sequence), header_format, "+1", append_to_id)
            translated_fasta.sequences[header] = sixframe.frame_translate(sequence, 0, True)
            header = create_header(sequence_name, len(sequence), header_format, "+2", append_to_id)
            translated_fasta.sequences[header] = sixframe.frame_translate(sequence, 1, True)
            header = create_header(sequence_name, len(sequence), header_format, "+3", append_to_id)
            translated_fasta.sequences[header] = sixframe.frame_translate(sequence, 2, True)
            sequence = sixframe.reverse_complement(sequence)
            header = create_header(sequence_name, len(sequence), header_format, "-1", append_to_id)
            translated_fasta.sequences[header] = sixframe.frame_translate(sequence, 0, True)
            header = create_header(sequence_name, len(sequence), header_format, "-2", append_to_id)
            translated_fasta.sequences[header] = sixframe.frame_translate(sequence, 1, True)
            header = create_header(sequence_name, len(sequence), header_format, "-3", append_to_id)
            translated_fasta.sequences[header] = sixframe.frame_translate(sequence, 2, True)
        return translated_fasta

    def write(self, line_length=50, dest=sys.stdout, sort=False):
        """
        Prints this Fasta out in standard FASTA format.

        @param line_length The length of each line of data
        @param dest The file to print to; the default is standard output
        """
        sequence_names = self.sequences.keys()
        if sort:
            sequence_names.sort()
        for sequence_name in sequence_names:
            sequence = self.sequences[sequence_name]
            print >> dest, ">%s" % sequence_name
            pos = 0
            while pos < len(sequence):
                print >> dest, sequence[pos:pos+line_length]
                pos += line_length;


def fasta_load(filename=None, file_handle=None, match=None, ensembl_detect=False):
    if filename is not None:
        file_handle = open(filename, 'r')
    if file_handle is not None:
        fasta = fasta_extract_sequences(file_handle, match=match, ensembl_detect=ensembl_detect)
        if filename is not None:
            file_handle.close()
        return fasta
    else:
        return None

    
def fasta_extract_sequences_print(file_handle, match=None):
    current_sequence_name = None
    current_sequence_print = False
    counter = 0
    for line in file_handle:
        line = line.strip()
        if not len(line):
            continue
        if line[0] == '>':
            if current_sequence_name is not None and current_sequence_print:
                print
            current_sequence_name = line[1:]
            if match is None:
                current_sequence_print = True
            elif isinstance(match, str):
                current_sequence_print = re.match(match, current_sequence_name) is not None
            elif isinstance(match, re.RegexObject):
                current_sequence_print = match.match(current_sequence_name) is not None
            else:
                current_sequence_print = False

            if current_sequence_print:
                counter += 1
                print >> sys.stderr, "Matched %d sequences..." % counter
                print line
            else:
                print >> sys.stderr, "Discarded sequence %s..." % current_sequence_name
                pass
        elif current_sequence_print:
            print line


def fasta_extract_sequences(file_handle, match=None, ensembl_detect=False):
    current_sequence_name = None
    current_sequence_match = False
    current_sequence = ''
    counter = 0
    fasta = Fasta()
    for line in file_handle:
        line = line.strip()
        if not len(line):
            continue
        if line[0] == '>':
            if current_sequence_name is not None and current_sequence_match:
                fasta.sequences[current_sequence_name] = current_sequence
                current_sequence = ''
            current_sequence_name = line[1:]
            if match is None:
                current_sequence_match = True
            elif isinstance(match, str):
                current_sequence_match = re.match(match, current_sequence_name) is not None
            elif isinstance(match, re.RegexObject):
                current_sequence_match = match.match(current_sequence_name) is not None
            else:
                current_sequence_match = False

            if current_sequence_match:
                # print >> sys.stderr, "Matched sequence (%d total): %s" % (len(fasta.sequences)+1, current_sequence_name)
                counter = 0
                if ensembl_detect:
                    ematch = ensembl_match(current_sequence_name)
                    if ematch is not None:
                        fasta.ensembl[current_sequence_name] = ematch
            else:
                print >> sys.stderr, "Discarded sequence: %s" % current_sequence_name

        elif current_sequence_match:
            """
            counter += len(line)
            if counter % 1000000 == 0:
                print "\t%d" % counter
            """
            current_sequence += line
    # Store last sequence
    if current_sequence_name is not None and current_sequence_match:
        fasta.sequences[current_sequence_name] = current_sequence
    return fasta
            

def main():
    parser = optparse.OptionParser('%prog <input_fasta.fa>')
    parser.add_option("-m", "--match", help="Retrieve only sequences that match the given regular expression", dest="match", default=None, metavar="REGEX")
    parser.add_option("-e", "--ensembl", help="Detect and preserve ensembl header attributes", dest="ensembl", action="store_true", default=False)
    parser.add_option("-s", "--split", help="Split all sequences into segments of length SIZE; appends a 'split' field to each split sequence header", metavar="SIZE", dest="split", default=None)
    parser.add_option("-q", "--quick", help="Echo data quickly without storing it first", dest="quick", action="store_true", default=False)
    parser.add_option("-z", "--summarize", help="Print summary instead of outputting full FASTA data", dest="summarize", action="store_true", default=False)
    parser.add_option("-l", "--linelength", help="If -f is specified, this specifies the line length to use when printing FASTA output (default=50", dest="linelength", metavar="LENGTH", default=50)
    parser.add_option("-t", "--translate", help="Perform a 6-frame translation of the data; appends a 'translate' field to each translated sequence header", dest="translate", action="store_true", default=False)
    parser.add_option("-p", "--postsplit", help="Like --split, but the split is performed after the sixframe translation; adds a 'sixframesplit' field to each sequence header", dest="postsplit", default=None, metavar="SIZE")
    parser.add_option("-o", "--sort", help="Sort sequences by header before output", dest="sort", action="store_true", default=False)

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    if options.quick:
        fasta_file = open(args[0], 'r')
        fasta_extract_sequences_print(fasta_file, options.match)
        fasta_file.close()
        sys.exit(0)

    fasta_file = open(args[0], 'r')
    fasta = fasta_extract_sequences(fasta_file, options.match, options.ensembl)
    fasta_file.close()

    # print >> sys.stderr, "Matched sequences (%d total):" % len(fasta.sequences)
    for sequence_name in fasta.sequences:
        print >> sys.stderr, "\t%s" % sequence_name

    if options.split is not None:
        options.split = int(options.split)
        fasta = fasta.split(size=options.split, append_to_id=True)

    if options.translate:
        fasta = fasta.translate(append_to_id=True)

    if options.postsplit is not None:
        options.postsplit = int(options.postsplit)
        fasta = fasta.split(size=options.postsplit, header_format=postsplit_keyword+":%d:%d:%d", append_to_id=True)

    if options.summarize:
        fasta.print_summary(sort=options.sort)
    else:
        options.linelength = int(options.linelength)
        fasta.write(line_length=options.linelength, sort=options.sort)

    sys.exit(0)


if __name__ == "__main__":
    main()
