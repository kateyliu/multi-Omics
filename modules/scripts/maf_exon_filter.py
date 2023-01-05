#!/usr/bin/env python
"""
maf_filter.py- Given a maf file, will return a file with coding exon mutations.
NOTE: key's in on the BIOTYPE field in 63th col (63th in python).
"""

import os
import sys
from optparse import OptionParser

def mafExonFilter(infile, outfile = "outfile.maf"):
    ifile = open(infile)
    ofile = open(outfile, "w")

    for l in ifile:
        if l.startswith("#") or l.startswith("Hugo_Symbol"):
            ofile.write(l)
            continue
        tmp = l.strip().split("\t")
        #NOTE the 63th fields can identify exon mutations by "protein_coding"
        if (tmp[63] == "protein_coding"):
            ofile.write(l)
            
def main():
    usage = "USAGE: %prog -m [maf_file.maf] -o [output file.maf]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--maf", help="maf file to filter")
    optparser.add_option("-o", "--out", help="output file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.maf or not options.out:
        optparser.print_help()
        sys.exit(-1)

    mafExonFilter(options.maf, options.out)

if __name__=='__main__':
    main()
