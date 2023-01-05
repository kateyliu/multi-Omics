#!/usr/bin/env python
"""
neoantigen_addSample.py- Given a neoantigen filtered.condensed.ranked.tsv 
file, will return a file with an added first col, Sample that is set to the
run name
"""

import os
import sys
from optparse import OptionParser

            
def main():
    usage = "USAGE: %prog -f [file.tsv] -n [run_name] -o [output file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--infile", help="file to filter")
    optparser.add_option("-n", "--name", help="run name")
    optparser.add_option("-o", "--out", help="output file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.infile or not options.name or not options.out:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.infile)
    out = open(options.out, "w")
    hdr = f.readline().strip().split("\t")
    hdr.insert(0,"Sample")
    out.write("%s\n" % "\t".join(hdr))

    #rest of the file
    for l in f:
        tmp = l.strip().split("\t")
        tmp.insert(0, options.name)
        out.write("%s\n" % "\t".join(tmp))
    f.close()
    out.close()

if __name__=='__main__':
    main()
