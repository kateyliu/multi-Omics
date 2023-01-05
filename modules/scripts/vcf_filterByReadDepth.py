#!/usr/bin/env python
"""
vcf_filterByDepth.py- Given a VCF file (.vcf or vcf.gz), and an READ DEPTH/COVERAGE threshold (int),
will return a file with READ DEPTH (DP > threshold)

NOTE: key's in on the DP field int eh 8th and 9th cols
"""

import os
import sys
from optparse import OptionParser
import gzip

def readDepthFilter(infile, threshold, field, outfile):
    """infile - .vcf or .vcf.gz"""
    gzfile = infile.endswith(".gz")
    
    if gzfile:
        ifile = gzip.open(infile, "rb")
    else:
        ifile = open(infile)
    
    ofile = open(outfile, "w")

    for l in ifile:
        if gzfile:
            #must convert from binary to string
            l = l.decode("utf-8")
            
        if l.startswith("#"): 
            ofile.write(l)
            continue

        tmp = l.strip().split("\t") 
        #NOTE: the 8th fields are the "keys" and the 9th are the values
        tags = tmp[8].split(":") #get list of keys
        vals = tmp[9].split(":")
        #generate a dictionary of key:val paris
        d = dict(zip(tags, vals))

        #filter
        if (field in d and int(d[field]) >= int(threshold)):
            #Write to out
            ofile.write("%s\n" % "\t".join(tmp))

def main():
    usage = "USAGE: %prog -v [vcf_file.vcf] -t [threshold (double)] -o [output file.vcf]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-v", "--vcf", help=".vcf/.vcf.gz file to filter")
    optparser.add_option("-t", "--threshold", help="read depth/coverage threshold", default=50)
    optparser.add_option("-f", "--field", help="the read depth field e.g. 'DP' (default) or 'AFDP'", default='DP')
    optparser.add_option("-o", "--out", help="output file .vcf")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.vcf or not options.threshold or not options.out:
        optparser.print_help()
        sys.exit(-1)

    readDepthFilter(options.vcf, options.threshold, options.field, options.out)

if __name__=='__main__':
    main()

