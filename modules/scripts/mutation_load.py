#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to calculate mutation load / MB
OUTPUT: txt file
"""

import os
import sys
from optparse import OptionParser

def Calculate_Nmutations(infile, outfile = "outfile.txt", size=30000000):
    infile_name = infile.lower()
    ofile = open(outfile, "w")
    i=0
    with open(infile,'r') as fin:
        for line in fin:
            if(line[0]!='#') or (line[0]!='Hugo_Symbol'):
                i=i+1
    rate=float(i)/float(size)
    ofile.write("the mutation in protein coding\n" + "the mutation load/MB : " + str(rate))
    ofile.close()


def main():
    usage = "USAGE: %prog -v [maf_file.maf] -o [output file.txt] -s [size]"
    optparser = OptionParser(usage=usage)                    
    optparser.add_option("-v", "--maf", help="maf file to filter")
    optparser.add_option("-o", "--out", help="output file")
    optparser.add_option("-s", "--size", help="size",default=3000000)

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.maf or not options.out:
        optparser.print_help()
        sys.exit(-1)
        
    Calculate_Nmutations(options.maf, options.out, options.size)

if __name__=='__main__':
    main()
