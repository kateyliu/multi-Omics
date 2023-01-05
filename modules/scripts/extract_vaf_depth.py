#!/usr/bin/env python
"""Script to extract VAF and DEPTH information from a vcf file
BASED on Jingxin Fu's script extract_vaf_depth.py

This script will: 
1. Extract each row's tumor information (col 10, where first col is column 1)
2. parse out the attributes, separated by :
3. take the first and second attribute and calculate VAF and DEPTH from those
   a. VAF = attrib2/sum(attrib1, attrib2)
   b. DEPTH = sum(attrib1, attrib2)

OUTPUT: two column file, VAF and DEPTH to std out
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -v [vcf file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-v", "--vcf", help="vcf file to extract VAF and DEPTH information", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.vcf:
        optparser.print_help()
        sys.exit(-1)

    f=open(options.vcf)
    print("VAF\tDEPTH") #pring new header
    for l in f:
        if l.startswith("#"): #skip vcf header
            continue
        tmp = l.strip().split('\t')
        attribs = tmp[9].split(":") #take tumor column
        x = attribs[1].split(",") #take the 2nd attrib and split the parts
        depth = int(x[0]) + int(x[1])
        if depth > 0:
            vaf = float(int(x[1])/depth)
            print("%s\t%s" % (vaf, depth))
    f.close()

if __name__=='__main__':
    main()
