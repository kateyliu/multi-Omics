#!/usr/bin/env python
"""Script to collect the sample summaries of each sample's _target_metrics.txt.sample_summary file

OUTPUTS to stdout:
HEADER: 
sample_id	total	mean	granular_Q1	granular_median	granular_Q3   %_bases_above_50

AND then the second row for each of the files given
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FPKM FILE_1] -f [FPKM FILE_2] ...-f [FPKM FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--cov_files", action="append", help="list of _target_metrics.txt.sample_summary files")
    optparser.add_option("-a", "--align_files", action="append", help="list of _mapping.txt files")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.cov_files or not options.align_files:
        optparser.print_help()
        sys.exit(-1)
    #check that they're the same # of files in each list
    if len(options.cov_files) != len(options.align_files):
        print("The number of coverage files %s, do not match the number of align_files %s" % (len(options.cov_files),len(options.align_files)))
        sys.exit(-1)

    hdr = "\t".join(["sample_id","total","mean","granular_Q1","granular_median","granular_Q3","%_bases_above_50"])
    table = []
    for (i, ff) in enumerate(options.cov_files):
        f = open(ff)
        #WE only want the info from the second line
        tmp = f.readline() #read in first line aka header
        tmp = f.readline().strip().split("\t")
        f.close()
        table.append(tmp)
    #BUT now we have to replace the totals b/c those are TOTAL DEPTH, not
    #TOTAL reads, which is from that mapping.txt files

    for (i, ff) in enumerate(options.align_files):
        f = open(ff)
        total = f.readline().strip().split()[0]
        f.close()

        #REPLACE total value in the table
        table[i][1] = total

    #PRINT output
    print(hdr)
    for row in table:
        print("\t".join(row))
        
if __name__=='__main__':
    main()


