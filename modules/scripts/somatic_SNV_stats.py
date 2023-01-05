#!/usr/bin/env python
"""
generate somatic mutation statistics. output .csv file
"""

import os
import sys
from string import Template
from optparse import OptionParser
            
def main():
    usage = "USAGE: %prog -v [filter.vep.vcf file] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-v", "--vcf", help="filter.vep.vcf file")
    optparser.add_option("-o", "--out", help="output .csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.vcf or not options.out:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.vcf)
    #HACK: get run name from file
    run_name = options.vcf.split("/")[-1].split("_")[0]
    #print(run_name)
    run_ct = {"A>C": 0, "A>G": 0, "A>T": 0,
              "C>A": 0, "C>G": 0, "C>T": 0,
              "G>A": 0, "G>C": 0, "G>T": 0,
              "T>A": 0, "T>G": 0, "T>C": 0,
              "A>A": 0, "C>C":0, "G>G":0, "T>T":0, #diagonal should be 0s
    }
    for l in f:
        if l.startswith("#"): #skip headers
            continue
        tmp = l.strip().split("\t")
        #print(tmp)
        mutation = "%s>%s" % (tmp[3],tmp[4])
        if (mutation in run_ct):
            run_ct[mutation] +=1
        #else: #NOT a snp
            #print(mutation)
    f.close()
    #print(run_ct)

    out = open(options.out, "w")
    col_sums=[0,0,0,0,0]
    out.write("%s\n" % ",".join(["Ref/Alt","A", "C","G","T","Total"]))
    for r in ["A","C","G","T"]:
        cts = [run_ct["%s>%s" % (r,c)] for c in ["A","C","G","T"]]
        total = sum(cts)
        cts.append(total)
        #Add cts to each col's running total in col_sums
        for (i, val) in enumerate(cts):
            col_sums[i] += val
        tmp = [str(i) for i in cts]
        tmp.insert(0, r) #add row name
        out.write("%s\n" % ",".join(tmp))
    out.write("Total,%s\n" % ",".join([str(s) for s in col_sums]))
    out.close()

if __name__=='__main__':
    main()
