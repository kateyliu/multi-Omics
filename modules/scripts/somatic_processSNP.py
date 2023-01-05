#!/usr/bin/env python
"""
generate somatic mutation statistics. output .csv file
"""

import os
import sys
from string import Template
from optparse import OptionParser
            
def main():
    usage = "USAGE: %prog -m [filter.maf file] (-w window size: default 500000)"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--maf", help="filter.maf files")
    optparser.add_option("-w", "--window", help="window size (500000)", default=500000)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.maf:
        optparser.print_help()
        sys.exit(-1)

    window = int(options.window)
    chromLens = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
                 'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
                 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
                 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
                 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                 'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
                 'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
                 'chr22': 50818468 , 'chrX': 156040895, 'chrY': 57227415}
    order = ['chr1', 'chr2', 'chr3',  'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
             'chr9','chr10','chr11', 'chr12','chr13','chr14','chr15','chr16',
             'chr17','chr18', 'chr19',  'chr20','chr21','chr22','chrX','chrY']
    
    #STORE them as a dictionary (chroms) lists of lists
    #initialize the dictionary of array of 0s
    indel = {}
    for chrom in chromLens:
        n = int(chromLens[chrom]/window) + 1
        indel[chrom] = [0 for i in range(n)]
        
    f = open(options.maf)
    for l in f:
        if l.startswith("#"):
            continue
        tmp = l.strip().split("\t")
        #Classification is 10th col
        chrom = tmp[4]
        start = tmp[5]
        end = tmp[6]
        label = tmp[9]
        
        #calculate which bin it should belong
        if chrom in chromLens and label == "SNP":
            start_bin = int(int(start) / window)
            indel[chrom][start_bin] += 1
    f.close()
    
    #convert this into a circos output
    #hdr
    #print("\t".join([hdr[0],hdr[1],hdr[2],hdr[4]]))
    for chrom in order:
        for i in range(int(chromLens[chrom]/window) + 1):
            #convert to hs
            tmp = "hs%s" % chrom[3:]
            start = i*window
            end = start + window
            if (end > chromLens[chrom]):
                 end =chromLens[chrom]
            
            print("\t".join([tmp, str(start), str(end), "%s" % indel[chrom][i]]))
    
if __name__=='__main__':
    main()
