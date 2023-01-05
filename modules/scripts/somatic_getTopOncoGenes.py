#!/usr/bin/env python
"""
Len Taing 2021 (TGBTG)
Given a MAF and a list of top onco driving genes (as a TSV file with one gene
per row), generate an output (TSV--one gene per row) of the genes represented
in the maf (SORTED by gene name)
"""

import os
import sys
from string import Template
from optparse import OptionParser


def getOncoGeneList(maf_f, oncoGeneList, maf_isString=False):
    """Given a maf file (as a string or as a filepath), and a list of top onco
    driving genes, returns a list of the genes represented in the maf
    """

    genes_in_maf = []
    if not maf_isString:
        f = open(maf_f)
    else:
        f = maf_f
        
    for l in f:
        if l.startswith("#") or l.startswith('Hugo_Symbol'):
            continue
        
        tmp = l.strip().split("\t")
        #GENE name is first elm
        gene = tmp[0]
        if gene in oncoGeneList and gene not in genes_in_maf:
            genes_in_maf.append(gene)
    if not maf_isString:
        f.close()

    #print(sorted(genes_in_maf))
    return sorted(genes_in_maf)

def main():
    usage = "USAGE: %prog -m [filter.maf files list] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--maf", help="filter.maf file")
    optparser.add_option("-l", "--oncoGeneList", help="list of top onco driving genes")
    optparser.add_option("-o", "--out", help="output mutation variant type count .csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.maf or not options.oncoGeneList or not options.out:
        optparser.print_help()
        sys.exit(-1)

    #READ in the cancerGeneList
    oncoGeneList = []
    f = open(options.oncoGeneList)
    for l in f:
        oncoGeneList.append(l.strip())
    f.close()

    genes = getOncoGeneList(options.maf, oncoGeneList)

    #WRITE to output
    out = open(options.out, "w")
    for g in genes:
        out.write("%s\n" % g)
    out.close()

if __name__=='__main__':
    main()
