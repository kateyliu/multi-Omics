#!/usr/bin/env python3
"""
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [filtered.tsv] -o [output file.txt]"
    optparser = OptionParser(usage=usage)                    
    optparser.add_option("-f", "--file", help="neoantigen filter.tsv file")
    optparser.add_option("-r", "--rna", help="include rna evidence (default: False)", action="store_true", default=False)
    optparser.add_option("-o", "--out", help="output file")

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.out:
        optparser.print_help()
        sys.exit(-1)
        
    f = open(options.file)
    #Fields we want to pick out from the pvacseq filtered output file
    fields = ['Gene Name', 'Mutation', 'Variant Type', 'HGVSc', 'HGVSp',
              'HLA Allele', 'MT Epitope Seq', 
              'Tumor DNA Depth', 'Tumor DNA VAF',
              'Best MT Score', 'Corresponding WT Score', 'Corresponding Fold Change']
    #ADD  'Tumor RNA Depth', 'Tumor RNA VAF'
    if options.rna:
        fields.insert(9, 'Tumor RNA Depth')
        fields.insert(10, 'Tumor RNA VAF')
    hdr = f.readline().strip().split("\t")
    results = []
    for l in f:
        tmp = dict(zip(hdr, l.strip().split("\t")))
        ls = [tmp[fld] for fld in fields]
        #print(ls)
        results.append(ls)
    f.close()

    out = open(options.out, 'w')
    #print out the header
    #SHORTEN the long fields names:
    _map = {'MT Epitope Seq': 'Epitope Seq',
            'Tumor DNA Depth': 'DNA Depth',
            'Tumor DNA VAF': 'DNA VAF',
            'Tumor RNA Depth': 'RNA Depth',
            'Tumor RNA VAF': 'RNA VAF',
            'Best MT Score': 'MT Score',
            'Corresponding WT Score': 'WT Score',
            'Corresponding Fold Change': 'Fold Change'}
    new_fields = list(map(lambda x: _map[x] if x in _map else x, fields))
    out.write("%s\n" % "\t".join(new_fields))
    for row in results:
        out.write("%s\n" % "\t".join(row))
    out.close()
    
if __name__=='__main__':
    main()
