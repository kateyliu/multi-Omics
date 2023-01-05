#!/usr/bin/env python
"""
Len Taing 2020 (TGBTG)
generate somatic mutation statistics. output .csv file
"""

import os
import sys
from string import Template
from optparse import OptionParser

def countBedBases(ffile):
    """Count the number of base pairs in a bed file"""
    count = 0
    f = open(ffile)
    for l in f:
        tmp = l.split("\t")
        (start, end) = (int(tmp[1]), int(tmp[2]))
        count = count + end - start
    f.close()
    return count

def getVariantStats(maf_f, run_name, maf_isString=False):
    """Calculate summary statistics for SNP, INS, DEL, also 
    Functional annotations
    
    Usuaully take in a maf file, but can also accept maf as string--
    see 3rd parm
    """
    cts = {}
    annot_cts = {'Missense_Mutation':{"SNP": 0, "INS": 0, "DEL": 0},
                 'Nonsense_Mutation':{"SNP": 0, "INS": 0, "DEL": 0},
                 'Silent':{"SNP": 0, "INS": 0, "DEL": 0},
    }

    if not maf_isString:
        f = open(maf_f)
    else:
        f = maf_f
        
    #print(run_name)
    run_ct = {"SNP": 0, "INS": 0, "DEL": 0}
    for l in f:
        if l.startswith("#"):
            continue
        tmp = l.strip().split("\t")
        #variant classification is 9th col, e.g. Missense_Mutation, Nonsense_Mutation
        #variant type is 10th col, e.g. SNP/INS/DEL
        classification = tmp[8]
        label = tmp[9]
        if label in run_ct:
            #print(label)
            run_ct[label] += 1
            #HANDLE annotation- i.e. count Missense,Nonsense,Silent
            if classification in annot_cts:
                annot_cts[classification][label] +=1
    if not maf_isString:
        f.close()
    run_ct['TOTAL']= sum([run_ct['SNP'], run_ct['INS'], run_ct['DEL']])
    for k in annot_cts.keys():
        annot_cts[k]['TOTAL']= sum([v for v in annot_cts[k].values()])
        
    cts[run_name] = run_ct

    return (cts, annot_cts)

def main():
    usage = "USAGE: %prog -m [filter.maf files list] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--maf", help="filter.maf file")
    optparser.add_option("-t", "--target", help="target region bed file")
    optparser.add_option("-o", "--out", help="output mutation variant type count .csv file")
    optparser.add_option("-a", "--annot", help="output functional annotations count .csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.maf or not options.target or not options.out or not options.annot:
        optparser.print_help()
        sys.exit(-1)

    #Calculate the number of target_Mbs
    target_mb = round(countBedBases(options.target)/float(10**6), 1)
    #HACK: get run name from file
    run_name = options.maf.split("/")[-1].split("_")[0]
    (cts, annot_cts) = getVariantStats(options.maf, run_name)

    out = open(options.out, "w")
    out.write("%s\n" % ",".join(["Run","TOTAL", "SNP","INSERT","DELETE","TMB (mut/Mb)"]))
    for r in sorted(cts.keys()):
        tmp = cts[r]
        #NOTE: TMB = total tumor variants (tumor_tmb)/number of target Mb
        TMB = round(tmp['TOTAL']/target_mb, 1)
        
        out.write("%s\n" % ",".join([r, str(tmp['TOTAL']), str(tmp["SNP"]),
                                     str(tmp["INS"]), str(tmp["DEL"]), str(TMB)]))
    out.close()

    #WRITE out annotation counts
    annot_out = open(options.annot, "w")
    annot_out.write("%s\n" % ",".join(["Classification","TOTAL", "SNP","INSERT","DELETE"]))
    for (classification, counts) in annot_cts.items():
        annot_out.write("%s\n" % ",".join([classification,
                                             str(counts['TOTAL']),
                                             str(counts["SNP"]),
                                             str(counts["INS"]),
                                             str(counts["DEL"])]))
    annot_out.close()

if __name__=='__main__':
    main()
