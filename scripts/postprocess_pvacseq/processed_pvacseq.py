#!/usr/bin/env python

import os
import sys
import subprocess
import yaml
from optparse import OptionParser
import pandas as pd
import numpy as np

def main():
   usage = "USAGE: multiomics % python processed_pvacseq.py -f Pt_3986_NL_AS.all_epitopes.tsv -s Pt1 -o Pt1_processed_output.csv -d /Users/aashna/desktop/multiomics/"
   optparser = OptionParser(usage=usage)
   optparser.add_option("-f", "--file", help="input pvacseq file", default=None)
   optparser.add_option("-o", "--output", help="output file", default=None)
   optparser.add_option("-d", "--destination", help="output path", default=None)
   optparser.add_option("-s", "--samplename", help="name of the sample being processed",default=None)
   (options, args) = optparser.parse_args(sys.argv)

   if not options.file or not options.samplename or not options.destination:
       optparser.print_help()
       sys.exit(-1)
   input_file = options.file
   output_file = options.output
   destination = options.destination
   sample_name= options.samplename

   input = pd.read_csv(input_file,sep='\t')
   input_df = pd.DataFrame(input)
   input_df['WT Epitope Seq'] = input_df['WT Epitope Seq'].replace('', np.nan)
   input_df= input_df.dropna(subset=['WT Epitope Seq'])
   input_df = input_df.assign(Sample = sample_name)
   input_df["HLA Allele"] = input_df["HLA Allele"].str[4:]
   #print(input_df["HLA"])
   input_df["HLA"]=input_df["HLA Allele"].str.replace('[*,:]', '')
   print(input_df["HLA"].map(str))
   input_df["MUTATION_ID"] = input_df["Gene Name"] + "_" + input_df["Chromosome"].map(str) + "_"+ input_df["Peptide Length"].map(str) + "_" +  input_df["Sample"].map(str) +"_" + input_df["HLA"].map(str)
   print(input_df["MUTATION_ID"].map(str))
   #id_values =pd.Series(np.arange(-1,len(input_df)))
   #input_df["ID"]=id_values
   input_df.insert(0, 'ID', range(0, 0 + len(input_df)))
   output_df = input_df[["ID","MUTATION_ID","Sample","WT Epitope Seq","MT Epitope Seq","Corresponding WT Score","Best MT Score","HLA","Peptide Length"]]
   output_df.columns = ['ID','MUTATION_ID','Sample','WT.Peptide','MT.Peptide','WT.Score','MT.Score','HLA','Peptide-Length']
   output_df.to_csv(destination+sample_name+'_processed_output.txt',index=False,sep ='\t')
   return output_df

if __name__=='__main__':
    main()
