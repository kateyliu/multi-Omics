#!/usr/bin/env python 
"""
This Script is for evaluation of germline matching situation of samples, 
The input parameter is the vcf compare file which is generated from vcf-compare.
"""
import os
import sys
import re 
import numpy as np
from optparse import OptionParser

def main():
	usage = "USAGE: %prog -i [fcvcompare] "
	optparser = OptionParser(usage=usage)
	optparser.add_option("-i", "--input", help="path to vcfcompare file")
	(options, args) = optparser.parse_args(sys.argv)
	if not options.input :
		optparser.print_help()
		sys.exit(-1)
	Input_file = options.input
	#print(Input_file)  
	f = open(Input_file, 'r')
	lines = f.readlines()
	for line in lines:
	    line  = line.strip()
	    p1 = re.compile(r'[(](.*?)[)]', re.S) 
	    bracket = re.findall(p1, line)
	    if line.startswith('VN') and len(bracket) == 2:
	        bracket_list = [float(i.rstrip('%'))/100 for i in bracket]
	        #print(bracket_list)
	        bracket_array = np.array(bracket_list)
	        if sum(bracket_array > 0.9) > 1:
	            print( 'match')
	        else:
	            print('mismatch')
	            sys.exit()

if __name__ == '__main__':
	#print('exe')
    main()