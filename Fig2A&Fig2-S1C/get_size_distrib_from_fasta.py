#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script takes fasta files as input, and summarize the length distribution.
To get the size distribution of a specific repeat element (eg. Rsp), use 
"extract_sequence_by_feature_gff.py" to output fastq format file, and then 
reformate it into fasta file as input for this script.

use -h for the help page with different optional arguments you can choose with this script.
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        list of input fasta file
  -o OUTPUT, --output OUTPUT
                        output file name, output file is tab deliminated file
                        with two columns of size and its percentage

@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

from argparse import ArgumentParser

p=ArgumentParser(description="This script summarize sequence length from fasta files.")
p.add_argument("-i", "--input", nargs='+', help="list of input fasta file")
p.add_argument("-o", "--output", help="output file name, output file is tab deliminated file with two columns of size and its percentage")

args=p.parse_args()
infile=args.input
outfile=args.output

summary={}
total_count={}
for file in infile:
    name=file.split('_',)[0]
    total=0
    summary[name]={}
    total_count[name]=0
    with open(file,'r') as ifile:
        for line in ifile:
            if not line.startswith('>'):
                line=line.strip('\n')
                length=len(line)
                if length not in summary[name].keys():
                    summary[name][length]=0
                summary[name][length]+=1
                total_count[name]+=1

#sort size key and input file names
keys=[]
iput=[]
for name in summary.keys():
    iput.append(name)
    for key in summary[name].keys():
        if key not in keys:
            keys.append(key)
keys.sort()
iput.sort()

#keep record of the order
print (iput)

with open(outfile,'w') as ofile:
    #ofile.write("%s\t%s\n" %('#size(nt)','#percentage(%)'))
    for name in iput:
        for key in keys:
            if key in summary[name].keys():
                count=100*summary[name][key]/total_count[name]
            else:
                count=0
            ofile.write("%d\t%.2f\n" %(key,count))
    
