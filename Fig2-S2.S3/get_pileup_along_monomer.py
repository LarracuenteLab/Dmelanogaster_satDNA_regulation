#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script gets the read distribution along a repeat monomer from blast output file.
example usage: python3 get_pileup_along_monomer.py -l repeat_monomer_length -i "Rsp_vs_consensus.blast" -o "Rsp_vs_monomer.distrib"
use -h for the help page with different optional arguments you can choose with this script.
  -h, --help            show this help message and exit
  -l LENGTH, --length LENGTH
                        repeat monomer length
  -n TOTAL, --total TOTAL
                        total mapped reads
  -i INPUT, --input INPUT
                        input .blast file
  -o OUTPUT, --output OUTPUT
                        output .distrib file
@author: xiaolu_wei@urmc.rochester.edu
"""

from argparse import ArgumentParser

p=ArgumentParser(description="This script summarize count distribution (percentage of reads that are extracted by this repeat) along a monomer from blastn output.")
p.add_argument("-l", "--length", type=int, help="repeat monomer length")
p.add_argument("-n", "--total", type=int, help="total mapped reads")
p.add_argument("-i", "--input", help="input .blast file")
p.add_argument("-o", "--output", help="output .distrib file")

args=p.parse_args()
length=args.length
total=args.total
infile=args.input
outfile=args.output

#pre-set count to 0 at all positions
pipos={}
sipos={}
allpos={}
for i in range(1,length+1):
    pipos[i]=0
    sipos[i]=0
    allpos[i]=0

#strand info
strand={}
strand['antisense']=0
strand['sense']=0
    
#count pile up at positions
summary={}
check={}
with open(infile,'r') as ifile:
    for line in ifile:
        line=line.strip('\n')
        info=line.split('\t')
        #only count each read once
        if not info[0] in check.keys():
            check[info[0]]='counted'
            #record strand information
            if int(info[8]) > int(info[9]):
                start=int(info[9])
                end=int(info[8])
                strand['antisense']+=1
            else:
                start=int(info[8])
                end=int(info[9])
                strand['sense']+=1
            
            #aligned to di/trimer, now convert alignment position to monomer position
            align=[]
            #align to the first monomer
            if end <= length:
                for i in range(start,end+1):
                    align.append(i)
            #align to junction
            elif start < length:
                for i in range(start, length+1):
                    align.append(i)
                for j in range(1, end-length):
                    align.append(j)
            #align to the second monomer
            elif start > length:
                for i in range(start-length, end-length+1):
                    align.append(i)
            
            #define pi/siRNA based on length
            if int(info[3]) <= 23:
                feature='siRNA'
                #count
                for i in align:
                    sipos[i]+=1
            else:
                feature='piRNA'
                #count
                for i in align:
                    pipos[i]+=1
            #count all
            for i in align:
                allpos[i]+=1
        else:
            print (info[0], "read already counted before, ignore now...")

#record strand information
#total number of rads from input
allReads = len(check)
sense = strand['sense'] / allReads
antisense = strand['antisense'] / allReads

#normalize
if total:
    #ouput RPM value normalized by total mapped reads
    norm=total/1000000
    note="RPM"
else: 
    #output percentage of reads from input fasta
    #normalize count by total number of reads that are mapped to this element
    norm=len(check)/100
    note="% of reads"

#write to output
with open(outfile,'w') as ofile:
    ofile.write("#sense strand proportion: %.2f, antisense strand proportion: %.2f\n" %(sense, antisense))
    ofile.write("#normalization method: %s\n" %(note))
    ofile.write("%s\t%s\t%s\t%s\n" %('position', 'piRNA_count', 'siRNA_count', 'total_count'))
    for i in range(1,length+1):
        position=i
        #norm
        piRNA_count=pipos[i]/norm
        siRNA_count=sipos[i]/norm
        total_count=allpos[i]/norm
        ofile.write("%d\t%.2f\t%.2f\t%.2f\n" %(position, piRNA_count, siRNA_count, total_count))


               
        
