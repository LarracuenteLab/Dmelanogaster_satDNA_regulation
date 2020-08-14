#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script extract reads that map to a repeat feature of interest, output fastq format.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

"""
bam file or sam file need to be sorted by leftmost coordinates
"""

from argparse import ArgumentParser
import pysam


p=ArgumentParser(description="This script output feature counts calculated from aignment file and gff file.")
p.add_argument("-f", "--format", help="alignment file is bam or sam format", choices=["bam", "sam"], default="bam")
p.add_argument("-i", "--attribute", help="GFF attribute to be used as feature ID. Provide it only if there are multiple ID attributes in the last column and the ID attributes are separated by ';'.", default="")
p.add_argument("-t","--target",help="target feature which you are interested in extracting sequence reads mapped to it")
p.add_argument("-s","--stranded",help="if the sequencing is stranded or not", choices=["RF","FR","R","F"])
p.add_argument("alignment_file", help="alignment file, sam or bam")
p.add_argument("gff_file", help="gff file")
p.add_argument("output", help="output file")


args=p.parse_args()
form=args.format
attribute=args.attribute
target=args.target
stranded=args.stranded
alignment_file=args.alignment_file
gff=args.gff_file
out=args.output


def reverse_complementary(sequence):
    """
    """
    Rsequence=sequence[::-1]
    RCsequence=[]
    for i in Rsequence:
        if i=="A":
            newi="T"
        elif i=="T" or i=="U":
            newi="A"
        elif i=="C":
            newi="G"
        elif i=="G":
            newi="C"
        elif i=="N":
            newi="N"
        else:
            print ("Error: can't recognize sequence: "+sequence)
        RCsequence.append(newi)
    newsequence=''.join(RCsequence)
    
    return(newsequence)
                

def get_gff_ID(entry, attribute):
    
    """
    Function: parse the last column of gff file to get ID of the entry
              for gff from repeatmasker 

    Args: entry: the last column of a gff line
          attribute: attribute to be used as feature ID
    Returns: ID
    """
    
    # if attribute is not provided
    
    if attribute=="": 
        # repeatmasker gff2 (Target "Motif:TAHRE" 4817 10457)
        if '"' in entry:
            ID=entry.split(' ',)[1].strip('"').split(':',)[1]
        
        # rmOutToGFF3 gff3
        else:
            # (TAHRE 4817 10457)
            if '=' not in entry:
                ID=entry.split(' ',)[0]
        
            #(name=TAHRE 4817 10457) strip off "name=" in annotations
            else:
                ID=entry.split(' ',)[0].split('=',)[1]
    
    # attribute is provided, parse based on it
    
    else:
        # (ID=TAHRE.2L_1-2-5132;geneID=TAHRE)
        IDs=entry.split(';',)
        for item in IDs:
            item=item.strip(' ')
            # only take the attribute user provided
            if item.startswith(attribute):
                # strip off =, white space, or double quote.
                ID=item.strip(attribute).strip(' ="')
                    
    return ID
 
    

def extract_reads_bam(alignment_file, gff, target, stranded):
    
    """
    Function: read in alignment bam file, repeat annotation gff3 file and target feature, extract reads mapped to the target feature.
              overlap>=50%
    Args:
        alignment_file 
        gff
        target feature
        stranded info
    Returns:
        summary: dictionay of read sequences
    """
    readfile=pysam.AlignmentFile(alignment_file, "rb")
    
    # total mapped reads number
    #total=readfile.mapped
    
    summary={}
    with open(gff, 'r') as gfile:
        lines=gfile.readlines()
        for line in lines:
            if not line.startswith('#'):
                line=line.strip('\n')
                #skip empty lines
                if not line=='':
                    info=line.split('\t',)
                    contig=info[0]
                    #typ=info[2]
                    start=int(info[3])
                    end=int(info[4])
                    
                    # get ID name for each line
                    ID=get_gff_ID(info[8].strip(' '), attribute)
                    
                    if target in ID:
                        name='.'.join([ID,contig,str(start),str(end)])
                        summary[name]=[]
                        
                        # fetch reads overlapped with this entry
                        # pysam fetch back reads overlapped to the feature
                        for read in readfile.fetch(contig, start, end):
                            # read: pysam AlignedSegment object
                            read_length=read.query_length
                            read_start=(read.reference_start)+1
                            read_end=read.reference_end
                            
                            #if stranded sequence, get the sequence with the same orientation as the actual RNAs
                            if stranded:
                                #if first strand
                                if stranded=='RF':
                                    if read.is_read1:
                                        sequence=reverse_complementary(read.query_sequence)
                                        quality=read.query_qualities[::-1]
                                    elif read.is_read2:
                                        sequence=read.query_sequence
                                        quality=read.query_qualities
                                    else:
                                        print ("Error: not read1 or read2")
                                elif stranded=='FR':
                                    if read.is_read2:
                                        sequence=reverse_complementary(read.query_sequence)
                                        quality=read.query_qualities[::-1]
                                    elif read.is_read1:
                                        sequence=read.query_sequence
                                        quality=read.query_qualities
                                    else:
                                        print ("Error: not read1 or read2")
                                elif stranded=='R':
                                    sequence=reverse_complementary(read.query_sequence)
                                    quality=read.query_qualities[::-1]
                                elif stranded=='F':
                                    sequence=read.query_sequence
                                    quality=read.query_qualities
                            else:
                                sequence=read.query_sequence
                                quality=read.query_qualities
                            
                            #convert quality array to string
                            tmp=[]
                            for q in quality:
                                tmp.append(str(q))
                            newquality=','.join(tmp)
                            
                            #assign reads
                            if read_start >= start and read_end <= end:
                                summary[name].append((sequence,newquality))
                            else:
                                if read_start > start:
                                    overlap_start=read_start
                                else:
                                    overlap_start=start
                                if read_end < end:
                                    overlap_end=read_end
                                else:
                                    overlap_end=end                        
                                overlap_length=overlap_end-overlap_start+1
                                if overlap_length/read_length>=0.5:
                                    summary[name].append((sequence,newquality))

    readfile.close()
    
    return summary

# main

summary=extract_reads_bam(alignment_file, gff, target, stranded)

with open(out,'w') as ofile:
    for key in summary:
        for pair in summary[key]:
            index=summary[key].index(pair)+1
            header='.'.join([key,str(index)])
            seq=pair[0]
            qual=pair[1]
            ofile.write("@%s\n%s\n+%s\n%s\n" %(header,seq, header, qual))
