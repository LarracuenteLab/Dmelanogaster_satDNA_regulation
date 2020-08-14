#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script count reads mapped to repeat features from alignment files.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

"""
bam file or sam file need to be sorted by leftmost coordinates
"""

from argparse import ArgumentParser
import pysam


p=ArgumentParser(description="This script output feature counts calculated from aignment file and gff file.")
p.add_argument("-p", "--paired", help="sequencing reads are paired end", action="store_true")
p.add_argument("-f", "--format", help="alignment file is bam or sam format", choices=["bam", "sam"], default="bam")
p.add_argument("-i", "--attribute", help="GFF attribute to be used as feature ID. Provide it only if there are multiple ID attributes in the last column and the ID attributes are separated by ';'.", default="")
p.add_argument("-s", "--sep_size", help="output two files for si/piRNA based on read length", action="store_true")
p.add_argument("alignment_file", help="alignment file, sam or bam")
p.add_argument("gff_file", help="gff file")
p.add_argument("output", help="output file")


args=p.parse_args()
paired=args.paired
form=args.format
attribute=args.attribute
sep_size=args.sep_size
alignment_file=args.alignment_file
gff=args.gff_file
out=args.output


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
        # (gene_id "MSTRG.1"; transcript_id "FBtr0078103"; ref_gene_id "FBgn0031216";)
        IDs=entry.split(';',)
        for item in IDs:
            item=item.strip(' ')
            # only take the attribute user provided
            if item.startswith(attribute):
                # strip off =, white space, or double quote.
                ID=item.strip(attribute).strip(' ="')
                    
    return ID

 
def get_read_identity(entry):
    
    """
    Function: get read identity based on read length 
    Args: entry: read length
    Returns: identity
    """
    
    # if attribute is not provided
    
    if entry>23:
        identity='piRNA'
    else:
        identity='siRNA'                  
    
    return identity
 
    
def summary_single_end_bam_sep_size(alignment_file, gff):
    
    """
    Function: read in alignment bam file (single end reads) and repeat annotation gff3 file, count reads mapped to repeats.
              count proportional to overlap
    Args:
        alignment_file 
        gff
    Returns:
        count: dictionay of ID and its count
    """
    readfile=pysam.AlignmentFile(alignment_file, "rb")
    
    # total mapped reads number
    total=readfile.mapped
    
    count={}
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
                    
                    if not ID in count.keys():
                        count[ID]={}
                        count[ID]['siRNA']=0
                        count[ID]['piRNA']=0
                    
                    # fetch reads overlapped with this entry
                    tmp_count={}
                    tmp_count['siRNA']=0
                    tmp_count['piRNA']=0
                             
                    # pysam fetch back reads overlapped to the feature
                    for read in readfile.fetch(contig, int(start), int(end)):
                        # read: pysam AlignedSegment object
                        
                        # infer read identity based on read length
                        alignment_length=read.query_alignment_length
                        identity=get_read_identity(alignment_length)
                        
                        # infer count proportional to overlap length
                        read_length=read.query_length
                        """
                        if read.is_reverse:
                            strandness='reverse'
                        else:
                            strandness='forward'
                        """
                        read_start=(read.reference_start)+1
                        read_end=read.reference_end
                        if read_start >= start and read_end <= end:
                            tmp_count[identity]+=1
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
                            tmp_count[identity]+=(overlap_length/read_length)
                    
                    count[ID]['siRNA']+=tmp_count['siRNA']
                    count[ID]['piRNA']+=tmp_count['piRNA']

    readfile.close()
    
    return count,total
    

def summary_single_end_bam(alignment_file, gff):
    
    """
    Function: read in alignment bam file (single end reads) and repeat annotation gff3 file, count reads mapped to repeats.
              count proportional to overlap
    Args:
        alignment_file 
        gff
    Returns:
        count: dictionay of ID and its count
    """
    readfile=pysam.AlignmentFile(alignment_file, "rb")
    
    # total mapped reads number
    total=readfile.mapped
    
    count={}
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
                    
                    if not ID in count.keys():
                        count[ID]=0
                    
                    # fetch reads overlapped with this entry
                    tmp_count=0
                    # pysam fetch back reads overlapped to the feature
                    for read in readfile.fetch(contig, int(start), int(end)):
                        # read: pysam AlignedSegment object
                        read_length=read.query_length
                        """
                        if read.is_reverse:
                            strandness='reverse'
                        else:
                            strandness='forward'
                        """
                        read_start=(read.reference_start)+1
                        read_end=read.reference_end
                        if read_start >= start and read_end <= end:
                            tmp_count+=1
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
                            tmp_count+=(overlap_length/read_length)
                    count[ID]+=tmp_count

    readfile.close()
    
    return count,total


def summary_paired_end_bam(alignment_file, gff):
    
    """
    Function: read in alignment bam file (paired end reads) and repeat annotation gff3 file, count reads mapped to repeats.
              count proportional to overlap
    Args:
        alignment_file 
        gff
    Returns:
        count:  
    """
    readfile=pysam.AlignmentFile(alignment_file, "rb")
    
    # total mapped reads number
    total=readfile.mapped

    count={}
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
                    
                    if not ID in count.keys():
                        count[ID]=0
                    
                    # pysam fetch back reads overlapped to the feature
                    tmp_count=0
                    for read in readfile.fetch(contig, int(start), int(end)):
                        # read: pysam AlignedSegment object                    
                        
                        if read.is_proper_pair:
                            number=0.5
                        else:
                            number=1
                        
                        read_length=read.query_length                    
                        read_start=(read.reference_start)+1
                        read_end=read.reference_end
                        if read_start >= start and read_end <= end:
                            tmp_count+=number
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
                            tmp_count+=number*(overlap_length/read_length)   
                            
                    count[ID]+=tmp_count                     

    readfile.close()    
    return count,total

# main
#output two files for si/piRNA separately
if sep_size:
    count=summary_single_end_bam_sep_size(alignment_file, gff)[0]
    total=summary_single_end_bam_sep_size(alignment_file, gff)[1]
    
    out_1=out+'_siRNA'
    with open(out_1,'w') as ofile:
        ofile.write("%s\t%s\t%s\n" %('#repeat','#count','#RPM')) 
        for key in count.keys():
            repeat=key
            raw=count[key]['siRNA']
            rpm=1000000*raw/total 
            ofile.write("%s\t%.2f\t%.2f\n" %(repeat,raw,rpm))
    out_2=out+'_piRNA'
    with open(out_2,'w') as ofile:
        ofile.write("%s\t%s\t%s\n" %('#repeat','#count','#RPM')) 
        for key in count.keys():
            repeat=key
            raw=count[key]['piRNA']
            rpm=1000000*raw/total 
            ofile.write("%s\t%.2f\t%.2f\n" %(repeat,raw,rpm))

#output one file for all
else:
    if paired:
        count=summary_paired_end_bam(alignment_file, gff)[0]
        total=summary_paired_end_bam(alignment_file, gff)[1]
    else:
        count=summary_single_end_bam(alignment_file, gff)[0]
        total=summary_single_end_bam(alignment_file, gff)[1]    
    
    with open(out,'w') as ofile:
        ofile.write("%s\t%s\t%s\n" %('#repeat','#count','#RPM')) 
        for key in count:
            repeat=key
            raw=count[key]
            rpm=1000000*raw/total 
            ofile.write("%s\t%.2f\t%.2f\n" %(repeat,raw,rpm))
