#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

from argparse import ArgumentParser

p=ArgumentParser(description="This script summarize count output files from different genotypes together, all features, for DESeq2 input.")
#p.add_argument("-f","--feature",nargs='+', help="target feature which you are interested")
p.add_argument("-n", "--name", help="target protein in the study")
p.add_argument("-s", "--study", help="tab delaminated file listing the study coresponding to each data, used when data come from different studies. 1st column is study name, 2nd column is data name separated by comma.")
p.add_argument("-c", "--control", nargs='*', help="list of files from control condition")
p.add_argument("-t", "--treated", nargs='+', help="list of files from treated condition")
p.add_argument("-o", "--output", help="output file")

args=p.parse_args()
#feature=args.feature
name=args.name
study_file=args.study
control=args.control
treatment=args.treated
out=args.output


def extract_info(file):
    
    """
    Function: extract information about target features from .out file from count_proportional to a list.
    Args:
        .out file 
        target: list of features
    Returns:
        summary: dict of lists of information
    """
    
    summary={}           
    with open(file,'r') as ifile:
        for line in ifile.readlines():
            if not line.startswith('#'):
                line=line.strip('\n')
                info=line.split('\t')
                feature=info[0]
                count=float(info[1])
                #rpm=float(info[2])
                summary[feature]=count

    return summary


def extract_name(file):

    """
    Function: extract file name from full path.
    Args:
        file_path
    Returns:
        file_name
    """
    if '/' in file:
        file_name=file.rsplit('/',1)[1]
        name=file_name.split('.',)[0]
    else:
        name=file.split('.',)[0]
    
    return name

 
    
#main

#extract information about data-study association
study_info={}
with open(study_file,'r') as ifile:
    for line in ifile.readlines():
        line=line.strip('\n')
        study=line.split('\t',)[0]
        for item in line.split('\t',)[1].split(',',):
            study_info[item]=study

#extract count information
summary={}

study_c=[]
for c in control:
    c_name=extract_name(c)
    #keep track of controls
    study_c.append(c_name)
    #summarize count
    wildtype=extract_info(c)
    for feature in wildtype.keys():
        if feature not in summary.keys():
            summary[feature]={}
        summary[feature][c_name]=wildtype[feature]
        
study_t=[] 
for t in treatment:
    t_name=extract_name(t)
    #keep track of controls
    study_t.append(t_name)
    #summarize count
    treated=extract_info(t)
    for feature in treated.keys():
        if feature not in summary.keys():
            summary[feature]={}
        summary[feature][t_name]=treated[feature]

#write count output
out_count=out+".countData"
with open(out_count,'w') as ofile:
    clist='\t'.join(study_c)
    tlist='\t'.join(study_t)
    ofile.write("%s\t%s\n" %(clist,tlist))
    for feature in summary.keys():
        count_c=[]
        count_t=[]
        for study in study_c:
            count=summary[feature][study]
            count_c.append(str(round(count)))
        count_c_list='\t'.join(count_c)
        for study in study_t:
            count=summary[feature][study]
            count_t.append(str(round(count)))
        count_t_list='\t'.join(count_t)
        ofile.write("%s\t%s\t%s\n" %(feature,count_c_list,count_t_list))

#write col info output    
out_col=out+".colData"
with open(out_col,'w') as ofile:
    ofile.write("%s\t%s\t%s\n" %('condition','study','name'))
    for study in study_c:
        info=study_info[study] 
        ofile.write("%s\t%s\t%s\t%s\n" %(study,'wildtype',info, name))
    for study in study_t:
        info=study_info[study]
        ofile.write("%s\t%s\t%s\t%s\n" %(study,'mutant',info, name))
    
