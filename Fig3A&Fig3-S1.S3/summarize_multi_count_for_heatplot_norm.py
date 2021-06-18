#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculate log2fc from input count files from different genotypes, designed for 
features including ('Rsp','1.688','20A','flamenco','42AB','38C1','38C2','80F').
optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  name of the protein targeted
  -s STUDY, --study STUDY
                        name of the file with study design info
  -c [CONTROL [CONTROL ...]], --control [CONTROL [CONTROL ...]]
                        list of files from control condition
  -t TREATMENT [TREATMENT ...], --treatment TREATMENT [TREATMENT ...]
                        list of files from treatment
  -o OUTPUT, --output OUTPUT
                        output log2FC file
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""


from argparse import ArgumentParser
from scipy import stats
import math
import numpy as np

p=ArgumentParser(description="This script summarize count output files from different genotypes together, designed for features including ('Rsp','bp','20A','flamenco','42AB','38C1','38C2','80F'), for R analysis.")
#p.add_argument("-f","--feature",nargs='+', help="target feature which you are interested")
p.add_argument("-n", "--name", help="name of the protein targeted")
p.add_argument("-s", "--study", help="name of the file with study design info")
p.add_argument("-c", "--control", nargs='*', help="list of files from control condition")
p.add_argument("-t", "--treatment", nargs='+', help="list of files from treatment")
p.add_argument("-o", "--output", help="output log2FC file")


args=p.parse_args()
#feature=args.feature
control=args.control
treatment=args.treatment
name=args.name
study=args.study
outfile=args.output+'.summary'
countfile=args.output+'.count'

target=['Rsp_SAT','359bp_SAT','353bp_SAT','356bp_SAT','260bp_SAT','piRNA_cluster_20A','piRNA_cluster_flamenco','piRNA_cluster_42AB','piRNA_cluster_38C1','piRNA_cluster_38C2','piRNA_cluster_80F']

def extract_info(file,target):
    
    """
    Function: extract information about target features from .out file from count_proportional to a list.
    Args:
        .out file 
        target: list of features
    Returns:
        summary: dict of lists of information
    """
    summary={}   
    for t in target:
        summary[t]={}
        summary[t]['rpm']=0
               
    with open(file,'r') as ifile:
        for line in ifile.readlines():
            if not line.startswith('#'):
                line=line.strip('\n')
                info=line.split('\t')
                feature=info[0]
                rpm=float(info[1])
                for t in target:
                    if t in feature:
                        summary[t]['rpm']+=rpm
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
with open(study,'r') as ifile:
    for line in ifile.readlines():
        line=line.strip('\n')
        study=line.split('\t',)[0]
        for item in line.split('\t',)[1].split(',',):
            study_info[item]=study

#
study_c=[]
wildtype={}
for c in control:
    c_name=extract_name(c)
    study_c.append(c_name)    
    wildtype[c_name]=extract_info(c,target)

study_t=[] 
mutant={}
for t in treatment:
    t_name=extract_name(t)
    study_t.append(t_name) 
    mutant[t_name]=extract_info(t,target)


summary={}
summary['wildtype']={}
summary['mutant']={}

for i in target:
    summary['wildtype'][i]={}
    summary['mutant'][i]={}
   
    #for WT
    rpm_sum=[]
    for key in wildtype.keys():
        rpm_sum.append(wildtype[key][i]['rpm'])
    rpm_ave=np.mean(rpm_sum)    #mean
    if len(rpm_sum)>1:
        #rpm_sem=np.std(rpm_sum, ddof=1) / math.sqrt(len(rpm_sum))
        rpm_sem=stats.sem(rpm_sum, axis=None, ddof=1)   #standard error
    elif len(rpm_sum)==1:
        rpm_sem='nan'
    summary['wildtype'][i]['rpm']=(rpm_ave, rpm_sem)
    
    #for mutant
    count_sum=[]
    rpm_sum=[]
    for key in mutant.keys():
        rpm_sum.append(mutant[key][i]['rpm'])
    rpm_ave=np.mean(rpm_sum)    #mean
    if len(rpm_sum)>1:
        #rpm_sem=np.std(rpm_sum, ddof=1) / math.sqrt(len(rpm_sum))
        rpm_sem=stats.sem(rpm_sum, axis=None, ddof=1)   #standard error
    elif len(rpm_sum)==1:
        rpm_sem='nan'
    summary['mutant'][i]['rpm']=(rpm_ave, rpm_sem)


with open(outfile,'w') as ofile:
    ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('data_WT','data_mutant','protein','repeat','WT_rpm_average','WT_rpm_sem','Mutant_rpm_average','Mutant_rpm_sem','log2_fold_change'))
    for key in target:
        WT_rpm_average=summary['wildtype'][key]['rpm'][0]
        WT_rpm_sem=summary['wildtype'][key]['rpm'][1]
        Mutant_rpm_average=summary['mutant'][key]['rpm'][0]
        Mutant_rpm_sem=summary['mutant'][key]['rpm'][1]
        fold_change=math.log2((Mutant_rpm_average+0.01)/(WT_rpm_average+0.01))
        data_WT='.'.join(study_c)
        data_mutant='.'.join(study_t)
        ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(data_WT,data_mutant,name,key,WT_rpm_average,WT_rpm_sem,Mutant_rpm_average,Mutant_rpm_sem,fold_change))
        
with open(countfile,'w') as ofile:
    ofile.write("%s\t%s\t%s\t%s\t%s\n" %('study','genotype','dataID','repeat','RPM'))
    for key in target:
        for c_name in study_c:
            genotype='wildtype'
            count=wildtype[c_name][key]['rpm']
            ofile.write("%s\t%s\t%s\t%s\t%s\n" %(study_info[c_name], genotype, c_name, key, count))
        for t_name in study_t:
            genotype=name
            count=mutant[t_name][key]['rpm']
            ofile.write("%s\t%s\t%s\t%s\t%s\n" %(study_info[t_name], genotype, t_name, key, count))
        