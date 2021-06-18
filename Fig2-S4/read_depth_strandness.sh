#!/bin/bash

"""
This script get reads alignment depth for each strand of Rsp major locus and 1.688 2L locus (260-bp).
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

module load samtools

workingdir=$(pwd)
echo "Working in dir: $workingdir"

###
#Rsp
###

target="Rsp_major"

input=("SRR..."
	   "SRR..."
	  )
	  
for item in "${input[@]}"
do
	#pull out chromosome region of interest
	samtools view -h "${item}.mapped.sorted.bam" 2R_16:1-168170 > "${item}.${target}.mapped.sorted.bam"
	
	#get depth for one strand
	samtools view -h -b -f 16 "${item}.${target}.mapped.sorted.bam" > "${item}.${target}.16.mapped.sorted.bam"
	samtools index "${item}.${target}.16.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.16.mapped.sorted.bam" -r 2R_16:1-168170 > "${item}.${target}.16.depth"
	
	#get depth for the other strand
	samtools view -h -b -F 16 "${item}.${target}.mapped.sorted.bam" > "${item}.${target}.0.mapped.sorted.bam"
	samtools index "${item}.${target}.0.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.0.mapped.sorted.bam" -r 2R_16:1-168170 > "${item}.${target}.0.depth"
	
	#combine
	awk '{print $1, $2, $3, "plus"}' "${item}.${target}.16.depth" > "${item}.${target}.depth"
	awk '{print $1, $2, "-"$3, "minus"}' "${item}.${target}.0.depth" >> "${item}.${target}.depth"

	###
	#unique
	###
	#pull out chromosome region of interest
	samtools view -h "${item}.unique.mapped.sorted.bam" 2R_16:1-168170 > "${item}.${target}.unique.mapped.sorted.bam"
	
	#get depth for one strand
	samtools view -h -b -f 16 "${item}.${target}.unique.mapped.sorted.bam" > "${item}.${target}.16.unique.mapped.sorted.bam"
	samtools index "${item}.${target}.16.unique.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.16.unique.mapped.sorted.bam" -r 2R_16:1-168170 > "${item}.${target}.16.unique.depth"
	
	#get depth for the other strand
	samtools view -h -b -F 16 "${item}.${target}.unique.mapped.sorted.bam" > "${item}.${target}.0.unique.mapped.sorted.bam"
	samtools index "${item}.${target}.0.unique.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.0.unique.mapped.sorted.bam" -r 2R_16:1-168170 > "${item}.${target}.0.unique.depth"
	
	#combine
	awk '{print $1, $2, $3, "plus"}' "${item}.${target}.16.unique.depth" > "${item}.${target}.unique.depth"
	awk '{print $1, $2, "-"$3, "minus"}' "${item}.${target}.0.unique.depth" >> "${item}.${target}.unique.depth"

done

###
#260bp
###

target="260bp"

input=("SRR..."
	   "SRR..."
	  )
	  
for item in "${input[@]}"
do 
	#pull out chromosome region of interest
	samtools view -h "${item}.mapped.sorted.bam" 2L_2:402701-460225 > "${item}.${target}.mapped.sorted.bam"
	
	#get depth for one strand	
	samtools view -h -b -f 16 "${item}.${target}.mapped.sorted.bam" > "${item}.${target}.16.mapped.sorted.bam"
	samtools index "${item}.${target}.16.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.16.mapped.sorted.bam" -r 2L_2:402701-460225 > "${item}.${target}.16.depth"
	
	#get depth for the other strand
	samtools view -h -b -F 16 "${item}.${target}.mapped.sorted.bam" > "${item}.${target}.0.mapped.sorted.bam"
	samtools index "${item}.${target}.0.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.0.mapped.sorted.bam" -r 2L_2:402701-460225 > "${item}.${target}.0.depth"
	
	#combine
	awk '{print $1, $2, $3, "plus"}' "${item}.${target}.16.depth" > "${item}.${target}.depth"
	awk '{print $1, $2, "-"$3, "minus"}' "${item}.${target}.0.depth" >> "${item}.${target}.depth"

	###
	#unique
	###
	#pull out chromosome region of interest
	samtools view -h ../"${item}.unique.mapped.sorted.bam" 2L_2:402701-460225 > "${item}.${target}.unique.mapped.sorted.bam"
	
	#get depth for one strand	
	samtools view -h -b -f 16 "${item}.${target}.unique.mapped.sorted.bam" > "${item}.${target}.16.unique.mapped.sorted.bam"
	samtools index "${item}.${target}.16.unique.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.16.unique.mapped.sorted.bam" -r 2L_2:402701-460225 > "${item}.${target}.16.unique.depth"
	
	#get depth for the other strand	
	samtools view -h -b -F 16 "${item}.${target}.unique.mapped.sorted.bam" > "${item}.${target}.0.unique.mapped.sorted.bam"
	samtools index "${item}.${target}.0.unique.mapped.sorted.bam"
	samtools depth -a "${item}.${target}.0.unique.mapped.sorted.bam" -r 2L_2:402701-460225 > "${item}.${target}.0.unique.depth"
	
	#combine
	awk '{print $1, $2, $3, "plus"}' "${item}.${target}.16.unique.depth" > "${item}.${target}.unique.depth"
	awk '{print $1, $2, "-"$3, "minus"}' "${item}.${target}.0.unique.depth" >> "${item}.${target}.unique.depth"

done
