#!/bin/bash

"""
This script count the number of reads mapped to each strand of Rsp and 1.688 satellite DNA.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

module load samtools
module load anaconda3/5.2.0

#define your working directory
workingdir=$(pwd)
echo "Working in dir: $workingdir"

REFERENCEGFF="dmel_scaffold2_plus0310_rm_for_htseq.gff3"

target="Rsp"

input=("SRR..."
	   "SRR..."
	  )
	  
for item in "${input[@]}"
do 
	echo "working on ${item}..."
	
	#pull out reads from one strand
	samtools view -h -b -f 16 "${item}.mapped.sorted.bam" > "${item}.16.mapped.sorted.bam"
	samtools index "${item}.16.mapped.sorted.bam"
	
	python3 extract_sequence_by_feature_gff.py -t "Rsp_SAT" -s RF "${item}.16.mapped.sorted.bam" "$REFERENCEGFF" "${item}_Rsp_16_trimmed.fq"
	python3 extract_sequence_by_feature_gff.py -t "bp_SAT" -s RF "${item}.16.mapped.sorted.bam" "$REFERENCEGFF" "${item}_1.688_16_trimmed.fq"	

	N=$(($(cat "${item}_Rsp_16_trimmed.fq" | wc -l) /4))
	M=$(($(cat "${item}_1.688_16_trimmed.fq" | wc -l) /4))
	echo "Rsp plus strand: $N; 1.688 plus strand: $M"
	
	#pull out reads from the other strand
	samtools view -h -b -F 16 "${item}.mapped.sorted.bam" > "${item}.0.mapped.sorted.bam"
	samtools index "${item}.0.mapped.sorted.bam"

	python3 extract_sequence_by_feature_gff.py -t "Rsp_SAT" -s RF "${item}.0.mapped.sorted.bam" "$REFERENCEGFF" "${item}_Rsp_0_trimmed.fq"
	python3 extract_sequence_by_feature_gff.py -t "bp_SAT" -s RF "${item}.0.mapped.sorted.bam" "$REFERENCEGFF" "${item}_1.688_0_trimmed.fq"	

	Nm=$(($(cat "${item}_Rsp_0_trimmed.fq" | wc -l) /4))
	Mm=$(($(cat "${item}_1.688_0_trimmed.fq" | wc -l) /4))
	echo "Rsp minus strand: $Nm; 1.688 minus strand: $Mm"
	
	###
	#unique
	###
	#pull out reads from one strand
	samtools view -h -b -f 16 "${item}.unique.mapped.sorted.bam" > "${item}.16.unique.mapped.sorted.bam"
	samtools index "${item}.16.unique.mapped.sorted.bam"

	python3 extract_sequence_by_feature_gff.py -t "Rsp_SAT" -s RF "${item}.16.unique.mapped.sorted.bam" "$REFERENCEGFF" "${item}_Rsp_16_trimmed.unique.fq"
	python3 extract_sequence_by_feature_gff.py -t "bp_SAT" -s RF "${item}.16.unique.mapped.sorted.bam" "$REFERENCEGFF" "${item}_1.688_16_trimmed.unique.fq"	

	N=$(($(cat "${item}_Rsp_16_trimmed.unique.fq" | wc -l) /4))
	M=$(($(cat "${item}_1.688_16_trimmed.unique.fq" | wc -l) /4))
	echo "Rsp plus strand unique: $N; 1.688 plus strand unique: $M"
	
	#pull out reads from the other strand
	samtools view -h -b -F 16 "${item}.unique.mapped.sorted.bam" > "${item}.0.unique.mapped.sorted.bam"
	samtools index "${item}.0.unique.mapped.sorted.bam"

	python3 extract_sequence_by_feature_gff.py -t "Rsp_SAT" -s RF "${item}.0.unique.mapped.sorted.bam" "$REFERENCEGFF" "${item}_Rsp_0_trimmed.unique.fq"
	python3 extract_sequence_by_feature_gff.py -t "bp_SAT" -s RF "${item}.0.unique.mapped.sorted.bam" "$REFERENCEGFF" "${item}_1.688_0_trimmed.unique.fq"	

	Nm=$(($(cat "${item}_Rsp_0_trimmed.unique.fq" | wc -l) /4))
	Mm=$(($(cat "${item}_1.688_0_trimmed.unique.fq" | wc -l) /4))
	echo "Rsp minus strand unique: $Nm; 1.688 minus strand unique: $Mm"

done