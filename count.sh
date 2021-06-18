#!/bin/bash

"""
This script count and extract reads that mapped to repeats in a annotation gff file (eg. Rsp, 1.688 satellite) from alignment bam files.
example usage: count.sh 
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

module load samtools/1.5
module load anaconda3/5.2.0

#define your working directory
workingdir=""
cd $workingdir || ( echo "can't open the $workingdir"; exit )
echo "Working in dir: $workingdir"

#this is your list of read files with full paths to reads
#only put the prefix, not full read file names
#reads="SRR..."
FASTQLIST=$workingdir/FASTQLIST.txt

#this is the annotation gff file
REFERENCEGFF="dmel_scaffold2_plus0310_rm_for_htseq.gff3"

mkdir $workingdir/extract_reads
mkdir $workingdir/count

#read in each alignment file
while FASTQLIST='' read -r reads || [[ -n "$reads" ]]; do
	
	name=$(echo "$reads" | awk -F'/' '{print $(NF-1)}')
	samplename=$(echo "$reads" | awk -F'/' '{print $NF}') 
	dir_path=$(echo "$reads" | awk -F'/' '{OFS="/"; $NF=""; print $0}')

	echo "Working on prefix $name and read files $samplename"
	
	echo "read path is: $reads"
	
	output=$workingdir/${samplename}
			
	N=$(samtools view "${output}.mapped.sorted.bam" | wc -l)
	echo "total mapped reads: $N."


###
#count number of reads mapped to repeats in the gff file using customized python script
###
	
	count_out=$workingdir/count/${samplename}
	
	python3 htseq_bam_count_proportional.py "${output}.mapped.sorted.bam" "$REFERENCEGFF" "${count_out}.out"


###
#normalize read counts to miRNA/flam
###
	
	#normalize to miRNA
	norm_factor="miRNA"
	norm=$(samtools view "${output}.miRNA.bam" | wc -l)
	echo "normalize to $norm_factor: $norm"
	awk -v var="$norm" -F"\t" 'NR==1 {OFS="\t"; print$1, $2"_norm_by_miRNA"}; NR>1 {OFS="\t"; print $1, $2*1000000/var}' "${count_out}.out" >  "${count_out}.out_norm_${norm_factor}"
	
	#normalize to flam
	norm_factor="flam"
	norm=$(egrep "flamenco" "${count_out}.out" | awk 'BEGIN{sum=0} {sum+=$2} END {print sum}')
	echo "normalize to $norm_factor: $norm"
	awk -v var="$norm" -F"\t" 'NR==1 {OFS="\t"; print$1, $2"_norm_by_flam"}; NR>1 {OFS="\t"; print $1, $2*1000000/var}' "${count_out}.out" >  "${count_out}.out_norm_${norm_factor}"

	
###
#extract reads mapped to a specific repeat using customized python script
###

	extract_out=$workingdir/extract_reads/${samplename}
	
	python3 extract_sequence_by_feature_gff.py -t "Rsp_SAT" "${output}.mapped.sorted.bam" "$REFERENCEGFF" "${extract_out}_Rsp_trimmed.fq"
	python3 extract_sequence_by_feature_gff.py -t "bp_SAT" "${output}.mapped.sorted.bam" "$REFERENCEGFF" "${extract_out}_1.688_trimmed.fq"	
	gzip -r $workingdir/extract_reads/
	
	#convert to fasta format
	zcat "${extract_out}_Rsp_trimmed.fq.gz" | grep -A1 "^[>@]" | sed 's/@/>/g'| grep -v "^--" > "${extract_out}_Rsp.fasta"
	zcat "${extract_out}_1.688_trimmed.fq.gz" | grep -A1 "^[>@]" | sed 's/@/>/g' | grep -v "^--" > "${extract_out}_1.688.fasta"

done < "$FASTQLIST"
echo "DONE!"

