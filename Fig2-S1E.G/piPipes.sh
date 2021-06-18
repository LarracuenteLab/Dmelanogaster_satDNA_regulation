#!/bin/bash

"""
This script uses piPipe ((Han, et al. 2015)) to calculate ping-pong score for Rsp and 1.688 satellite DNAs for Figure 2-figure supplement 1E&G.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

module load r/3.5.0/b1
module load bedtools
module load samtools/1.5
module load bowtie/1.2.1.1 
module load anaconda3/5.2.0

#define your working directory
workingdir=""
cd $workingdir || ( echo "can't open the $workingdir"; exit )
echo "Working in dir: $workingdir"

REFERENCEINDEX="dmel_scaffold2_plus0310.fasta"
REFERENCEGFF="dmel_scaffold2_plus0310_rm_for_htseq.gff3"

#This is the folder with pingpong calculation scripts
cal_pingpong_dir="/softwares/piPipes/cal_pingpong_scripts"

#this is your list of read files with full paths to reads
#only put the prefix, not full read file names
#reads="SRR5445299"
FASTQLIST=$workingdir/FASTQLIST.txt


if [ -f "$REFERENCEINDEX.1.ebwt" ]; then
	echo "Reference index exists. Moving on..."
else
	echo "making index of reference: $REFERENCEINDEX"
	#index the reference genome if you haven't already
	bowtie-build -f $REFERENCEINDEX $REFERENCEINDEX
fi


targets=('Rsp'
		'1.688'
		)

mkdir ${workingdir}/extract_reads
mkdir ${workingdir}/pingpong_analysis

while FASTQLIST='' read -r reads || [[ -n "$reads" ]]; do

	name=$(echo "$reads" | awk -F'/' '{print $(NF-1)}')
	samplename=$(echo "$reads" | awk -F'/' '{print $NF}')

	echo "Working on prefix $name and read files $samplename"
	
	###
	#extract reads mapped to a feature
	###
	
	mapping_out=$workingdir/${samplename}
	extract_out=$workingdir/extract_reads/${samplename}
	
	python3 extract_sequence_by_feature_gff.py -t "Rsp_SAT" "${mapping_out}.mapped.sorted.bam" "$REFERENCEGFF" "${extract_out}_Rsp_trimmed.fq"
	python3 extract_sequence_by_feature_gff.py -t "bp_SAT" "${mapping_out}.mapped.sorted.bam" "$REFERENCEGFF" "${extract_out}_1.688_trimmed.fq"	
	
	gzip -r $workingdir/extract_reads/

		
	for target in "${targets[@]}"
	do 		
	
		###
		#use mapped bam file to calculate pingpong score
		###
		
		#calculate pingpong score	
		extract_out=$workingdir/extract_reads/${samplename}_${target}
		output=${workingdir}/pingpong_analysis/${samplename}_${target}
							
		#go to cal_pingpong_scripts folder to calculate pingpong score
		cd "$cal_pingpong_dir" || ( echo "can't open the $cal_pingpong_dir"; exit )
	
		#convert fastq to insert file	
		./piPipes_fastq_to_insert "${extract_out}_trimmed.fq.gz" "${output}_trimmed.insert"

		#mapping
		bowtie -r -p 6 -v 2 -a --best --strata --sam $REFERENCEINDEX "${output}_trimmed.insert" | samtools view -@ 6 -b -h -F 4 - | samtools sort -@ 6 - | bedtools bamtobed -i - > "${output}.mapped.sorted.bed"
	
		#convert bed to bed2 format
		./piPipes_insertBed_to_bed2 "${output}_trimmed.insert" "${output}.mapped.sorted.bed" >  "${output}.mapped.sorted.bed2"
	
		#calculate pingpong score
		./CalandDraw_10th_ping_pong.sh "${output}.mapped.sorted.bed2"

	done
	
done < "$FASTQLIST"

echo "DONE!"