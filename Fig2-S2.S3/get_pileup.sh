#!/bin/bash

"""
This script extract reads from alignment files, blast to repeat consensus sequence to get pileup along the consensus.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

module load blast/2.6.0+ 
module load anaconda3/5.2.0

#define your working directory
workingdir=""
cd $workingdir || ( echo "can't open the $workingdir"; exit )
echo "Working in dir: $workingdir"

#this is your list of read files with full paths to reads
#only put the prefix, not full read file names
FASTQLIST=$workingdir/FASTQLIST.txt


REFERENCEINDEX="dmel_scaffold2_plus0310.fasta"
REFERENCEGFF="dmel_scaffold2_plus0310_rm_for_htseq.gff3"
db="dmel_sat_canonical_blastdb"

mkdir $workingdir/extract_reads
mkdir $workingdir/pile_up_consensus


while FASTQLIST='' read -r reads || [[ -n "$reads" ]]; do
	
	name=$(echo "$reads" | awk -F'/' '{print $(NF-1)}')
	samplename=$(echo "$reads" | awk -F'/' '{print $NF}') 
	dir_path=$(echo "$reads" | awk -F'/' '{OFS="/"; $NF=""; print $0}')

	echo "Working on prefix $name and read files $samplename"
	
	output=$workingdir/${samplename}
	
	N=$(samtools view "${output}.mapped.sorted.bam" | wc -l)
	echo "got mapped using bowtie. total: $N."


###
#extract reads mapped to a feature (Rsp as an example)
###
	extract_out=$workingdir/extract_reads/${samplename}
	
	python3 extract_sequence_by_feature_gff.py -t "Rsp_SAT" "${output}.mapped.sorted.bam" "$REFERENCEGFF" "${extract_out}_Rsp_trimmed.fq"

	gzip -r $workingdir/extract_reads/
	
###
#blast reads to consensus dimer/trimer, get reads pile up along monomer (Rsp as an example)
###
	
	extract_out=$workingdir/extract_reads/${samplename}
	pile_out=$workingdir/pile_up_consensus/${samplename}	
	
	#extract reads for 1U/10A bias
	zcat "${extract_out}_Rsp_trimmed.fq.gz" | grep -A1 "^[>@]" | sed 's/@/>/g'| grep -v "^--" | awk 'NR%2==0 {OFS="\n"; split($1,a,""); if (a[10]=="A" || a[1]=="T") print b,$0} {b=$0}' > "${pile_out}_Rsp.fasta"

	#get number for all smallRNA and smallRNA with 1U/10A 
	n=$(( $(zcat "${extract_out}_Rsp_trimmed.fq.gz" | grep -A1 "^[>@]" | sed 's/@/>/g'| grep -v "^--" | wc -l) /2 ))
	p=$(( $(less "${pile_out}_Rsp.fasta" | wc -l) /2 ))
	perc=$((100*$p/$n))
	echo "Rsp: small RNA: $n; piRNA with 1T/10A: $p (${perc}%)"	

	#blast
	blastn -word_size 10 -query "${pile_out}_Rsp.fasta" -db "${db}/dmel_Rsp_trimer" -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge | awk '$4==$8' > "${pile_out}_Rsp_vs_consensus.blast" 
	#get pile up
	python3 get_pileup_along_monomer.py -l 233 -i "${pile_out}_Rsp_vs_consensus.blast" -o "${pile_out}_Rsp_vs_monomer_perc.distrib"
	
done < "$FASTQLIST"
echo "DONE!"
