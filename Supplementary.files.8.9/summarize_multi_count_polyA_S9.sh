#!/bin/bash

"""
This script run DESeq2 with read count files for Supllementary file 9 (poly-A RNA).
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

module load anaconda3/5.2.0
module load r/3.5.1

workingdir=$(pwd)
echo "Working in dir: $workingdir"

target="target.txt"
header="header.csv"
study="study_design.txt"

#output from count_proportional python script
suffix="out"

#combine Rsp_L and Rsp_R in one
#for ifile in ${workingdir}/*.${suffix}
for ifile in ${workingdir}/SRR*.${suffix}
do
        echo "working on file: ${ifile}..."
        rm "${ifile}_Rsp"
        grep "Rsp" $ifile | awk 'BEGIN{OFS="\t";sum1=0;sum2=0} {sum1+=$2;sum2+=$3} END{print "Rsp_SAT",sum1,sum2}' > "tmp"
        grep "bp" $ifile | awk 'BEGIN{OFS="\t";sum1=0;sum2=0} {sum1+=$2;sum2+=$3} END{print "1pt688_SAT",sum1,sum2}' >> "tmp"
        cat $ifile "tmp" | grep -v "bp" | grep -v "\-Rsp" > "${ifile}_sat"
        echo "generated file: ${ifile}_sat"
done

suffix="out_sat"

echo "working on files with suffix: $suffix"

name="Rhino"
echo "Running DESeq2 for $name..."
python3 /scratch/alarracu_lab/Xiaolu/pipelines/htseq-count/summarize_multi_count_for_DESeq2.py -n "$name" -s "$study" -c "SRR8578630.${suffix}" "SRR8578631.${suffix}" "SRR8578632.${suffix}" "SRR8578633.${suffix}" -t "SRR8578634.${suffix}" "SRR8578635.${suffix}" -o "$name" 
Rscript /scratch/alarracu_lab/Xiaolu/pipelines/htseq-count/DESeq2.R "${name}.countData" "${name}.colData" "${name}.csv"


#summarize
echo "finished analyzing all data, starting to summarize..."

output=$(echo "$workingdir" | awk -F'/' '{print $NF}')
echo "generating summary file: ${output}.csv.summary"

#print header first 
cat "$header" | awk -F',' '{OFS=","} {print "\"protein\"","\"repeatID\"", $2, $3, $4, $5, $6, $7}' > "${output}.csv.summary"
#extract info about target repeat and add one column about protein studied
for infile in *.csv; do
	echo "working on file: $infile"
	name=$(echo "$infile" | awk -F'.' '{print $1}')
	echo "$name"
	grep -f "$target" "$infile" | awk -v var="$name" -F',' '{OFS=","} {print "\""var"\"", $1, $2, $3, $4, $5, $6, $7}' >> "${output}.csv.summary"
done


echo "DONE!"


