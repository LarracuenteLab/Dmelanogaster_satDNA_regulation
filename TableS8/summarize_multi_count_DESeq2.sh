#!/bin/bash

"""
This script run DESeq2 with count files to get differential expression results.
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
echo "working on files with suffix: $suffix"

name="Rhino"
echo "Running DESeq2 for $name..."
python3 summarize_multi_count_for_DESeq2.py -n "$name" -s "$study" -c "SRR1187952.${suffix}" "SRR1187953.${suffix}" "SRR5445322.${suffix}" "SRR5445323.${suffix}" -t "SRR1187954.${suffix}" "SRR5445324.${suffix}" "SRR5445325.${suffix}" -o "$name" 
Rscript DESeq2.R "${name}.countData" "${name}.colData" "${name}.csv"

name="Cutoff"
echo "Running DESeq2 for $name..."
python3 summarize_multi_count_for_DESeq2.py -n "$name" -s "$study" -c "SRR1187952.${suffix}" "SRR1187953.${suffix}" -t "SRR1187955.${suffix}" -o "$name" 
Rscript DESeq2.R "${name}.countData" "${name}.colData" "${name}.csv"

name="Deadlock"
echo "Running DESeq2 for $name..."
python3 summarize_multi_count_for_DESeq2.py -n "$name" -s "$study" -c "SRR1187952.${suffix}" "SRR1187953.${suffix}" -t "SRR1187956.${suffix}" -o "$name" 
Rscript DESeq2.R "${name}.countData" "${name}.colData" "${name}.csv"

name="Moonshiner"
echo "Running DESeq2 for $name..."
python3 summarize_multi_count_for_DESeq2.py -n "$name" -s "$study" -c "SRR5445322.${suffix}" "SRR5445323.${suffix}" -t "SRR5445326.${suffix}" "SRR5445327.${suffix}" -o "$name" 
Rscript DESeq2.R "${name}.countData" "${name}.colData" "${name}.csv"


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
