#!/bin/bash

"""
This script calculate log2FC from count files.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

module load anaconda3/5.2.0

workingdir=$(pwd)
echo "Working in dir: $workingdir"

study="study_design.txt"

###
#RDC/Moon normalized by miRNA analysis as example
###

suffix="out_norm_miRNA"
echo "working on files with suffix: $suffix"

#wildtype=("SRR1187947" "SRR5445299" "SRR5445307" "SRR891254" "SRR891255" "SRR5803096" "SRR1524976" "SRR037091" "SRR827714")
#Rhino=("SRR1187948" "SRR5445300" "SRR5803097" "SRR037092")
#Cutoff=("SRR1187949" "SRR891256" "SRR891257" "SRR1524977")
#Deadlock=("SRR1187950" "SRR827720")
#Moonshiner=("SRR5445301" "SRR5445318")

###############

name="Rhino_1"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR1187947.${suffix}" -t "SRR1187948.${suffix}" -o "${name}" 
echo "finish generating input for $name."

name="Rhino_2"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR5445299.${suffix}" "SRR5445307.${suffix}" -t "SRR5445300.${suffix}" -o "$name" 
echo "finish generating input for $name."

name="Rhino_3"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR5803096.${suffix}" -t "SRR5803097.${suffix}" -o "$name" 
echo "finish generating input for $name."

name="Rhino_4"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR037091.${suffix}" -t "SRR037092.${suffix}" -o "$name" 
echo "finish generating input for $name."

###############

name="Cutoff_1"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR1187947.${suffix}" -t "SRR1187949.${suffix}" -o "$name" 
echo "finish generating input for $name."

name="Cutoff_2"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR891254.${suffix}" "SRR891255.${suffix}" -t "SRR891256.${suffix}" "SRR891257.${suffix}" -o "$name" 
echo "finish generating input for $name."

name="Cutoff_3"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR1524976.${suffix}" -t "SRR1524977.${suffix}" -o "$name" 
echo "finish generating input for $name."


################

name="Deadlock_1"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR1187947.${suffix}" -t "SRR1187950.${suffix}" -o "$name" 
echo "finish generating input for $name."

name="Deadlock_2"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR827714.${suffix}" -t "SRR827720.${suffix}" -o "$name" 
echo "finish generating input for $name."


##################

name="Moonshiner"
study="study_design.txt"
python3 summarize_multi_count_for_heatplot_norm.py -n "$name" -s "$study" -c "SRR5445299.${suffix}" "SRR5445307.${suffix}" -t "SRR5445301.${suffix}" "SRR5445318.${suffix}" -o "$name" 
echo "finish generating input for $name."


##################


echo "finished analyzing all data, starting to summarize..."

output=$(echo "$workingdir" | awk -F'/' '{print $NF}')
echo "generating summary file: ${output}.log2FC.txt"

echo -n > "${output}.log2FC.txt"
for infile in *.summary; do
	awk -F'\t' 'NR==1 {next}; NR > 1 {OFS="\t"; split($3,a,"_"); print a[1],$3, $4, $9}' "$infile" >> "${output}.log2FC.txt"
done

echo "DONE!"

