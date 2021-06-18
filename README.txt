The following repository contains code for data analysis and figures in:

Heterochromatin-dependent transcription of satellite DNAs in the Drosophila female germline. BioRxiv https://doi.org/10.1101/2020.08.26.268920
Xiaolu Wei, Danna G. Eickbush, Iain Speece, Amanda M. Larracuente
University of Rochester

In this paper, we show that complex satDNAs (>100-bp repeat units) are transcribed into long noncoding RNAs and processed into piRNAs (PIWI interacting RNAs). This satDNA piRNA production depends on the Rhino-Deadlock-Cutoff complex and the transcription factor Moonshiner—a previously-described non-canonical pathway that licenses heterochromatin-dependent transcription of dual-strand piRNA clusters. We show that this pathway is important for establishing heterochromatin at satDNAs. Therefore, satDNAs are regulated by piRNAs originating from their own genomic loci.  This novel mechanism of satDNA regulation provides insight into the role of piRNA pathways in heterochromatin formation and genome stability.

Please direct questions to Xiaolu Wei at xiaolu_wei@urmc.rochester.edu or Amanda Larracuente at alarracu@bio.rochester.edu, or code author in each script.

We organize the data files and scripts in this repository according to the figures or supplementary files in Wei et al. 2020. 

We also include scripts and general data files that are not associated with specific figures, including:
#count.sh: bash script used to count and extract reads that mapped to repeats from alignment files.
#extract_sequence_by_feature_gff.py: python script used to extract reads that map to a repeat feature of interest and output fastq format.
#htseq_bam_count_proportional.py: python script used to count reads mapped to repeat features from alignment files.
#dmel_scaffold2_plus0310_rm_for_htseq.gff3: repeat annotation file for our heterochromatin enriched assembly from Chang and Larracuente 2019 (https://doi.org/10.1534/genetics.118.301765).
#different.1pt688.anaylsis: We performed analysis of the 1.688 satDNA repeats for ChIP-seq, RIP-seq and RNA-seq data in different ways to make sure that different analysis methods do not affect our conclusions. The results for the different analysis methods are organized by data type in this folder.

Usage: Please see scripts for usage instructions.


