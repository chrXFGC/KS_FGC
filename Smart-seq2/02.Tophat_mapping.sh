#!/bin/bash
samp_name=$1
brief_name=$2  #with modify the sample names 
cln_dir=/clean_data/RNA
tophat_dir=/mapping/RNA/hg19
genome=/Data_Base/smrt-seq2/hg19/hg19_ERCC92_RGC # Genome reference
gtf_file=/Data_Base/smrt-seq2/hg19/refGene.gtf # GTF file 

# Mapping RNA-seq data with Tophat
tophat   -p 8 -G $gtf_file --library-type fr-unstranded                   \
   -o $tophat_dir/$brief_name  $genome                              \
   $cln_dir/$samp_name/$samp_name.R1.clean.fq.gz

