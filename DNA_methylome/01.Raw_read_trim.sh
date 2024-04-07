#!/bin/bash

fq1=$1
fq2=$2
cln_dir=$3

/trim_galore --path_to_cutadapt /cutadapt --clip_R1 9 --clip_R2 9 --paired --quality 20 --phred33 --stringency 3 --gzip --length 50  --trim1 --paired $fq1 $fq2 --output_dir $cln_dir 

