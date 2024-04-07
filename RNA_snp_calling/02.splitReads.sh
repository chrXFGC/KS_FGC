#!/bin/bash

sample=$1
dir=/output/${sample}

java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /GenomeAnalysisTK.jar \
-T SplitNCigarReads \
-R /genome.fa \
-I $dir/${sample}_marked.bam \
-o $dir/${sample}_sort_split_marked.bam \
-rf ReassignOneMappingQuality \
-RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
