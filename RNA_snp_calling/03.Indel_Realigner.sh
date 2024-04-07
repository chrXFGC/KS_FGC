#!/bin/bash

sample=$1
dir=/output/${sample}

java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /GenomeAnalysisTK.jar -T RealignerTargetCreator -R /genome.fa -I $dir/${sample}_sort_split_marked.bam -known \
indel.vcf    -o $dir/${sample}_realignment_targets.list

java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /GenomeAnalysisTK.jar -T IndelRealigner -R /genome.fa -I $dir/${sample}_sort_split_marked.bam -targetIntervals \
$dir/${sample}_realignment_targets.list -known indel.vcf  -o $dir/${sample}_realigned_reads.bam

