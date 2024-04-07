#!/bin/bash

sample=$1
dir=/output/${sample}
##--------------------
java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /GenomeAnalysisTK.jar -T BaseRecalibrator -R /genome.fa -I $dir/${sample}_realigned_reads.bam -knownSites \
dbsnp.vcf -knownSites indel.vcf  -o $dir/${sample}_recal_data.table
#
java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /GenomeAnalysisTK.jar -T BaseRecalibrator -R /genome.fa -I $dir/${sample}_realigned_reads.bam -knownSites \
dbsnp.vcf -knownSites indel.vcf -BQSR $dir/${sample}_recal_data.table -o $dir/${sample}_post_recal_data.table

java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /GenomeAnalysisTK.jar -T PrintReads -R /genome.fa -I $dir/${sample}_realigned_reads.bam -BQSR $dir/${sample}_recal_data.table -o $dir/${sample}_recal_reads.bam

