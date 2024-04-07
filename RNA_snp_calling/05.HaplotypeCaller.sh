#!/bin/bash

sample=$1
dir=/output/${sample}

java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /GenomeAnalysisTK.jar -T HaplotypeCaller \
        -R /genome.fa \
        -I $dir/${sample}_recal_reads.bam \
        -dontUseSoftClippedBases \
        -stand_call_conf 20.0 \
        -o  $dir/${sample}.raw_snps.indels.vcf 
