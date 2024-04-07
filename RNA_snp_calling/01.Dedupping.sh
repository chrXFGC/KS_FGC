#!/bin/bash

sample=$1
dir=/output/${sample}

java -Djava.io.tmpdir=`pwd`/tmp -Xmx60g -jar /picard.jar MarkDuplicates \
I=$dir/${sample}.bam  \
O=$dir/${sample}_marked.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=$dir/${sample}.metrics

