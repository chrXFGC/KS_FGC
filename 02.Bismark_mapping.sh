#!/bin/bash

sample=$1
cln_fq1=$2
cln_fq2=$3
ref=/hg19/Bismark2

/Bismark-0.22.3/bismark   --multicore 4   --score_min L,0,-0.6  -quiet --non_directional -o /outdir/$sample  $ref  $cln_fq1
/Bismark-0.22.3/bismark   --multicore 4   --score_min L,0,-0.6  -quiet --non_directional -o /outdir/$sample  $ref  $cln_fq2


