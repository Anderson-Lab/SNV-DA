#!/bin/bash

sample=$1
snv_loci=$2
coverage=$3
mkdir coverages
bedtools genomecov -ibam $sample -bg | awk '$4 > ${$coverage-1}' | intersectBed -a stdin -b $snv_loci -wb > coverages/${sample::-4}.snv_loci_cov.bed
