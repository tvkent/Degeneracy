#!/bin/bash
#take gff, choose annotation type, output bed file
#only argument is annotation type as a string

awk -v var=$2 'BEGIN {OFS="\t"}{if ($3==var) print $1,($4-1),$5,$7,$8}' $1 |\
sort -k1,1 -k2,2n
