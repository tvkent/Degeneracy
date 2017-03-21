#!/bin/bash
#pipeline for getting 4fold degenerate sites
#from gff and fasta, using python, awk, and bedtools
#Tyler Kent
#14 March 2017

###################################
# SET UP PATHS
#
# This is the only portion of the
# pipeline that needs to be
# adjusted. Assume gff is gzipped.
###################################

gff='/ohta/tyler.kent/Storage/Crubella_183_v1.0.gene_exons.gff3.gz'
fasta='/ohta/tyler.kent/Storage/Capsella_rubella_v1.0_combined.fasta'
CDSbedout='/ohta/tyler.kent/Storage/Crub_CDS.bed'
fastaCDSout='/ohta/tyler.kent/Storage/Crub_CDS.tab'
fourfoldbedout='/ohta/tyler.kent/Storage/Crub_4fold'

###################################
# STEP 1: GET BED FILE OF CDS AND
# SHIFT TO MATCH PHASE
#
# CDS sequence in GFF format
# contains sections of translated
# sequence, with phase info, which
# indicates the start of the first
# codon.
###################################

bash gff2bed.sh <(zcat ${gff}) CDS | awk -f gffphaseshift.awk - > ${CDSbedout}

###################################
# STEP 2: USE BED FILE AND FASTA
# TO GET FILE OF POS AND SEQUENCE
#
# Use bedtools to get relevant 
# fasta sequence into useable
# format.
###################################

bedtools getfasta -s -tab -fi ${fasta} -bed ${CDSbedout} > ${fastaCDSout}

###################################
# STEP 3: CONVERT FASTA DNA
# SEQUENCE INTO CODONS, FLIP FOR
# PHASE, AND REPORT 4FOLD SITES
#
# Using python and a codon table
###################################

python degeneracy.py -i ${fastaCDSout} -o ${fourfoldbedout}.bed 

###################################
# STEP 4: SORT OUTPUT
#
# Just need to sort the python 
# output like you would a normal
# bed file, and drop mistake dups.
###################################

sort -k1,1 -k2,2n ${fourfoldbedout}.bed | uniq > ${fourfoldbedout}.sorted.bed
