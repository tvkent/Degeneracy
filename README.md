# Degeneracy
degeneracy.py will take a bedtools getfasta tab output file and return 0 and 4fold sites.
Annotated coding sequence that overlaps other CDS annotations will be filtered to keep only sites that are in agreement.
CDS that overlaps in position but exists on opposite strands will be filtered such that 4fold sites will only be called if they don't overlap with 0fold sites on the opposite strand and if they overlap between strands (basically a double check).

get_4fold_sites.sh will run the full pipeline starting with a fasta and gff annotation.

#Dependencies
bedtools, python3
