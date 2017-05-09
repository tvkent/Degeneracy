# Degeneracy
degeneracy.py will take a bedtools getfasta tab output file and return 0 and 4fold sites.
Annotated coding sequence that overlaps other CDS annotations will be filtered to keep only sites that are in agreement.
CDS that overlaps in position but exists on opposite strands will be filtered such that 4fold sites will only be called if they don't overlap with 0fold sites on the opposite strand and if they overlap between strands (basically a double check).

get_4fold_sites.sh will run the full pipeline starting with a fasta and gff annotation.

# Dependencies
bedtools, python3

# Running the pipeline

get_4fold_sites.sh has a config section with the relevant paths uncommented.
Scroll through the script if the filenames are unclear.

# Output

Output will be two files and few intermediates.
0fold sites will be written to ${output path}.0fold.bed, 4fold to ${output path}.4fold.bed

**Note**: the name column in the output bed only reflects *one* gene name in a group of overlapping annotations.
