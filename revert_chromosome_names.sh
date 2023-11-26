#!/bin/bash

# Input genome files
genomeA="$1"
genomeB="$2"

# Create temporary files for reverting changes
reverted_genomeA=$(mktemp)
reverted_genomeB=$(mktemp)

# Extract the original chromosome names from genome A
grep -o ">.*" "${genomeA}" | tr -d ">" > genomeA_original_chr_names.txt

# Extract the original chromosome names from genome B
grep -o ">.*" "${genomeB}" | tr -d ">" > genomeB_original_chr_names.txt

# Function to revert chromosome names
revert_chromosome_names() {
    input_genome="$1"
    output_genome="$2"
    while read -r original_name; do
        read -r new_name
        sed "s/$new_name/$original_name/g" "$input_genome" >> "$output_genome"
    done < genomeA_original_chr_names.txt
}

# Revert chromosome names in genome A
revert_chromosome_names "${genomeA}" "${reverted_genomeA}"

# Revert chromosome names in genome B
revert_chromosome_names "${genomeB}" "${reverted_genomeB}"

# Rename the reverted genomes to their original filenames
mv "${reverted_genomeA}" "${genomeA}"
mv "${reverted_genomeB}" "${genomeB}"

# Clean up temporary files
#rm genomeA_original_chr_names.txt
#rm genomeB_original_chr_names.txt
