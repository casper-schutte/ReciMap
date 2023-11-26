#!/bin/env bash

# shellcheck disable=SC2086
echo
echo "------------------------------------------------------------------------- "
echo " ReciMap "
echo "-------------------------------------------------------------------------"
echo " Select two genomes in FASTA (.fa) format and ReciMap will identify the borders of "
echo " rearrangement events between them. ReciMap will also attempt to match up the synteny blocks "
echo " between the genomes. "
echo " ReciMap requires that the chromosomes in the FASTA files are named in the following format:"
echo " >chr1, >chr2, >chr3, etc. "
echo " ReciMap will not work if the chromosomes are named in any other way. "
echo " This script takes no arguements, you will be promted to select the two genomes you wish to "
echo " compare by selecting numbers from a list in the terminal. "
echo " If you made a mistake during selection, please press 'Ctrl + c' to exit the program."
echo " The borders of rearrangement events will be written to files called genA.txt and genB.txt"
echo " The final output showing the synteny blocks between the genomes will be written to text "
echo " files with the names blocks_(name_of_FASTA_file).txt "
echo " For more information, please see the README in repository for this project at: "
echo " https://github.com/casper-schutte/ReciMap "
echo "-------------------------------------------------------------------------"
echo
echo "Please select the first genome (Genome A)"
select file in $(ls *.f*)
do
	genomeA="$file"
	python3 make_reads_fa.py "$genomeA" readsA.fa
	echo "First genome: $genomeA"
	readsA="readsA.fa"
	echo
	break
done

echo "Please select the second genome (Genome B)"
select file in $(ls *.f*)
do
	genomeB="$file"
	echo "Second genome: $genomeB"
	python3 make_reads_fa.py "$genomeB" readsB.fa
	readsB="readsB.fa"
	echo
	break
done

bowtie2-build "$genomeA" "$genomeA"
bowtie2-build "$genomeB" "$genomeB"


single_report_align(){
	bowtie2 --local -f -p 8 -x $1 -U $2 -S $3.sam
}

# Map reads from genome B onto genome A
single_report_align "${genomeA}" ${readsB} SBA
# Map reads from genome A onto genome B
single_report_align "${genomeB}" ${readsA} SAB

multiple_report_align(){
	bowtie2 --local -f -k 2 -p 8 -x $1 -U $2 -S $3.sam
}

# Map reads from genome B onto genome A. Report the 2 best mapping positions.
multiple_report_align "${genomeA}" ${readsB} MBA
# Map reads from genome A onto genome B. Report the 2 best mapping positions.
multiple_report_align "${genomeB}" ${readsA} MAB

sort_sams(){
	samtools sort "$1" -O sam > $2
}
# Sort SAM files
sort_sams SBA.sam SBA.sorted.sam
sort_sams SAB.sam SAB.sorted.sam
sort_sams MBA.sam MBA.sorted.sam
sort_sams MAB.sam MAB.sorted.sam

# Record border positions in s**.txt for single-report mode, m** for multi-report mode
python3 find_borders_bash.py SBA.sorted.sam sBA.txt S $genomeA
python3 find_borders_bash.py SAB.sorted.sam sAB.txt S $genomeB
python3 find_borders_bash.py MBA.sorted.sam mBA.txt M $genomeA
python3 find_borders_bash.py MAB.sorted.sam mAB.txt M $genomeB

echo

python3 analyse_borders.py
rm *.bt2 # Remove temp files created by Bowtie2

# Identify synteny blocks between the genomes.
python3 get_synteny_blocks.py "$genomeA" B mAB.txt "$genomeB"
echo
python3 get_synteny_blocks.py "$genomeB" A mBA.txt "$genomeA"

# Clean up ugly temp files
rm mAB.txt mBA.txt sAB.txt sBA.txt


echo "Done..."

# Usage notes: