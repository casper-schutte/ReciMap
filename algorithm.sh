#!/bin/env bash

# shellcheck disable=SC2086
echo
echo " Please make sure that the genome FASTA files you wish to use start with 'ref_'"
echo " and that the FASTA files containing the reads start with 'reads_'."
echo " If you made a mistake during selection, please press 'Ctrl + c' to exit the program."
echo
echo "Please select the first genome (Genome A)"
select file in $(ls *.fna)
do
	genomeA=$file;
	python3 make_reads_fa.py "$genomeA" readsA.fa
	echo "First genome: $genomeA"
	readsA=readsA.fa
	echo
break;
done

echo "Please select the second genome (Genome B)"
select file in $(ls *.fna)
do
	genomeB=$file;
	echo "Second genome: $genomeB"
	python3 make_reads_fa.py "$genomeB" readsB.fa
	readsB=readsB.fa
	echo
break;
done

bowtie2-build "$genomeA" "$genomeA"
bowtie2-build "$genomeB" "$genomeB"


single_report_align(){
	bowtie2 --local -f -p 8 -x $1 -U $2 -S $3.sam
}

single_report_align "${genomeA}" ${readsB} SBA
single_report_align "${genomeB}" ${readsA} SAB

multiple_report_align(){
	bowtie2 --local -f -k 2 -p 8 -x $1 -U $2 -S $3.sam
}

multiple_report_align "${genomeA}" ${readsB} MBA

multiple_report_align ${genomeB} ${readsA} MAB

sort_sams(){
	samtools sort "$1" -O sam > $2
}

sort_sams SBA.sam SBA.sorted.sam
sort_sams SAB.sam SAB.sorted.sam
sort_sams MBA.sam MBA.sorted.sam
sort_sams MAB.sam MAB.sorted.sam

python3 find_borders.py SBA.sorted.sam sBA.txt S
python3 find_borders.py SAB.sorted.sam sAB.txt S
python3 find_borders.py MBA.sorted.sam mBA.txt M
python3 find_borders.py MAB.sorted.sam mAB.txt M

echo

python3 analyse_borders.py
echo "Done..."
