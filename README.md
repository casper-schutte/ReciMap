# ReciMap

ReciMap (a reciprocal mapping tool) was developed as a bioinformatics command-line tool/pipeline to find rearrangements
breakpoints between two closely related genomes. It uses Burrows-Wheeler read mapping (with Bowtie2) to map synthetic 
reads from one genome onto the other, and vice versa (reciprocally). The reads are created from the genomes and 
partially mapping reads with high MAPQ scores are used to identify the borders of rearrangement events
(breakpoints). The functionality/scope of the pipeline was extended to include the ability to identify synteny blocks 
between the two genomes.

ReciMap was created by Casper Schutte as part of an M.Sc. project at the University of Stellenbosch, South Africa. 

## Usage:
This pipeline is best suited for use between very closely related genomes and was not designed to handle complex
(overlapping) rearrangement events.
Important Note: In its current state, the pipeline only works when ALL of the chromosome names in BOTH FASTA files are 
in the following format: (This is temporary, and a fix is in the works)
```
>Chr1
(sequence for Chr1)
>Chr2
(sequence for Chr2)
>Chr3
(sequence for Chr3)
etc
```
Install the required libraries and software tools from the (ENVIRONMENT FILE). 
The two genomes need to be in FASTA format and the files need to be in the same directory as all the scripts. 
Running the pipeline is as simple as running the following command in the terminal:

```
./recimap.sh
```

You will be prompted to select the first genome:
```
------------------------------------------------------------------------- 
 ReciMap 
-------------------------------------------------------------------------
 Select two genomes in FASTA (.fa) format and ReciMap will identify the borders of 
 rearrangement events between them. ReciMap will also attempt to match up the synteny blocks 
 between the genomes. 
 This script takes no arguements, you will be promted to select the two genomes you wish to 
 compare by selecting numbers from a list in the terminal. 
 If you made a mistake during selection, please press 'Ctrl + c' to exit the program.
 The borders of rearrangement events will be written to files called genA.txt and genB.txt
 The final output showing the synteny blocks between the genomes will be written to text 
 files in the form blocks_(name_of_FASTA_file).txt 
 For more information, please see the repository for this project at: 
 https://github.com/casper-schutte/ReciMap 
-------------------------------------------------------------------------

Please select the first genome (Genome A)
1) reference_genome.fa
2) rearranged_genome.fa
#? 
```
The order in which the genomes are selected makes no difference. The borders will be in the form of 
synteny blocks. The rearrangement borders are where one block ends and another begins. 
This output will be written to a text file called blocks_(name_of_FASTA_file).txt. One of 
these files will be created for each of the two genomes (see an example of the format under 
the "Format Example" heading below):

## Format Example:
The format is as follows
```
n - (chromosome name, start position of block, end position of block)
```
Where "n" is the numerical label of the block, representing the original order of the blocks in the other genome.
For example, the contents of the file "blocks_rearranged_genome.txt":
```
Original order of blocks in genome reference_genome.fa
1 - ('Chr1', 1, 15400)
2 - ('Chr2', 1, 210)
3 - ('Chr2', 211, 14140)
4 - ('Chr3', 1, 700)
5 - ('Chr3', 701, 1401)
6 - ('Chr3', 1401, 5880)
7 - ('Chr3', 5882, 6790)
8 - ('Chr3', 6791, 8540)
9 - ('Chr4', 1, 487)
10 - ('Chr4', 493, 4480)
11 - ('Chr5', 1, 2240)
12 - ('Chr5', 2241, 3361)
13 - ('Chr5', 3361, 6292)

Order of blocks in genome rearranged_genome.fa
n - ('Chr1', 1, 700)
n - ('Chr1', 701, 1400)
n - ('Chr1', 1401, 14700)
n - ('Chr2', 1, 210)
3 - ('Chr2', 211, 14490)
5 - ('Chr3', 1, 700)
4 - ('Chr3', 701, 1401)
6 - ('Chr3', 1401, 5880)
8 - ('Chr3', 5882, 7630)
7 - ('Chr3', 7631, 8540)
9 - ('Chr4', 1, 490)
12 - ('Chr4', 493, 1608)
10 - ('Chr4', 1611, 5600)
11 - ('Chr5', 1, 2240)
13 - ('Chr5', 2241, 5172)
```
The numerical labels of the blocks represent the order of the blocks in the OTHER genome. For example,
the block labelled "12" is located between blocks 9 and 10. This is due to a rearrangement event as in the reference 
genome it was the 12th block, located between blocks 11 and 13. 

The other file ("blocks_reference_genome.txt") will be in the same format, but with the order of the blocks in the 
rearranged genome shown as the "correct" (numerical) order, and the order in which those synteny blocks appear in the 
reference genome.

The blocks labelled "n" are blocks whose synteny could not be established. This is often due to duplication events, 
which do not leave a clear border. Block ('Chr2', 1, 210) in this example could not be identified as being syntenous 
between the genomes, this is due to the method used to identify synteny blocks and the fact that the borders used to 
identify the presence of a border resulted from an inversion of this block (see the "Method" section in my thesis for more
information on the inner workings of this pipeline).

For more information on the rearrangement border, an 
output file called "output.txt" will be left by the pipeline, giving more detail on each identified rearrangement border