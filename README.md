# ReciMap

Recimap (reciprocal mapping tool) was developed as a bioinformatics command-line tool/pipeline to find rearrangements breakpoints 
between two closely related genomes. It uses Burrows-Wheeler read mapping (with Bowtie2) to map synthetic 
reads from one genome onto the other, and vice versa (reciprocally). The reads are created from the genomes and 
partially mapping reads with high MAPQ scores are used to identify the borders of rearrangement events
(breakpoints). 

Recimap was created by Casper Schutte as part of an M.Sc. project at the University of Stellenbosch, South Africa. 

## Usage:
Install the required libraries and software tools from the (ENVIRONMENT FILE). 
The two genomes need to be in FASTA format and the files need to be in the same directory as all the scripts. 
Running the pipeline is as simple as running the following command in the terminal:

```
./recimap.sh
```

