import sys

filepath = sys.argv[1]
read_file_name = sys.argv[2]


def get_chrom_names(genome):
    chrom_list = []
    with open(genome, "r") as a:
        for line in a.readlines():
            if line.startswith(">"):
                chrom_list.append(line.replace(">", "").replace("r", "").replace("\n", ""))
    a.close()
    # print(chrom_list)
    return chrom_list


def get_seq(genome):
    """
    Takes a FASTA file and returns the full sequence sans the title and newline characters
    input: FASTA (.fna) file
    output: string
    """

    with open(genome, "r") as a:
        b = (a.read()).split(">")
        z = []
        for x in b[1:]:
            z.append(x[x.index("\n") + 1:])
        d = []
        chrom_num = 0
        for x in z:
            d.append(x.replace("\n", ""))
            print(f"{chrom_num + 1}")
            print(len(x.replace("\n", "")))
            chrom_num += 1
    # print(len(d))
    # This method returns a list of the sequences in the FASTA file, where each chromosome is an element. This is
    # important because if this wasn't done, reads would be created that span 2 chromosomes (not biologically accurate).
    return d


def get_reads(seq_file):
    """
    Creates reads according to the parameters in the argument
    input: The full sequence as a string
    output: list of reads
    """
    read_len = 200
    reads = []
    read_pos = 0
    overlap = 190
    while read_pos < len(seq_file):
        reads.append(seq_file[read_pos:read_pos + read_len])
        read_pos += read_len - overlap
    return reads


def reads_to_file(reads, pos):
    """
    Writes a fastA file with the list of reads
    input: a list of reads
    output: a fastA file containing the reads and their number
    """
    my_file = open(read_file_name, "w")
    # reads = sorted(reads)
    # print(reads)
    for y in reads:
        for x in y:
            # my_file.write(">r" + str(pos) + "\n")
            my_file.write(">r" + chroms[chrom_num - 1] + "r" + str(pos) + "\n")
            my_file.write(str(x) + "\n")
            pos += (len(x) - 190)
    my_file.close()
    print(f"done writing reads to file: {read_file_name}")


if __name__ == "__main__":
    """
    This script creates reads from one of the genomes. The reads are 200bp in length and they are created in order of 
    their appearance in the genome. They are created at positions that start 10bp apart, such that the coverage is 
    approx 20x. The names are in the form:
    r + chromosome name + r + starting position relative to chromosome.
    Eg. rCh4r311
    NOTE that this script erases any instances of the letter "r" in the chromosome name, this is done because I haven't 
    found a better character to split by. 
    The lines that print are for debugging. 
    """
    chroms = get_chrom_names(filepath)
    # print(f"chroms: {chroms}")
    h = get_seq(filepath)
    # print(f"len(h) = {len(h)}")
    j = []
    chrom_num = 0
    for i in h:
        j.append(get_reads(i))
        # chrom_num += 1

    my_file = open(read_file_name, "w")
    # reads = sorted(reads)
    # print(reads)
    pos = 1
    prev_chrom = []
    for y in j:
        # print(f"chrom: {chroms[chrom_num]}")

        for x in y:
            if prev_chrom == [chroms[chrom_num]]:
                # my_file.write(">r" + str(pos) + "\n")
                my_file.write(">r" + chroms[chrom_num] + "r" + str(pos) + "\n")
                my_file.write(str(x) + "\n")
                pos += 10
                # prev_chrom = [chroms[chrom_num]]
            else:
                pos = 1
                prev_chrom = [chroms[chrom_num]]
        chrom_num += 1
    my_file.close()
    print(f"done writing reads to file: {read_file_name}")
