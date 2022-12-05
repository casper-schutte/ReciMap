import sys

# This python script successfully takes a FASTA file and creates reads according to set parameters
# these reads are written to a FASTA (.fa) file with the sequence name correlating numerically to the
# order of the reads.
filepath = sys.argv[1]
read_file_name = sys.argv[2]


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
        for x in z:
            d.append(x.replace("\n", ""))
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


def reads_to_file(reads):
    """
    Writes a fastA file with the list of reads
    input: a list of reads
    output: a fastA file containing the reads and their number
    """
    my_file = open(read_file_name, "w")
    reads = sorted(reads)
    read_num = 1
    for y in reads:
        for x in y:
            my_file.write(">r" + str(read_num) + "\n")
            my_file.write(str(x) + "\n")
            read_num += 1
    my_file.close()
    print(f"done writing reads to file: {read_file_name}")


if __name__ == "__main__":
    h = get_seq(filepath)
    j = []
    new_pos = 0
    for i in h:
        j.append(get_reads(i))
        new_pos += len(i)
    reads_to_file(j)
