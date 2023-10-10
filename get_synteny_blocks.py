import sys
import csv
import ast

fasta_name = sys.argv[1]
out_name = sys.argv[2]
connected_border_file = sys.argv[3]
other_genome = sys.argv[4]


def get_chrom_info(filename):
    """
    This function gets the lengths of each chromosome and stores them in a dictionary in the form:
    {chromosome name: length}
    This allows the blocks to be created correctly.
    :param filename:
    :return:
    """
    chromosome_lengths = {}
    with open(filename, "r") as my_fasta:
        current_chrom = ''
        current_seq_len = 0
        for current_line in my_fasta:
            if current_line.startswith('>'):
                # This indicates the start of a chromosome
                if current_chrom:
                    chromosome_lengths[current_chrom] = current_seq_len
                    current_seq_len = 0
                current_chrom = current_line.strip()[1:]
            else:
                current_seq_len += len(current_line.strip())

        if current_chrom:
            chromosome_lengths[current_chrom] = current_seq_len
    return chromosome_lengths


chrom_lens = get_chrom_info(fasta_name)
other_len = get_chrom_info(other_genome)

# Define a function to find synteny blocks given the list of borders
def get_synteny_blocks(border_file, lens):
    """
    This function creates the synteny blocks from the borders. It does this by first creating the chromosomes as blocks
    and then adding the borders to the appropriate chromosomes. The function also extracts the read information from
    the file containing the read information from the other genome. This information is used to correlate the blocks
    between the genomes.
    :param border_file:
    :return: syn_blocks
    """
    my_read_info = []
    if out_name == "A":
        my_read_info = get_read_info("sBA.txt")
    elif out_name == "B":
        my_read_info = get_read_info("sAB.txt")
    syn_blocks = []
    with open(border_file, "r") as file:
        lines = file.readlines()
        my_lines = []
        for i in lines:
            my_lines.append(eval(i))
        borders = []
        for j in enumerate(my_lines):
            borders.append((j[1][0], int(j[1][1]), f"r{j[1][2]}"))
            # This appends the chromosome name, border position, and read name to the borders list as a tuple.
        chromosomes = []
        chromosomes_with_borders = []
        for chrom in lens:
            for chrom_name in borders:
                if chrom == chrom_name[1]:
                    chromosomes_with_borders.append(chrom)
                    # If there reported borders in a chromosome, append the chromosome to the appropriate list.
        for chrom in lens:
            chromosomes.append((
                chrom, 1, lens.get(chrom), f"r0"
            ))
            # Create the full chromosomes as blocks, if there are no borders found then it gets appended as a block.
        # t is the list of temporary blocks
        t = []
        # The code above creates the chromosomes as single blocks. The code below iterates over them and
        # adds the borders
        for chrom in chromosomes:
            temp_list = []
            for b in borders:
                if chrom[0] == b[0]:
                    temp_list.append(b)
            temp_temp_list = [chrom]
            for x in temp_list:
                temp_temp_list.append(x)
            t.append(temp_temp_list)
        # The chromosomes are now in individual lists starting with an element that represents the entire chromosome.

        bwr = {}
        # Blocks with reads
        if my_read_info:
            for b in t:
                if len(b) < 2:
                    pass
                else:
                    for b2 in b:
                        to_append = []
                        # for each border in list b
                        if len(b2) == 3:
                            # if b2 is not representative of an entire chromosome (and therefore has borders)
                            for reads in my_read_info:
                                # reads is a list of reads mapping to the chromosome
                                for r2 in reads:
                                    # for each individual read in the list
                                    if r2[0] == b2[0] and int(b2[1]) - 5 <= int(r2[1]) <= int(b2[1]) + 5:
                                        to_append.append([r2[2][0], r2[2][1], r2[3]])

                            # Now the format of to_append is a list where each element is a list of 3 elements
                            # extracted from the read names: [chromosome name, position, +/-]. This
                            # information relates to the origin of the read, the position it is mapped to is already
                            # saved in b2. The problem is that to_append can contain multiple groups of reads eg:
                            # [['Ch3', '811', '-'], ['Ch3', '41', '+']]. They need to be extracted and saved to the
                            # dictionary as:
                            # {b2 : [chromosome name, min(pos), max(pos), +/-] }.

                            clean_list = []
                            current = []
                            pos_range = []
                            for info in to_append:
                                if not current:
                                    current.append(info)
                                    pos_range.append(int(info[1]))
                                elif info[2] == current[0][2]:
                                    current.append(info)
                                    pos_range.append(int(info[1]))
                                else:
                                    clean_list.append(
                                        [current[0][0], min(pos_range), max(pos_range), current[0][2]]
                                    )
                                    current = [info]
                                    pos_range = [int(info[1])]
                            # print(current)
                            if len(current) > 0:
                                clean_list.append([current[0][0], min(pos_range), max(pos_range), current[0][2]])
                            bwr[b2] = clean_list

        # print("Blocks:")
        for x in range(len(t)):
            # This iterates over the list where each element is a list comprised of the chromosome block and any
            # border blocks.
            rb = []  # real blocks
            if len(t[x]) == 1:
                rb.append(t[x][0])
                # If the list for a specific chromosome has length of 1, it means no borders are present and the
                # whole chromosome gets added as a block to the list of real blocks.
            else:
                # Here we iterate over the elements in the sub-lists. So this loop is repeated for each chromosome
                # list.
                curr_pos = 0
                counter = 1
                for y in range(len(t[x]) - 1):
                    if len(t[x][y]) == 4:
                        # a length of 4 indicates the start of a chromosome in the form ('Chr1', 1, 14700, 'r0')
                        # Chromosome, start, end, read name. The first two elements are appended as normal, followed
                        # by the start of the NEXT border. The "- 1" makes sure the blocks do not overlap.
                        if bwr.get(t[x][y]) is not None:
                            rb.append(
                                (
                                    t[x][y][0], t[x][y][1], int(t[x][y + 1][1]) - 1, bwr.get(t[x][y])
                                )
                            )
                        else:
                            rb.append((t[x][y][0], t[x][y][1], int(t[x][y + 1][1]) - 1, t[x][y][3]))
                        counter += 1
                        curr_pos = t[x][y + 1][1]
                        # The position is adjusted such that the next block will start in the correct place.

                    elif len(t[x][y]) == 3 and counter != y:
                        # If this does not represent the start OR end of a chromosome:
                        # append the same stuff as above (the read name is in a different position now).
                        rb.append(
                            (
                                t[x][y][0], t[x][y][1], int(t[x][y + 1][1]) - 1, bwr.get(t[x][y])
                            )
                        )
                        curr_pos = t[x][y + 1][1]

                    elif t[x][y][2] != curr_pos:
                        rb.append(
                            (
                                t[x][y][0], t[x][y][1], int(t[x][y + 1][1]), bwr.get(t[x][y])
                            )
                        )
                # Then append the last block that goes from the end of the last border and the end of the chromosome.
                rb.append(
                    (
                        t[x][y][0], curr_pos, lens.get(t[x][y][0]), bwr.get(t[x][y + 1])
                        # Suspect this is where a mistake is creeping in wrt to chromosome lengths.
                    )
                )

            for block in rb:
                # print(block)
                if block[2] <= block[1]:
                    pass
                    # Skip "blocks" that have the same start and end position, as these are not blocks.
                else:
                    syn_blocks.append(block)

        # print(syn_blocks)
        return syn_blocks


def get_read_info(file_with_reads):
    """
    This function extracts the information regarding the original positions of the reads that were partially mapped
    and thus indicate rearrangement borders. This information is later used to correlate the synteny blocks between
    the genomes.
    :param file_with_reads:
    :return:
    """
    file_name = file_with_reads
    lists_between_gt = []
    current_list = []
    with open(file_name, "r") as file:
        for my_line in file:
            my_line = my_line.strip()
            if my_line.startswith(">"):
                if current_list:
                    lists_between_gt.append(current_list)
                    current_list = []
            else:
                columns = my_line.split("\t")
                current_list.append(

                    [columns[0], columns[1], columns[2].split("r")[1:], columns[5]]
                )
    # Example of format:
    # [Chromosome name, position it maps to, [read name info from original genome: Chrom, start position of read], +/-]
    # The current form of the read data in lists_between_gt is: [1, 2, [3a, 3b], 4]
    # 1: Chromosome name,
    # 2: position read mapped to,
    # 3: [Chromosome read originated from, position of start of read in original chromosome]
    # 4: +/- signifying which side of the read mapped correctly to the other genome
    lists_between_gt.append(current_list)
    return lists_between_gt


class Block:
    """
    This class is not used, but may provide future functionality.
    """
    def __init__(self, num, block_info, read_info, block_range):
        self.num = num
        self.block_info = block_info
        self.read_info = read_info
        self.block_range = block_range


def correlate_blocks(genome):
    """
    This function correlates the synteny blocks in the genome with the read information from the other genome. Allowing
    the blocks to be labelled correctly. This function also creates a dictionary that contains the information for each
    block in the genome. The format of the dictionary is:
    {block info: [read info, read info, ...]}
    :param genome:
    :return:
    """
    gen_a = genome
    gen_b = []

    # Remember the 2 rules:
    # 1) Read info refers to the border at the START of the block (x[1] in other words).
    # 2) "+" means the read belongs to the END of the PREVIOUS block. "-" means it belongs
    #       at the current block.

    block_dict = {}
    prev_key = None
    all_keys = set()

    # This looks complicated, but I am just extracting fields I care about. This is semi vestigial because I used
    # to write a lot of the output to temp files just to have to parse the files, so I cut out the middle man and
    # this is the result.
    for sublist in gen_a[0:]:
        key = tuple(sublist[:3])
        sublist_4th = sublist[3]
        all_keys.add(key)

        if isinstance(sublist_4th, list):
            for subsublist in sublist_4th:
                if len(subsublist) >= 4:
                    if prev_key is not None and subsublist[3] == '+':
                        if prev_key in block_dict:
                            block_dict[prev_key].append(subsublist)
                        else:
                            block_dict[prev_key] = [subsublist]
                    else:
                        block_dict.setdefault(key, []).append(subsublist)
        if key not in block_dict.keys():
            block_dict[key] = []
        prev_key = key

    keys_to_remove = []
    # Remove the keys for the dict if the difference between the start and end blocks is less than 10.
    for b in block_dict:
        if abs(b[1] - b[2]) <= 10:
            keys_to_remove.append(b)
    for key in keys_to_remove:
        block_dict.pop(key)

    my_blocks = []
    counter = 1
    for b in block_dict:
        block = Block(counter, b, block_dict.get(b), (b[1], b[2]))
        my_blocks.append(block)
        counter += 1

    # for x in my_blocks:
    #     print(x.num, x.block_info, x.read_info, x.block_range)
    return my_blocks, block_dict


def get_connected_borders(filename):
    """
    This function extracts the 2 best mapping positions of reads from the alignment that was done using Bowtie2's
    multi-report mode (mAB.txt or mBA.txt). This information allows the function assign_numbers_to_n_blocks() to
    label the unlabelled synteny blocks.
    :param filename:
    :return:
    """
    with open(filename, "r") as f:
        connected_borders = []
        temp_list = f.read()
        temp_list = temp_list.split("\n")
        for x in temp_list:
            if len(x.split("\t")) > 1:
                x = x.split("\t")
                y = ast.literal_eval(x[1])
                for a in y:
                    if a[3] != a[7]:
                        if a[3] == "+":
                            border = [[a[0], a[1], a[3]], [a[4], a[5], a[7]]]
                            if border not in connected_borders:
                                connected_borders.append(border)
                        else:
                            border = [[a[4], a[5], a[7]], [a[0], a[1], a[3]]]
                            if border not in connected_borders:
                                connected_borders.append(border)
    return connected_borders


def assign_numbers_to_n_blocks(my_borders, connected_borders):
    """
    This function updates the unlabelled ("n") synteny blocks in the list my_borders by determining the correct number
    of the block from the connected borders.
    :param my_borders:
    :param connected_borders:
    :return:
    """

    for (chr1, pos1, sign1), (chr2, pos2, sign2) in connected_borders:
        right = None
        left = None
        for i in range(len(my_borders)):
            # Match the left side of cbs with right side (end) of blocks.
            if chr1 == my_borders[i][1][0] and abs(my_borders[i][1][2] - pos1) <= 10:
                left = my_borders[i][0]

            # Then match the right side of the cbs with the left side (start) of the blocks.
            if chr2 == my_borders[i][1][0] and abs(my_borders[i][1][1] - pos2) <= 10:
                right = my_borders[i][0]

            if right is not None and left is not None:
                if right == "n" or left == "n":
                    # print(f"{left, right}")
                    if left == "n" and right != "n":
                        left = right - 1
                        my_borders[i][0] = left
                    if right == "n" and left != "n":
                        right = left + 1
                        my_borders[i][0] = right


if __name__ == "__main__":
    """
    This is the main function that runs the code. It takes the synteny blocks from the two genomes and correlates them
    using the read information from the other genome. It then uses the connected borders to label the unlabelled blocks.
    The output is saved to a file in the form "blocks_{FASTA_NAME}.txt", the data format is:
    [block number, block info]
    """
    syntenyA = []
    syntenyB = []


    # for line in get_synteny_blocks(f"genA.txt", chrom_lens):
    #     syntenyA.append([line[0], line[1], line[2], line[3]])
    # for line in get_synteny_blocks(f"genB.txt", other_len):
    #     syntenyB.append([line[0], line[1], line[2], line[3]])

    # print(syntenyA)
    if out_name == "A":
        for line in get_synteny_blocks(f"genA.txt", chrom_lens):
            syntenyA.append([line[0], line[1], line[2], line[3]])
        for line in get_synteny_blocks(f"genB.txt", other_len):
            syntenyB.append([line[0], line[1], line[2], line[3]])
        genome1, g1_list = correlate_blocks(syntenyA)
        genome2, g2_list = correlate_blocks(syntenyB)
        for line in get_synteny_blocks(f"genA.txt", chrom_lens):
            syntenyA.append([line[0], line[1], line[2], line[3]])
        for line in get_synteny_blocks("genB.txt", other_len):
            syntenyB.append([line[0], line[1], line[2], line[3]])
    elif out_name == "B":
        for line in get_synteny_blocks(f"genA.txt", other_len):
            syntenyA.append([line[0], line[1], line[2], line[3]])
        for line in get_synteny_blocks(f"genB.txt", chrom_lens):
            syntenyB.append([line[0], line[1], line[2], line[3]])
        genome1, g1_list = correlate_blocks(syntenyB)
        genome2, g2_list = correlate_blocks(syntenyA)
        for line in get_synteny_blocks(f"genB.txt", chrom_lens):
            syntenyA.append([line[0], line[1], line[2], line[3]])
        for line in get_synteny_blocks("genA.txt", other_len):
            syntenyB.append([line[0], line[1], line[2], line[3]])
    connected_borders = get_connected_borders(connected_border_file)
    # The code above makes sure that the correct file is chosen from which to extract the block information.

    g1_ordered = []
    g2_ordered = []
    done = []
    # The code below correlates the blocks in the selected genome with the blocks in the "other" genome by taking the
    # read information regarding where the read originated. This can be done thanks to the way the reads are named
    # when they are created.
    for g1 in genome1:
        n1 = g1.num
        r1 = g1.block_range
        b1 = g1.block_info
        i1 = g1.read_info
        g1_ordered.append(
            (n1, b1)
        )
        for g2 in genome2:
            n2 = g2.num
            r2 = g2.block_range
            b2 = g2.block_info
            i2 = g2.read_info
            for read in i2:
                if b1[0].replace("r", "") == read[0]:

                    if read[3] == "-":
                        if r1[0] in range(read[1], read[2] + 150):
                            # Useful for debugging
                            # print(b1)
                            # print(read)
                            # print(n2)
                            if b2 not in done and (n2, b1) not in g2_ordered:
                                g2_ordered.append((
                                    n2, b1
                                ))
                            done.append(b1)
                    elif read[3] == "+":
                        if r1[1] in range(read[1], read[2] + 150):
                            # print(b1)
                            # print(read)
                            # print(n2)
                            if b2 not in done and (n2, b1) not in g2_ordered:
                                g2_ordered.append((
                                    n2, b1
                                ))
                            done.append(b1)
    print("__________________________________________")
    my_borders = []
    position = 1
    for block_info in g1_list:
        found = False
        for number, info in g2_ordered:
            if info == block_info:
                found = True
                print((number, info))
                position = info[2]
                if [number, info] not in my_borders:
                    my_borders.append([number, info])
                break

        if not found:
            position = block_info[2]
            if ["n", block_info] not in my_borders:
                print(("n", block_info))
                my_borders.append(["n", block_info])

    assign_numbers_to_n_blocks(my_borders, connected_borders)
    # for num, block in my_borders:
    #     print([num, block])
    # Save output to files in the form "blocks_{FASTA_NAME}.txt"
    with open(f"blocks_{sys.argv[4].split('.')[0]}.txt", "w") as final_file:
        final_file.write(f"Original order of blocks in genome {sys.argv[1]}\n")
        # Write the original order of the blocks in the OTHER genome to the file.
        for x in genome2:
            final_file.write(f"{x.num} - {x.block_info}\n")
        final_file.write(f"\n")
        final_file.write(f"Order of blocks in genome {sys.argv[4]}\n")
        for num, block in my_borders:
            final_file.write(f"{num} - {block}\n")

# Go through the logic. I need to write about in the paper anyway. What reads are being mapped to what genome to
# create the files and stuff. There is also a mistake with the lengths of the 5th chromosome... it wasn't there
# before I changed stuff today.
# Will make a good figure too, can create a huge flow diagram showing the different data
