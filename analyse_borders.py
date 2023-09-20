"""
This script parses the 4 text files created with each run of the find_borders_bash.py script, the valid borders and
associated information is written to a text file output.txt.
"""

sAB = "sAB.txt"
sBA = "sBA.txt"
mAB = "mAB.txt"
mBA = "mBA.txt"
to_print = []


def extract_borders(file_name):
    """
    This function extracts the borders from the text files created by the find_borders_bash.py script. It also filters
    out the borders that are too short to be considered rearrangement breakpoints. Output format:
    [[chromosome, position, read_name, mapq, CIGAR string], ...]
    :param file_name:
    :return: checked_borders
    """
    with open(file_name, "r") as f:
        checked_borders = []
        borders = []
        my_lines = f.read()
        border_groups = my_lines.split(">")
        border_groups = border_groups[1:]
        border_groups = [x.split("\n") for x in border_groups]
        for i in border_groups:
            temp_list = []
            for j in i:
                if j != "" and j.isdigit() is False:
                    j = j.split("\t")
                    temp_list.append([
                        j[0], j[1], j[2].strip("r"), j[3], j[4], j[5]])
            borders.append(temp_list)
        for k in borders:
            threshold = 5
            # This threshold correlates to the "length" of the border. Since the reads are created 10bp apart, this
            # threshold x 10 (+- read_len) is the minimum "length" of a border.
            # The threshold in "find_borders.py" is more useful for distinguishing between SNPs/small
            # errors in alignment and actual rearrangement borders.
            start_pos = []
            read_names = []
            for border in k:
                start_pos.append((int(border[2].split("r")[1])))
                read_names.append(border[2])
            if max(start_pos) - min(start_pos) > threshold:
                checked_borders.append(k)

    to_print.append(f"{file_name}: total number of borders: {my_lines.count('>')}\n"
                    f" refined borders: {len(checked_borders)}")
    return checked_borders


def get_connected_borders(file_name):
    """
    This function extracts the connected borders from the text files created by the find_borders_bash.py script.
    :param file_name:
    :return: shared_borders
    """
    with open(file_name, "r") as f:
        connected_borders = []
        temp_list = f.readlines()
        for i in temp_list:
            i = i.strip("\n").split("\t")
            connected_borders.append(i[0])
    return connected_borders


def compare_borders(borders_a, borders_b):
    """
    This function is currently not used.
    :param borders_a:
    :param borders_b:
    :return:
    """
    shared = []
    a = []
    b = []
    for i in borders_a:
        a.append([i[0][0], i[0][1], int(i[0][2].strip("r"))])
    for i in borders_b:
        b.append(i[0][:2])
    for i in a:
        if i in b:
            shared.append(i)
    return shared


def output_all_borders(s1, s2, m1, m2):
    """
    This function writes the borders to a text file. The borders are written in the following format:
    [chromosome, position, read_name, mapq, CIGAR string]\n
    :param s1:
    :param s2:
    :param m1:
    :param m2:
    :return:
    """
    gen_a = []
    gen_b = []
    with open("genA.txt", "w") as gen_a_file, open("genB.txt", "w") as gen_b_file:
        for i in s1:
            temp = [i[0][0], i[0][1], i[0][2], i[0][5]]
            gen_a.append(temp)
        for i in s2:
            temp = [i[0][0], i[0][1], i[0][2], i[0][5]]
            gen_b.append(temp)
        to_print.append(f"\nGenome A:")
        for i in gen_a:
            to_print.append(i)
        to_print.append(f"connected borders: {len(get_connected_borders(mAB))}")
        for i in m1:
            to_print.append(i)
        to_print.append(f"\nGenome B:")
        for i in gen_b:
            to_print.append(i)
        to_print.append(f"connected borders: {len(get_connected_borders(mBA))}")
        for j in m2:
            to_print.append(j)
        for k in gen_a:
            gen_a_file.write(f"{k}\n")
        for l in gen_b:
            gen_b_file.write(f"{l}\n")


if __name__ == "__main__":
    sab = extract_borders(sAB)
    sba = extract_borders(sBA)
    mab = get_connected_borders(mAB)
    mba = get_connected_borders(mBA)
    output_all_borders(sab, sba, mab, mba)
    my_file = open("output.txt", "w")
    for item in to_print:
        my_file.write(f"{item}\n")
    print(f"Finished writing output to file: output.txt")
