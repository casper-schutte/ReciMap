sAB = "sAB.txt"
sBA = "sBA.txt"
mAB = "mAB.txt"
mBA = "mBA.txt"
to_print = []


def extract_borders(file_name):
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
            threshold = 7
            # This threshold correlates to the "length" of the border. Since the reads are created 10bp apart, this
            # threshold x 10 (+- read_len) is the minimum "length" of a border.
            # The threshold in "find_borders.py" is more useful for distinguishing between SNPs/small
            # errors in alignment and actual rearrangement borders.
            start_pos = []
            for border in k:
                start_pos.append(int(border[2]))
            if max(start_pos) - min(start_pos) > threshold:
                checked_borders.append(k)
    to_print.append(f"{file_name}: total number of borders: {my_lines.count('>')}\n"
                    f" refined borders: {len(checked_borders)}")
    return checked_borders


def get_connected_borders(file_name):
    with open(file_name, "r") as f:
        connected_borders = []
        temp_list = f.readlines()
        for i in temp_list:
            i = i.strip("\n").split("\t")
            connected_borders.append(i[1][:])
    return connected_borders


def compare_borders(borders_a, borders_b):
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
    gen_a = []
    gen_b = []
    for i in s1:
        temp = [i[0][0], i[0][1], i[0][5]]
        gen_a.append(temp)
    for i in s2:
        temp = [i[0][0], i[0][1], i[0][5]]
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
