from class_SAM import SAM
from itertools import groupby
import sys

file = SAM()
filepath = sys.argv[1]
outname = sys.argv[2]
number = file.ReadSAMFile(filepath)
mode = sys.argv[3]

"""
It is VITAL that the variable "read_len" is set to be the same as the actual length of the reads
"""
read_len = 200


def get_border_reads():
    """
    This function creates a list with various fields for the reads. The reads which do not map properly (indicated
    by the CIGAR string) are returned.
    :return:
    """
    pos = []
    mapq = []
    chrom = []
    read_name = []
    cigar = []
    rocigar = []
    for i in range(number):
        pos.append(int(file.GetField(i, "POS")))
        mapq.append(int(file.GetField(i, "MAPQ")))
        chrom.append(file.GetField(i, "RNAME"))
        read_name.append(file.GetField(i, "QNAME"))
        cigar.append(file.GetCIGARvalues(i))
        rocigar.append(file.GetField(i, "CIGAR"))
    borders = []
    # Now the reads that had CIGAR strings that did not show full matches are appended to the "borders" list, along
    # with their associated information.
    for read in enumerate(rocigar):
        if read[1] == f"{read_len}M" or read[1].count("S") == 0 or read[1].count("M") > 1:
            pass
        else:
            borders.append([chrom[read[0]], pos[read[0]], read_name[read[0]], mapq[read[0]], read[1]])
    return borders


def get_borders(borders):
    """
    This function takes the reads which are suspected to contain rearrangement breakpoints and extracts the CIGAR
    strings. Then, the reads are returned BUT with ADJUSTED positions. These positions refer to the last (or first)
    mapping base pair, depending on whether it matched initially or started matching somewhere along its length.
    :param borders:
    :return breakpoints:
    """
    bps = []
    for border in borders:
        cigar_str = border[4]
        temp_list = []
        if cigar_str.count("M") == 1:
            for character in cigar_str:
                if character == "M":
                    bps.append(temp_list)
                    break
                    # Honestly a pretty a cool piece of original code.
                elif character.isalpha() and character != "M":
                    temp_list = ["-"]
                    bps.append(temp_list)
                    break
                elif character.isdigit():
                    if not temp_list:
                        temp_list.append(character)
                    else:
                        temp_list.append("".join(character))
    matching = []
    for numbers in bps:
        matching.append("".join(numbers[0:]))
    breakpoints = []
    for m in enumerate(matching):
        if borders[m[0]][4].count("M") == 1 and borders[m[0]][4].count("S") < 2:
            if m[1].count("-") == 0:
                breakpoints.append([
                    borders[m[0]][0], borders[m[0]][1] + int(m[1]), borders[m[0]][2], borders[m[0]][3],
                    borders[m[0]][4], "+"
                ])
                # The "+ int(m[1])" shifts the position to the breakpoint instead of where the read starts.
                # For the reads that do match at the start (below), the position already refers to the potential
                # breakpoints.
            elif m[1].count("-") > 0:
                breakpoints.append([
                    borders[m[0]][0], borders[m[0]][1], borders[m[0]][2], borders[m[0]][3],
                    borders[m[0]][4], "-"
                ])
    return breakpoints


def refine_breakpoints(bp):
    """
    This function takes the breakpoints and groups them by whether they are within a certain distance from one another.
    Then, the breakpoints that are within a certain distance from the non-confident breakpoints are added to the
    group. The non-confident breakpoints are defined as those with a MAPQ score of less than 30. The distance threshold
    is 10bp for the confident breakpoints and 5bp for the non-confident breakpoints. Output format:
    [[chromosome, position, read_name, mapq, CIGAR string], ...]
    :param bp:
    :return: border_areas
    """

    conf_borders = []
    non_conf_borders = []
    mapq_threshold = 30
    for read in bp:
        if int(read[3]) >= int(mapq_threshold):
            conf_borders.append(read)
        else:
            non_conf_borders.append(read)
    border_areas = []
    allowance = 10
    counter = 0
    # counter to keep track of where we are in the list
    conf_borders = sorted(conf_borders, key=lambda x: x[1])
    for border in conf_borders:
        # This loop groups the most confident borders by whether they are within a certain distance from one another.
        if not border_areas:
            border_areas.append([border])
        elif border[1] - allowance <= border_areas[counter][0][1] <= border[1] + allowance \
                and border[0] == border_areas[counter][0][0] and border[5] == border_areas[counter][0][5]:
            border_areas[counter].append(border)
        else:
            border_areas.append([border])
            counter += 1
    for border in border_areas:
        if border_areas.count(border) > 1:
            print(f"{border} has multiple hits")
    done = []
    non_conf_allowance = 5
    for border in non_conf_borders:
        if border not in done:
            for cb in border_areas:
                if border[1] - non_conf_allowance <= cb[0][1] <= border[1] + non_conf_allowance and \
                        border[0] == cb[0][0] and border[5] == cb[0][5]:
                    border_areas[border_areas.index(cb)].append(border)
                    done.append(border)
    for border in done:
        if border_areas.count(border) > 1:
            print(f"{border} has multiple hits")
    return border_areas


def write_to_file(breakpoints, duplicated_reads):
    """
    This function writes either the breakpoints or the connected borders to a file, depending on the mode (M or S).
    If the mode is "S", the potential borders that had less than 5 reads supporting them are discarded. If the mode is
    "M", the duplicated reads are written to the file. The format for the output is:
    [chromosome, position, read_name, mapq, CIGAR string, orientation]
    :param breakpoints:
    :param duplicated_reads:
    :return None:
    """
    if mode not in ["M", "S"]:
        print("Invalid mode. Please provide either 'M' or 'S'.")
        return

    try:
        my_file = open(file_name, "w")
    except IOError as e:
        print(f"Error opening the file: {e}")
        return

    if mode == "M":
        if duplicated_reads is not None:
            for dup in sorted(duplicated_reads):
                my_file.write(f"{dup}\t{duplicated_reads.get(dup)}\n")
    elif mode == "S":
        breakpoints.sort(key=lambda x: x[0])
        for area in breakpoints:
            area = [i for i in area if len(i) >= 6 and i[2].count("r") >= 2]
            if len(area) <= 1:  # Ignore areas with less than two valid elements
                continue

            area.sort(key=lambda x: x[3], reverse=True)
            positions = [int(i[2].split("r")[2]) for i in area]

            if not positions:  # Skip if positions list is empty
                continue

            # Calculate the difference and threshold length
            diff = max(positions) - min(positions)
            threshold_len = 1  # You can change this threshold to your desired value

            if len(area) > threshold_len:
                my_file.write(f">{diff}")
                my_file.write("\n")
                for bp in area:
                    my_file.write(f"{bp[0]}\t{bp[1]}\t{bp[2]}\t{bp[3]}\t{bp[4]}\t{bp[5]}\n")

    my_file.close()


def find_duplicate_reads(breakpoints):
    """
    This function finds the reads that map to more than one place. The reads are then grouped by whether they map to
    the same place. The format for the output is:
    TODO: update this
    :param breakpoints:
    :return: border_dict
    """

    flat_borders = [border for x in breakpoints for border in x]
    read_names = [x[2] for x in flat_borders]
    duplicates = []
    for name in enumerate(read_names):
        if read_names.count(name[1]) > 1:
            duplicates.append(flat_borders[name[0]])
            # flat_borders[name[0]] format: ['Chr3', 1401, 'rCh3r641', 255, '50M150S', '+']
    duplicates.sort(key=lambda x: x[2])
    # The sorcery below checks for reads mapping to more than one place
    grouped_borders = [list(x) for y, x in groupby(duplicates, lambda x: x[2])]
    # grouped borders format example:
    # [[['Chr1', 701, 'rCh1r1231', 255, '160S40M', '-'], ['Chr1', 1401, 'rCh1r1231', 38, '160M40S', '+']],
    # [['Chr1', 701, 'rCh1r1241', 255, '150S50M', '-'], ['Chr1', 1401, 'rCh1r1241', 25, '150M50S', '+']],...]
    connected_borders = []
    connected_reads = []

    for border in grouped_borders:
        if border[0][0] != border[1][0] or border[0][1] != border[1][1]:
            # if either the positions or the chromosome names are different:
            connected_borders.append([border[0][0], border[0][1],
                                      border[1][0], border[1][1]])
            connected_reads.append([border[0][0], border[0][1], border[0][2], border[0][5],
                                    border[1][0], border[1][1], border[1][2], border[1][5]])
    border_dict = {}
    for x in connected_borders:
        if connected_borders.count(x) > 5:
            temp_list = []
            for y in connected_reads:
                if [y[i] for i in [0, 1, 4, 5]] == x:
                    temp_list.append(y)
            # print(temp_list)
            border_dict[(x[0], x[1], x[2], x[3])] = temp_list


    if not grouped_borders:
        return None
    else:
        return border_dict


if __name__ == "__main__":
    read_info = get_border_reads()
    find_breakpoints = get_borders(read_info)
    file_name = outname  # Desired name for output file
    refined_breakpoints = refine_breakpoints(find_breakpoints)
    dup_reads = find_duplicate_reads(refined_breakpoints)
    write_to_file(refined_breakpoints, dup_reads)
