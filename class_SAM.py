"""
This script is from Prof, I use it to get the fields (such as MAPQ) from the SAM files.
I did not write this.
"""


class SAM:

    def __init__(self):
        self.entries = []
        self.members = 0
        self.option_keys = {'AS': 0, 'XN': 1, 'XM': 2, 'XO': 3, 'XG': 4, 'NM': 5, 'MD': 6}
        self.fields = {'QNAME': 0, 'FLAG': 1, 'RNAME': 2, 'POS': 3, 'MAPQ': 4, 'CIGAR': 5, 'RNEXT': 6, 'PNEXT': 7,
                       'TLEN': 8, 'SEQ': 9, 'QUAL': 10, 'OPT': 11}
        self.cigar_categories = 'MIDNSHP=X'

    def __del__(self):
        pass

    def ReadSAMFile(self, filename):
        file_handle = open(filename, 'r')
        for line in file_handle:
            if (line[0] != '@'):
                temp_list = line.split('\t')
                pad = ' '
                temp_string = pad.join(temp_list[11:])
                del temp_list[11:]
                temp_list.append(temp_string.replace('\n', ''))
                self.entries.append(temp_list)
                self.members += 1
        return self.members

    def GetEntry(self, index):
        if (index >= 0 and index < self.members):
            return self.entries[index]
        else:
            return 0

    def GetField(self, index, field):
        if (self.GetEntry(index) != 0):
            if field in self.fields:
                return self.GetEntry(index)[self.fields[field]]
            else:
                return 0
        else:
            return 0

    def HasOPTField(self, index):
        if (self.GetEntry(index) != 0):
            if (len(self.GetEntry(index)) == 12):
                return True
            else:
                return False

    def GetOPTValue(self, index, tag):
        if (self.GetEntry(index) != 0):
            if (11 < len(self.GetEntry(index))):
                temp = self.GetEntry(index)[11]
                temp2 = temp.split(' ')
                temp3 = temp2[self.option_keys[tag]].split(':')
                return (temp3[2])
            else:
                return str('')

    def GetCIGARvalues(self, index):
        cigar_string = self.GetField(index, 'CIGAR')
        cigar_string_length = len(cigar_string)
        self.cigar_categories = 'MIDNSHP=X'
        cigar_dictionary = {}
        for category in self.cigar_categories:
            start = 0
            end = 0
            value = 0
            while (end < cigar_string_length):
                if (cigar_string[end].isalpha()):
                    if (cigar_string[end] == category):
                        value += int(cigar_string[start:end])
                    end += 1
                    start = end
                else:
                    end += 1
            cigar_dictionary[category] = value
        return cigar_dictionary

    def GetNumberOfNucleotidesInRead(self, index):
        cigar_dictionary = self.GetCIGARvalues(index)
        number_of_nucleotides = cigar_dictionary['M'] + cigar_dictionary['I'] + cigar_dictionary['D']
        return number_of_nucleotides

    def GetNumberOfMismatchesInRead(self, index):
        return self.GetOPTValue(self, index, self.options['NM'])

    def GetRolledOutCIGAR(self, index):
        cigar_string = self.GetField(index, 'CIGAR')
        length = self.GetNumberOfNucleotidesInRead(index)
        rolled_out_cigar = [0] * length

    def MapMDStringOntoReadSequence(self, index):
        rolled_out_cigar = self.GetRolledOutCIGAR(index)


"""
This script is from Prof, I use it to get the fields (such as MAPQ) from the SAM files.
I did not write this.
"""
