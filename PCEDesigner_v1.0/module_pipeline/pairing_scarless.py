#########################################
# Author: Xiaoli Shi
# 2024.12.04

import pandas as pd
import re
from collections import defaultdict

class SequencePairFinder:
    def __init__(self, pamInfo_dict, offDisplay_dict, pair_dic):
        self.pamInfo_dict = pamInfo_dict
        self.offDisplay_dict = offDisplay_dict
        self.pair_dic = pair_dic

        self.sorted_pamInfo = self.read_filter_ontarget_data()

    def extract_numeric(self, value):
        # using re to pick the numeric part from item
        match = re.search(r"([0-9]*\.?[0-9]+)", value)
        return float(match.group(0)) if match else None

    def read_filter_ontarget_data(self):
        nested_lst = []
        for tarID, tar_info_lst in self.pamInfo_dict.items():
            tarStrand, tarSeq, PAM, tarGC, tarPos = tar_info_lst
            nested_lst.append([tarID, tarPos, tarStrand, tarSeq, PAM, tarGC])
        # sort item based on their start position
        sorted_lst = sorted(nested_lst, key= lambda x:x[1])
        '''
        for item in sorted_lst:
            print(item)
        exit(0)
        '''
        return sorted_lst

    def find_ontarget_pairs(self):
        # initialize a list to deposit paired target sequence and OT sum
        pairs_on_dict = {}

        # iterate all items
        for i in range(len(self.sorted_pamInfo)):
            for j in range(i + 1, len(self.sorted_pamInfo)):
                tarID1 = self.sorted_pamInfo[i][0]
                tarID2 = self.sorted_pamInfo[j][0]

                seq1 = self.sorted_pamInfo[i][3]
                seq2 = self.sorted_pamInfo[j][3]

                pam1 = self.sorted_pamInfo[i][4]
                pam2 = self.sorted_pamInfo[j][4]

                gc1 = self.sorted_pamInfo[i][5]
                gc2 = self.sorted_pamInfo[j][5]

                strand1 = self.sorted_pamInfo[i][2]
                strand2 = self.sorted_pamInfo[j][2]

                start1 = self.sorted_pamInfo[i][1]
                end1 = start1 + 22

                start2 = self.sorted_pamInfo[j][1]
                end2 = start2 + 22
                '''
                print("tarID1: ", tarID1)
                print("offDisplay_dict")
                for k, v in self.offDisplay_dict.items():
                    print(k, v)
                '''
                off_score1 = self.offDisplay_dict[str(tarID1)]
                off_score2 = self.offDisplay_dict[str(tarID2)]


                cur_pair_info_dict = self.pair_dic[str(tarID1)]
                cur_pair_len_lst = sorted([len(i) for i in cur_pair_info_dict.keys()], reverse=True)
                if len(cur_pair_len_lst) != 0:
                    pair_len1 = max(cur_pair_len_lst)
                else:
                    pair_len1 = 0

                cur_pair_info_dict = self.pair_dic[str(tarID2)]
                cur_pair_len_lst = sorted([len(i) for i in cur_pair_info_dict.keys()], reverse=True)
                if len(cur_pair_len_lst) != 0:
                    pair_len2 = max(cur_pair_len_lst)
                else:
                    pair_len2 = 0

                # compute distance bw two sequences
                distance = start2 + 6 - (end1 -6) -1

                # paired or not
                if distance < 20 or distance > 70:
                    continue
                if strand1 == '+' and strand2 == '-':
                    # Sum of OT
                    # total_score = score1 + score2
                    # Append sequence pair and total_score into list
                    # pairs_loss.append((tarID1, seq1, str(start1), str(end1), strand1, str(score1), tarID2, seq2, str(start2), str(end2), strand2, str(score2)))
                    key = (tarID1, tarID2)
                    value = (seq1, start1, end1, strand1, pam1, gc1, off_score1, pair_len1, seq2, start2, end2, strand2, pam2, gc1, off_score2, pair_len2, distance)
                    pairs_on_dict[key] = value

        if len(pairs_on_dict) > 0:
            return pairs_on_dict
        else:
            return {}
    '''
    def write_output(self, pairs_on_dict,del_pairs_lst, output_file):
        # Sort the items based on the sum of OT
        sorted_pairs = sorted(pairs, key=lambda x: x[-1])
        with open(output_file, 'w') as f:
            f.write("seq1\tstart1\tend1\tstrand1\tOT_score1\tseq2\tstart2\tend2\tstrand2\tOT_score2\tsum_score\n")
            for pair in sorted_pairs:
                f.write("\t".join(pair) + "\n")
        with open(output_file, 'w') as f:
            f.write("seq1\tstart1\tend1\tstrand1\tmax_OT_score1\tseq2\tstart2\tend2\tstrand2\tmax_OT_score2\tsum_max_OT\toff_target_score\n")
            for key, value in pairs_on_dict.items():
                if key in del_pairs_lst:
                    continue
                f.write("\t".join(map(str, value))+'\n')
        print(f"Done, paired sequence has been deposited into {output_file} file.")
    '''

'''
if __name__ == "__main__":
    file_path = '/mnt/Data_diskE/Gene_Editing_Research_Server/LiHongChao_TKO_PE_gene_silence/XieXR-target_desgin/linux_tardesign_maize/result/results/'
    on_path = file_path+'Zm00001eb000010_1.result.txt'

    pattern = re.escape(file_path) + "(.*?)" + r'\.result\.txt'
    match = re.search(pattern, on_path)
    if match:
        extracted_string = match.group(1)
        output_file = file_path + extracted_string + ".TKOpair.txt"
        print(f"Output file set to: {output_file}")
    else:
        print("Pattern not found in on_path.")
        output_file = file_path+'TKOpair.txt'

    finder = SequencePairFinder(on_path, off_path)
    pairs_on_dict = finder.find_ontarget_pairs()
    finder.write_output(pairs_on_dict,pairs_off_dict, del_pairs_lst, output_file)
'''
