import re
import sys

def revSeq(sequence):
    base_info = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'N' : 'N',\
            'V' : 'B', 'D' : 'H', 'B' : 'V', 'H' : 'D', 'W' : 'W', 'S' : 'S',\
            'K' : 'M', 'M' : 'K', 'Y' : 'R', 'R' : 'Y'}
    seq_lst = []
    for base in sequence:
        seq_lst.append(base_info[base]) 
    # seq_lst = seq_lst[::-1]
    rev_seq = ''.join(seq_lst)
    return rev_seq

def pairSeq(tar, ref):
    pair_dic = {} 
    minPair = 8
    # print("revTarSeq: ", tar)
    # print("sgRNA scaffold: ", ref)
    for i in range(0, (len(tar)-minPair+1)):
        curSeq = tar[i:] 
        for k in range(0, len(curSeq)-minPair+1):
            childSeq = curSeq[0:len(curSeq)-k]
            # print("childSeq_orig: ", childSeq)
            childSeqPos = (i, i+len(childSeq))
            childSeq = re.sub("A", "[AG]", childSeq)
            childSeq = re.sub("C", "[CT]", childSeq)
            # print("childSeq_regu: ", childSeq)
            comChildSeq = re.compile(childSeq)
            mat = comChildSeq.finditer(ref)
            for j in mat:
                jStr = j.group()
                jPos = j.span()
                # print("jStr ", jStr)
                # print("jPos ", jPos)
                isExist = False
                for p in list(pair_dic.keys()):
                    if re.search(jStr, p):
                        isExist = True 
                    elif re.search(p, jStr):
                        del pair_dic[p]
                if not isExist:
                    pos_lst = [jPos, childSeqPos]
                    pair_dic.setdefault(jStr,[]).extend(pos_lst)
    '''
    for k, v in pair_dic.items():
        #key: string; value: [(ref_start, ref_end),(tar_start, tar_end)]
        print(k, v)
    exit(0)
    '''
    return pair_dic

def outPair(pamInfo_dict, sgRNA, PAM):
    out_dic = {}
    sgRNA = sgRNA[::-1] # reversed but not complemented
    # for l in open(targetFile, "r"):
    for tarID, curInfo_lst in pamInfo_dict.items():
        '''
        ls = l.rstrip()
        if re.search("^>", ls):
            tarID = ls[1:] 
            continue
        '''
        ls = curInfo_lst[1]
        # if PAM == "NGG": # cas9
        tarSeq = ls[:-(len(PAM))]
        '''
        else: # this is for cas12
            tarSeq = ls[len(PAM):]
        '''
        revTarSeq = revSeq(tarSeq) # complemented but not reversed
        # print ">%s\t%s" % (tarID, revTarSeq)
        pair_dic = pairSeq(revTarSeq, sgRNA) 
        out_dic.setdefault(str(tarID), pair_dic)
        # print ">%s\t%s\t%s" % (tarID, revTarSeq, pair_dic)   
    return out_dic

"""
(xx, targetFile)=sys.argv
sgRNA = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"

outPair(targetFile, sgRNA)
"""
