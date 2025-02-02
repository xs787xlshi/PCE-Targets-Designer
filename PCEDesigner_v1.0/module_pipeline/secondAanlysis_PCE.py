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
        # key: string; value: [(ref_start, ref_end),(tar_start, tar_end)]
        print(k, v)
    exit(0)
    '''
    return pair_dic

def outPair(pamInfo_dict, sgRNA_lst, PAM):
    out_dic = {}
    # reversed but not complemented
    sgRNA_cis = sgRNA_lst[0][::-1] 
    sgRNA_trans = sgRNA_lst[1][::-1]
    # for l in open(targetFile, "r"):
    for tarID, curInfo_lst in pamInfo_dict.items():
        '''
        ls = l.rstrip()
        if re.search("^>", ls):
            tarID = ls[1:] 
            continue
        '''
        strand = curInfo_lst[0]
        ls = curInfo_lst[1]
        if PAM == "NGG":
            tarSeq = ls[:-(len(PAM))]
        else:
            tarSeq = ls[len(PAM):]
        # complemented but not reversed
        revTarSeq = revSeq(tarSeq)
        # print ">%s\t%s" % (tarID, revTarSeq)
        if strand == "+":
            sgRNA = sgRNA_cis
        else:
            sgRNA = sgRNA_trans
        pair_dic = pairSeq(revTarSeq, sgRNA) 
        out_dic.setdefault(str(tarID), pair_dic)
        # print ">%s\t%s\t%s" % (tarID, revTarSeq, pair_dic)
    return out_dic


