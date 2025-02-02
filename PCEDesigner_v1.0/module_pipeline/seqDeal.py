#####################################################
# Program: seqDeal.py
# Function: to Deal the fasta sequence and outPut all base as ATCG etc.
# Author: Xianrong Xie
#####################################################
import re
import sys

def readFa(seqFile):
    seq_dict = {}
    seqName = "Seq"
    seq_lst = []
    for line in open(seqFile,"r"):
        line = line.strip()
        if line.startswith(">"):
            if len(seq_lst) != 0:
                genome_seq = ''.join(seq_lst)
                seq_dict[seqName] = genome_seq
            seqName = line[1:]
            seq_lst = []
            continue

        seqEachLine = line.upper()
        matchBase = re.findall("[ATCGRYMKSWHBVDN]", seqEachLine)
        seq_lst.append(''.join(matchBase))

    genome_seq = ''.join(seq_lst)
    seq_dict[seqName] = genome_seq

    return seq_dict

def readSeq(seqStr):
    seqStr = seqStr.upper()
    #matchBase = re.findall("[ATCGRYMKSWHBVDN]", seqStr)
    #seqFormat = ''.join(matchBase)
    return seqStr

def outSeq(fa):
    seq_lst = []
    for l in open(fa, "r"):
        ls = l.rstrip()
        if re.search(">", ls):
            continue
        matchBase = re.findall("[ATCGRYMKSWHBVDN]", ls)
        seq_lst.append(matchBase) 
    seq = "".join(seq_lst)
    return seq


def countSeqLen(seqIn):
    return len(seqIn)


def revSeq(refseq):
    base={}
    base['A']='T'
    base['T']='A'
    base['C']='G'
    base['G']='C'
    base['R']='Y'
    base['Y']='R'
    base['M']='K'
    base['K']='M'
    base['S']='S'
    base['W']='W'
    base['H']='D'
    base['B']='V'
    base['V']='B'
    base['D']='H'
    base['N']='N'
    base['a']='T'
    base['t']='A'
    base['c']='G'
    base['g']='C'
    base['r']='Y'
    base['y']='R'
    base['m']='K'
    base['k']='M'
    base['s']='S'
    base['w']='W'
    base['h']='D'
    base['b']='V'
    base['v']='B'
    base['d']='H'
    base['n']='N'

    rs=''
    for i in refseq:
        rs+=base[i]
    return rs[::-1]

def judge_bad(tarSeq, GC, pair_len_lst, off_val):
    warn = ""
    '''
    # search TTTT or longer ==> change TTTT to TTTTT or longer
    polyT_len = 0
    #polyT = re.search("T{4,}", tarSeq)
    polyT = re.search("T{5,}", tarSeq)##revised by Xiaoli Shi-20241203
    if polyT:
        polyT_len = len(polyT.group())
    ''' 
    # pair-sgRNA secondary structure
    if len(pair_len_lst) != 0:
        pair_len = max(pair_len_lst)
    else:
        pair_len = 0
    # GC 
    if polyT_len >= 5: 
        warn += "!"
    if float(GC) <= 20 or float(GC) >= 80:
        warn += "!"
    if pair_len >= 8:
        warn += "!"
    if float(off_val) >= 0.7:
        warn += "!"
    return warn

# test the program
"""
(xx, f) = sys.argv
seqDict = readFa(f)
print seqDict
"""
