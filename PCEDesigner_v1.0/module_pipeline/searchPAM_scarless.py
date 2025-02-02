###########################################################################
# Program: searchPAM_scaless.py
# Function: search the crispr PAM seq and position and each position also
# Authtor: Xiaoli Shi developed by modifying searchPAM from Xianrong Xie  2016.12.1
# Date: 2024.12
############################################################################
from __future__ import division
import re
import seqDeal


unbi_dict = {'N' : '[ATCGRYMKSWHBVDN]',\
             'R' : '[RAG]', \
             'Y' : '[YCT]', \
             'M' : '[MAC]', \
             'K' : '[KGT]', \
             'S' : '[SGC]', \
             'W' : '[WAT]', \
             'H' : '[HWMYATC]', \
             'B' : '[BKSYGTC]', \
             'V' : '[VRSMGAC]', \
             'D' : '[DRKWGAT]' \
            }

class targetPAM:
    def __init__(self, seq_str, pamStr, pam_direct, guide_len, target_strand, start_pos):
        self.seq_str = seq_str
        self.pamStr = pamStr 
        self.pam_direct = pam_direct
        self.guide_len = guide_len
        self.target_strand = target_strand
        self.start_pos = start_pos

        # create dict to record the founded PAM postion
        self.rawPamInfo_dict = {}

        self.pamPosIndex_lst = [] 
        self.pam_lst = []
        self.direct_lst = []
        self.targetPos_lst = []
        self.targetSeq_lst = []
        self.targetGC_lst = []
        self.seqID_lst = []

        # run the PAM position get 
        self.pamPos()

    def pamPos(self):
        singleSeq = self.seq_str
        if self.target_strand == "+":
            pam = self.pamStr
        else:
            pam = seqDeal.revSeq(self.pamStr)
        self.pamStr_re = ''
        for i in pam:
            if i in list(unbi_dict.keys()):
                repStr = unbi_dict[i]
            else:
                repStr = i
            self.pamStr_re += repStr

        """
        # the PAM to be compile
        if self.pamStr == "NGG":
            self.pamStr1 = "[ATCGRYMKSWHBVDN]GG" 
            self.pamStr2 = "CC[ATCGRYMKSWHBVDN]"
            self.directType = "upper" # based on the pamStr, juduge the direction of the target site, NGG's target is 5' upper, 20nt
        if self.pamStr == "TTN":
            self.pamStr1 = "TT[ATCGRYMKSWHBVDN]"
            self.pamStr2 = "[ATCGRYMKSWHBVDN]AA"
            self.directType = "down1" # Cpf1's target is 3' down, 24 nt
        if self.pamStr == "TTTN":
            self.pamStr1 = "TTT[ATCGRYMKSWHBVDN]"
            self.pamStr2 = "[ATCGRYMKSWHBVDN]AAA"
            self.directType = "down2"
        """
        pamCompile = re.compile(self.pamStr_re)
        # cut the seq to the pam with the length of base
        seqID = 0 
        self.seqID_dict = {}
        for st in range(len(singleSeq)-len(self.pamStr)): 
            shortSeq = singleSeq[st:st+len(self.pamStr)]
            pamFind = pamCompile.match(shortSeq)
            if pamFind:
                # getPam = pamFind.group()
                # get the target position, here is 1-based
                if self.target_strand == "+":
                    if self.pam_direct == "three":
                        targetStart = st-int(self.guide_len)
                        targetEnd = st + len(self.pamStr) 
                        if targetStart < 0:
                            continue
                    elif self.pam_direct == "five":
                        targetStart = st
                        targetEnd = st + int(self.guide_len) + len(self.pamStr) 
                        if targetEnd > len(singleSeq):
                            continue

                    targetSeq = singleSeq[targetStart:targetEnd] 
                    targetPos = self.start_pos+targetStart + 1
                    targetDirect = "+"
                else:
                    if self.pam_direct == "three":
                        targetStart = st
                        targetEnd = targetStart+int(self.guide_len) + len(self.pamStr) 
                        if targetEnd > len(singleSeq):
                            continue
                    elif self.pam_direct == "five":
                        targetStart = st-int(self.guide_len)
                        targetEnd = st + len(self.pamStr) 
                        if targetStart < 0:
                            continue

                    targetSeq = singleSeq[targetStart:targetEnd] 
                    targetSeq = seqDeal.revSeq(targetSeq)   # rev the target
                    targetPos = self.start_pos + targetStart + 1
                    targetDirect = "-" 
                if self.pam_direct == "three":
                    perGC = "%0.1f" % float((targetSeq[:-len(self.pamStr)].count("G")+targetSeq[:-len(self.pamStr)].count("C"))/len(targetSeq[:-len(self.pamStr)])*100)
                else:
                    perGC = "%0.1f" % float((targetSeq[len(self.pamStr):].count("G")+targetSeq[len(self.pamStr):].count("C"))/len(targetSeq[len(self.pamStr):])*100)

                info_lst = [targetDirect, targetSeq, self.pamStr, perGC]
                self.rawPamInfo_dict.setdefault(targetPos, []).extend(info_lst)
