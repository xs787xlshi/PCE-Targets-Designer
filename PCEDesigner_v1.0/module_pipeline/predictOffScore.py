from __future__ import division
import sys
import re


class offPredictor:
    def __init__(self, sgRNA, cigr):
        self.offTable = "../module_pipeline/offScore.table"
        self.baseRev_dic = { \
            'A':'T', \
            'T':'A', \
            'C':'G', \
            'G':'C', \
            'R':'Y', \
            'Y':'R', \
            'M':'K', \
            'K':'M', \
            'S':'S', \
            'W':'W', \
            'H':'D', \
            'B':'V', \
            'V':'B', \
            'D':'H', \
            'N':'N'}
        self.pamOffScore_dic = { \
            "A":0, \
            "T":0, \
            "C":0, \
            "G":1, \
            "AA":0, \
            "AC":0, \
            "AG":0.259259259, \
            "AT":0, \
            "CA":0, \
            "CC":0, \
            "CG":0.107142857, \
            "CT":0, \
            "GA":0.069444444, \
            "GC":0.022222222, \
            "GG":1, \
            "GT":0.016129032, \
            "TA":0, \
            "TC":0, \
            "TG":0.038961039, \
            "TT":0}

        self.sgRNA = sgRNA
        self.cigr = cigr
        self.readOffTable()
        self.scoreCal()
        self.scoreSum()

    def readOffTable(self):
        self.offScore_dict = {}
        for l in open(self.offTable, "r"):
            lst = l.rstrip().split("\t") 
            basePair = lst[0]
            pos = int(lst[1])
            score = float(lst[2])
            self.offScore_dict.setdefault(pos, {}).setdefault(basePair, score)

    def scoreCal(self):
        matchPos_lst = re.findall("\d+", self.cigr)
        un_matchBase_lst = re.findall("[A-Z]", self.cigr)
        # judge the unmatch base and position
        self.un_match_dict = {} 
        pos = 0

        pam22 = self.sgRNA[21]
        try:
            pam23 = self.sgRNA[22]
        except IndexError:
            pam23 = "" 

        for i in range(len(un_matchBase_lst)):
            curPos = int(matchPos_lst[i])
            # the base of the guide RNA corresponding to the pos
            rBase = self.sgRNA[curPos+pos]
            # the base of the genome base, note to reverse the base 
            dBase = un_matchBase_lst[i] 
            dBaseRev = self.baseRev_dic[dBase]

            pos = pos + curPos + 1 
            if pos <=20:
                basePair = "%s:%s" % (rBase, dBaseRev)
                self.un_match_dict.setdefault(pos, basePair)
            elif pos == 21:
                continue 
            elif pos == 22:
                pam22 = dBase
            elif pos == 23:
                pam23 = dBase
        self.pamDin = pam22 + pam23

    def scoreSum(self):
        self.score = 1
        for pos, basePair in self.un_match_dict.items():
             try:
                 self.score = self.score * self.offScore_dict[pos][basePair]
             except:
                 continue
        self.score = self.score * self.pamOffScore_dic[self.pamDin]
        # print(self.score)



"""
sgRNA = "TACAAGCGAAGAGGAAAGATTGG"
cigr = '18T5'

res = offPredictor(sgRNA, cigr)
print res.score
"""
