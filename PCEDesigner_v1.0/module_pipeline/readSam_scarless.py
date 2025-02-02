import re
import sys
import seqDeal
from predictOffScore import offPredictor
# from annotByPos_scarless import gffJudge 

strand_dict = {"0": "+", 
               "16": "-"}

region_seq = ["CDS", "five_prime_UTR", "three_prime_UTR", "exon", "pos"]

class samReader:
    def __init__(self, genome, guide_len, pamInfo_dict, sam, edit_type, terminal):
        self.genome = genome
        # self.hitSeqChrom = hitSeqChrom
        # self.seqDirect = "+"
        # if hitSeqStart > hitSeqEnd:
        # self.seqDirect = "-"

        # self.hitSeqStart = min(hitSeqStart, hitSeqEnd)
        # self.hitSeqEnd = max(hitSeqStart, hitSeqEnd)

        # self.pam = pam
        # self.tarSeq_len = int(guide_len) + len(self.pam) 
        self.pamInfo_dict = pamInfo_dict 
        self.sam = sam
        self.edit_type = edit_type
        self.terminal = terminal


        # self.targetAnnotePos_dic = {}
        # self.targetAnnoteInfo_dic = {}
        self.annotePos_dic = {}
        self.annotePosIn_dic = {}  # pos dict for annotation
        # self.annoteInfo_dic = {}
        self.alignInfo = {}
        self.offDisplay_dic = {}
        self.maxScore = "NA"    # the offScore to be displayed in HTML

        self.readFile()
        self.annotePos()

    def readFile(self):
        offID = 1
        targetID_bf = "" 
        self.eachOff_dic = {}
        self.offScore_dic = {}
        for l in open(self.sam, "r"):
            ls = l.rstrip()
            if re.search("^@", ls):
                continue
            lst = ls.split("\t")
            idx = lst[0]  # the seqID
            # start to align a new spacer
            if targetID_bf != idx:
                # save the alignment of before spacer
                if len(targetID_bf) != 0:
                    if len(list(self.eachOff_dic.keys()))!=0:
                        newEachOff_dic = self.sortOffByScore()
                        self.alignInfo.setdefault(self.tarID, newEachOff_dic)
                    else:
                        self.alignInfo.setdefault(self.tarID, {})
                offID = 1
                targetID_bf = idx
                self.eachOff_dic = {}
                self.offScore_dic = {}

            self.tarID = idx

            if lst[1]=='4': # some hit seqs have hit to a unplaced ref
                continue
            self.strand = strand_dict[lst[1]]  # the align seqStrand
            self.chrom = lst[2]
            self.pos = int(lst[3])
            self.cigr = lst[-1].split(":")[-1] # the cigr str

            # revserse the cigr
            if self.strand == "-":
                self.cigr = self.revCIGR()

            # cal the offTarget score  
            # used for multiple pams 
            sgRNA = self.pamInfo_dict[int(idx)][1] 
            sgRNA_pos = self.pamInfo_dict[int(idx)][-1]
            offScore = 'NA'
            # if self.pam[1:] == "GG" or self.pam == "NG": 
            res = offPredictor(sgRNA, self.cigr)
            offScore = res.score 

            # the position region (exon, cds, etc..) of off-target
            self.annotePosIn_dic.setdefault(self.chrom, []).append(self.pos)
            '''
            if self.chrom == self.hitSeqChrom:
                if self.seqDirect == "+":
                    tarPos = self.hitSeqStart + sgRNA_pos - 1
                else:
                    tarSeq_len = len(sgRNA)
                    tarPos = self.hitSeqEnd - sgRNA_pos + 1 - tarSeq_len + 1
                # print tarPos, self.pos
                if self.hitSeqEnd>=self.pos>=self.hitSeqStart:
                    self.targetAnnotePos_dic.setdefault(self.tarID, [self.chrom, self.pos])
            '''
            contLst = [self.chrom, self.pos, self.strand, self.cigr, offScore] 
            # must sort the offID by score first
            self.eachOff_dic.setdefault(offID, contLst) 
            self.offScore_dic.setdefault(offScore, []).append(offID)
            offID += 1
        # store the last spacer
        if len(list(self.eachOff_dic.keys()))!=0:
            newEachOff_dic = self.sortOffByScore()
            self.alignInfo.setdefault(self.tarID, newEachOff_dic)
        else:
            self.alignInfo.setdefault(self.tarID, {})


    def sortOffByScore(self):
        newEachOff_dic = {}
        offScoreRev_lst = []
        offScore_lst = reversed(sorted(list(self.offScore_dic.keys())))
        newID = 1
        for i in offScore_lst:
            offID_lst = self.offScore_dic[i]
            for offID in offID_lst:
                cont_lst = self.eachOff_dic[offID]

                # add the annote dic corresponding
                self.annotePos_dic.setdefault(cont_lst[0], {}).setdefault(cont_lst[1], {}).setdefault(newID, []).append(self.tarID)
                newEachOff_dic.setdefault(newID, []).extend(cont_lst)
                newID += 1
        return newEachOff_dic

    def revCIGR(self): 
        matBase = re.compile("[ATCGRYMKSWHBVDN]$") 
        # remove the 0, if the cigr is start with 0
        if re.search("^0", self.cigr):
            self.cigr = self.cigr[1:]
        # plus 0 to the end, if the cigr is end with [ATCG, etc]
        if re.search(matBase, self.cigr):
            self.cigr += "0"
        numLst = re.findall("\d+", self.cigr)
        alpLst = re.findall("[ATCGRYMKSWHBVDN]", self.cigr)
        numLst.reverse()
        alpLst.reverse()
        newCigrLst = []
        for i in range(len(numLst)):
            curNum = numLst[i]
            newCigrLst.append(curNum)
            if i <= (len(alpLst) - 1):
                curAlp = seqDeal.revSeq(alpLst[i])
                newCigrLst.append(curAlp)
        newcigr = ''.join(newCigrLst)
        return newcigr

    def filterUnPAM(self):
        filterPass = True
        if self.strand == "+":
            mat = re.search("\d+$", self.cigr)
            if not mat:
                filterPass = False
            else:
                matNum = mat.group()
                if int(matNum) < 3:
                    filterPass = False
        else:
            mat = re.search("^\d+", self.cigr)
            if not mat:
                filterPass = False
            else:
                matNum = mat.group()
                if int(matNum) < 3:
                    filterPass = False

        return filterPass

    # run the annote program
    def annotePos(self):
        # seqInfo = [self.hitSeqChrom, self.hitSeqStart, self.hitSeqEnd]
        # res = gffJudge(seqInfo, self.annotePosIn_dic, self.genome)
        # res = gffJudge(self.annotePosIn_dic, self.genome)
        # hitRegion_dic = res.hitRegion_dic
        # hitGeneInfo_dic = res.hitGene_info
        # self.annoteTarget(hitGeneInfo_dic)
        # do not annotate the on- and off-target region
        '''
        for chrom in list(hitRegion_dic.keys()):
            for pos in list(hitRegion_dic[chrom].keys()):
                region_info = hitRegion_dic[chrom][pos]
                curOff_dic = self.annotePos_dic[chrom][pos]
                for offID in list(curOff_dic.keys()):
                    tarID_lst = curOff_dic[offID]
                    for tarID in tarID_lst:
                        if tarID in list(self.targetAnnotePos_dic.keys()):
                            curTargetPos = self.targetAnnotePos_dic[tarID]
                            # if curTargetPos[0]==chrom and self.hitSeqEnd>=curTargetPos[1]>=self.hitSeqStart:
                            #    del self.alignInfo[tarID][offID]
                            #    continue
                        # self.annoteInfo_dic.setdefault(tarID, {}).setdefault(offID, region_info)
        '''
        # get the max off-target score of each tarID
        # on-target is also in self.alignInfo
        for tarID in list(self.alignInfo.keys()):
            tarStrand = self.pamInfo_dict[int(tarID)][0]
            maxScore = "NA"
            secondMaxScore = "NA"
            # if self.pam[1:]=="GG" or self.pam=="NG":
            if len(list(self.alignInfo[tarID].keys()))!=0:
                maxScore_offID = sorted(list(self.alignInfo[tarID].keys()))[0]
                maxScore = "%0.3f" % self.alignInfo[tarID][maxScore_offID][-1] 
                second_maxScore_offID = sorted(list(self.alignInfo[tarID].keys()))[1]
                secondMaxScore = "%0.3f" % self.alignInfo[tarID][second_maxScore_offID][-1]
            # insertion 5' scar: "+" spacer is in reference, "-" is in inserted segment
            # insertion 3' scar: "+" spacer is in inserted segment, "-" is in reference
            if self.edit_type == "Integration":
                if (self.terminal == "5'" and tarStrand == "-") or (self.terminal == "3'" and tarStrand == "+"):
                    self.offDisplay_dic.setdefault(tarID, maxScore)
                else:
                    self.offDisplay_dic.setdefault(tarID, secondMaxScore)
            else:
                self.offDisplay_dic.setdefault(tarID, secondMaxScore)
 
    def annoteTarget(self, hitGene_info):
        self.targetSeqGene = "NA"
        targetSeqGene_lst = []
        targetSeqGene_lst.extend(list(hitGene_info.keys()))
        if len(targetSeqGene_lst) != 0:
            self.targetSeqGene = '; '.join(targetSeqGene_lst)

        for tarID in list(self.pamInfo_dict.keys()):
            if self.seqDirect == "+":
                pos = self.pamInfo_dict[tarID][-1] + self.hitSeqStart - 1
            else:
                tarSeq_len = len(pamInfo_dict[tarID][1])
                pos = self.hitSeqEnd - self.pamInfo_dict[tarID][-1] + 1 - tarSeq_len + 1

            anote_lst = []
            for gene in list(hitGene_info.keys()):
                trans_dic = hitGene_info[gene]
                for trans_id in list(trans_dic.keys()):
                    info_dic = trans_dic[trans_id]
                    for t in list(info_dic.keys()):
                        pos_lst = info_dic[t]
                        for each_lst in pos_lst:
                            start = each_lst[0]
                            end = each_lst[1]
                            if start<=pos<=end:
                                anote_lst.append(t)
            '''
            # judge the hit region at last 
            if len(anote_lst)==0:
                region = "intergenic"
                self.targetAnnoteInfo_dic.setdefault(tarID, region)
                continue
            for k in region_seq:
                if k in anote_lst:
                    region = k
                    break
            if region == "pos":
                region = "intron"

            self.targetAnnoteInfo_dic.setdefault(tarID, region)
            '''

"""
(xx, f) = sys.argv
res = samReader(f)
res
"""
