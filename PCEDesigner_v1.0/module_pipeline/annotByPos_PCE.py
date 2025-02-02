import re
import sys
sys.path.append("../../db/IRGSP1.0/") 
from gffInfo_IRGSP import gffInfo
region_seq = ["CDS", "five_prime_UTR", "three_prime_UTR", "exon", "pos"]


class gffJudge:
    def __init__(self, seqInfo, pos_dic, genome):
        self.seqInfo = seqInfo
        self.pos_dic = pos_dic 
        self.genome = genome 
        gffRes = gffInfo()

        self.gff_dic = gffRes.gff_dic
        self.gene_type_dic = gffRes.gene_type_dic
        self.gene_pos_dic = gffRes.gene_pos_dic
        # self.gene_strand_dic = gffRes.gene_strand_dic
        self.bef_gene_dic = gffRes.bef_gene_dic
        # print(self.gff_dic) 

        if len(seqInfo)!=0:
            self.hitGene_info = {} 
            self.judgeSeq()

        self.hitRegion_dic = {} 
        self.judgePos()

    def gffReader(self):
        self.gff_dic = {}
        self.gene_type_dic = {}
        self.gene_pos_dic = {}
        # self.gene_strand_dic = {}
        # self.bef_gene_dic = {}
        bef_gene = ''
        # the start point of a new gene or ncRNA

        for l in open(self.gff, "r"):
            ls = l.rstrip()
            if not re.search(r"^([0-9])", ls):
                continue

            lst = ls.split("\t")
            chrom = lst[0]
            seq_type = lst[2]
            if seq_type == "chromosome":
                continue 
            seq_begin = int(lst[3])
            seq_end = int(lst[4])
            seq_strand = lst[6]

            if seq_type == "gene" or seq_type == "ncRNA_gene":
                gene_id = lst[8].split(";")[0].split(":")[1]
                self.gene_type_dic.setdefault(chrom , {}).setdefault(gene_id, seq_type) 
                self.gene_pos_dic.setdefault(chrom , {}).setdefault(gene_id, [seq_begin, seq_end])
                # self.gene_strand_dic.setdefault(chrom , {}).setdefault(gene_id, seq_strand)
                # self.bef_gene_dic.setdefault(chrom , {}).setdefault(gene_id, bef_gene) 
                bef_gene = gene_id
            elif seq_type == "mRNA" or seq_type == "transcript":
                trans_id = lst[8].split(";")[0].split(":")[1]
                self.gff_dic.setdefault(chrom , {}).setdefault(gene_id, {}).setdefault(trans_id, {}).setdefault("pos", []).append([seq_begin, seq_end])
            else:
                self.gff_dic.setdefault(chrom , {}).setdefault(gene_id, {}).setdefault(trans_id, {}).setdefault(seq_type, []).append([seq_begin, seq_end])

    def judgeSeq(self):
        chrom = self.seqInfo[0]
        # remove the "Chr" etc.. chrac in chrom 
        if self.genome =="IRGSP1.0":
            if re.search("^[cC]hr", chrom):
                chrom = chrom[3:]
            if re.search("^0", chrom):
                chrom = chrom[1:]

        seqSt = int(self.seqInfo[1])
        seqEnd = int(self.seqInfo[2])

        for gene_id in self.gff_dic[chrom].keys():
            gene_pos = self.gene_pos_dic[chrom][gene_id]
            if seqEnd < gene_pos[0] or seqSt > gene_pos[1]: 
                continue
            trans_dic = self.gff_dic[chrom][gene_id] 
            self.hitGene_info.setdefault(gene_id, trans_dic)  


    def judgePos(self):
        for chrom in self.pos_dic.keys():
            pos_lst = sorted(self.pos_dic[chrom])
            curRegion_dic = self.judgePosEach(chrom, pos_lst) 
            self.hitRegion_dic.setdefault(chrom, curRegion_dic)

    def judgePosEach(self, chrom, posQuery_lst):
        # remove the "Chr" etc.. chrac in chrom 
        if self.genome =="IRGSP1.0":
            if re.search("^[cC]hr", chrom):
                chrom = chrom[3:]
            if re.search("^0", chrom):
                chrom = chrom[1:]
        anote_dic = {}
        for i in posQuery_lst:
            anote_dic.setdefault(i, {})

        # sequence of some species have no annotation info
        # add 2017.12.11
        if chrom not in self.gff_dic.keys():
            cur_hit_region_dic = {} 
            for pos in posQuery_lst:
                gene = ''
                region = "intergenic"
                hit_region = [gene, region] 
                cur_hit_region_dic.setdefault(pos, []).extend(hit_region)
            return cur_hit_region_dic
 
        for gene_id in self.gff_dic[chrom].keys():
            # gene_type = self.gene_type_dic[chrom][gene_id]
            gene_pos = self.gene_pos_dic[chrom][gene_id]
            # gene_strand = self.gene_strand_dic[chrom][gene_id]
            # bef_gene = self.bef_gene_dic[chrom][gene_id]
            for pos in posQuery_lst:
                if pos < gene_pos[0] or pos > gene_pos[1]:
                    continue

                # the dict of trans 
                trans_dic = self.gff_dic[chrom][gene_id] 
                for trans_id in trans_dic.keys():
                    info_dic = trans_dic[trans_id] 
                    for t in info_dic.keys():
                        pos_lst = info_dic[t]
                        for each_lst in pos_lst:
                            start = each_lst[0]
                            end = each_lst[1]
                            if start<=pos<=end:
                                anote_dic.setdefault(pos, {}).setdefault(gene_id, {}).setdefault(trans_id, []).append(t)

        # judge the hit region at last 
        cur_hit_region_dic = {} 

        for pos in anote_dic.keys():
            each_anote_dic = anote_dic[pos]

            if len(each_anote_dic) == 0:
                gene = ''
                region = "intergenic"
            else:
                # put all hit trans info to a single lst
                each_region_lst = []
                gene = ''
                for i in each_anote_dic.keys():
                    gene = i
                    for trans in each_anote_dic[gene]:
                        each_region_lst.extend(each_anote_dic[gene][trans])
                # judge display most important hit region
                for k in region_seq:
                    if k in each_region_lst:
                        region = k
                        break
                if region == "pos":
                    region = "intron" 

            hit_region = [gene, region] 
            cur_hit_region_dic.setdefault(pos, []).extend(hit_region)
        return cur_hit_region_dic



"""
posDic = {"chr01": [10000, 2000300, 10275, 10430], "chr12": [5400,100,1090200,2682,4060004]} 

genome = "IRGSP"

# (xx, gff)=sys.argv
res = gffJudge(posDic, genome)
hitRegion_dic = res.hitRegion_dic
print hitRegion_dic 
"""

# print gff_dic
# print gene_type_dic
