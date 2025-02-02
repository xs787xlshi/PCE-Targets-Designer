#####################################################
# Program: blastToNCBI.py
# Function: blast the input sequence to the reference 
# Usage: blastToNCBI.py geneBankID seqFile
# Author: Xianrong Xie
#####################################################

import os 
import sys
import seqDeal


class ncbiBlast:
    def __init__(self, blast, usrid, genome, seqFile, blast_dir):
        self.blast = blast
        self.genome = genome
        self.seqFile = seqFile
        self.blast_dir = blast_dir
        self.blastFileName = "{}/{}.blst".format(self.blast_dir, usrid)
        self.runBlast()
        
        self.blastInfo_dict = {}
        self.readBlast()
     

    def runBlast(self): 
        genome_file = "../../db/%s/%s.fa" % (self.genome, self.genome)
        blast_shell = "%s -query %s -db %s -outfmt '6 std qlen' -num_threads 2 -out %s" % (self.blast, self.seqFile, genome_file, self.blastFileName)
        print(blast_shell)
        os.system(blast_shell)         

    def readBlast(self):
        for l in open(self.blastFileName, "r"):
            lst = l.rstrip().split("\t")
            # hit chr and position
            hit_Chr = lst[1]
            hit_st = int(lst[8])
            hit_end = int(lst[9])
            hit_length = int(lst[3])
            # query hit st and end
            query_st = int(lst[6])
            query_end = int(lst[7])
            query_length = int(lst[-1])
            # identity 
            ident = float(lst[2])
            mis = int(lst[4])
            gap = int(lst[5])

            hit_key = "%s:%s-%s" % (hit_Chr, hit_st, hit_end) 
            self.blastInfo_dict.setdefault(hit_key, []).extend([ident, mis, gap])
            break
