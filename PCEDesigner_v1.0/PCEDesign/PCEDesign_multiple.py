####################################################
# Written by Xiaoli Shi 2024.12
# Developed by modifying run_design_pipeline.py from Xianrong Xie
# Paring Twin pegRNAs
# Usage: python PCEDesign_multiple.py ./example/test_PCE.fa

import sys
sys.path.append("../module_pipeline")
sys.path.append("../../db/IRGSP1.0")
import re
import os

blast = \
       "/usr/bin/blastn"   # replaced by your own blastn
batmap = \
        "../BatMis-3.00/bin/batmap"   # replaced by your own batmap


import collections
import seqDeal
import secondAanlysis_PCE
from blastToNCBI import ncbiBlast
from searchPAM_PCE import targetPAM
from readSam_PCE import samReader
from pairing_PCE import SequencePairFinder

(xx, fa_fn) = sys.argv

### pre: DB dir
db_dir = "../../db/IRGSP1.0"  # replace your database directory

### step00: prepare directory for files 
res_dir = "result_PCE"
fa_dir = "{}/fasta".format(res_dir)
tar_dir = "{}/targets".format(res_dir)
blast_dir = "{}/blast".format(res_dir)
sam_dir = "{}/sam".format(res_dir)
out_dir = "{}/results".format(res_dir)

if not os.path.exists(res_dir):
    os.mkdir(res_dir)
if not os.path.exists(fa_dir):
    os.mkdir(fa_dir)
if not os.path.exists(tar_dir):
    os.mkdir(tar_dir)
if not os.path.exists(blast_dir):
    os.mkdir(blast_dir)
if not os.path.exists(sam_dir):
    os.mkdir(sam_dir)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

### step01: split multiple seq contained fasta 
seqName_lst = []
seq_dict = collections.OrderedDict() 
# read fasta
seqName = "Seq"
seq_lst = []
recomb_type_dict = {}

for line in open(fa_fn,"r"):
    line = line.strip()
    if line.startswith(">"):
        if len(seq_lst) != 0:
            seq_dict[seqName] = ''.join(seq_lst)
        word = line.split()
        seqName = word[0][1:]
        recomb_type_dict[seqName] = word[1]
        seqName_lst.append(seqName)
        seq_lst = []
        continue
    seqEachLine = line.upper()
    matchBase = re.findall("[ATCGRYMKSWHBVDN]", seqEachLine)
    seq_lst.append(''.join(matchBase))
seq_dict[seqName] = ''.join(seq_lst)

# save files
for seqName, seqStr in seq_dict.items():
    seq_fn = "{}/{}.fa".format(fa_dir, seqName)
    seq_op = open(seq_fn, "w")
    seq_op.write(">{}\n{}".format(seqName, seqStr))
    seq_op.close()

# sgRNA seq
sgRNA_file = "sgRNA_RT_seq.tab"
par_fn = "par_PCEdesign.conf"
sgRNA_dict = {}

if not os.path.exists(sgRNA_file):
    print(f"Error: file '{sgRNA_file}' does not exist!")
    sys.exit(1)

with open(sgRNA_file,'r') as f:
    lines = f.readlines()

for seqName, recomb_type in recomb_type_dict.items():
    sgRNA_dict[seqName] = []
    for line in lines:
        line = line.strip()
        word = line.split()
        if word[0] == recomb_type:
            # print(word[2].upper())
            sgRNA_dict[seqName].append(word[2].upper())

### step02: read parameter configure file
par_lst = []  #genome, pam_seq, pam_direct, target_len
for line in open(par_fn, "r"):
    if line.startswith("#"):
        continue
    lst = line.strip().split()
    par_lst.append(lst[0])
    par_lst.append(lst[1])
    par_lst.append(lst[2])
    par_lst.append(int(lst[3]))

### step03: design targets of each seq
genome = par_lst[0]
PAM = par_lst[1]
pam_direct = par_lst[2]
guide_len = par_lst[3]

for seqName in seqName_lst:
    seq_fn = "{}/{}.fa".format(fa_dir, seqName)
    ### step04: get the genomic position of the seq
    blastGo = ncbiBlast(blast, seqName, genome, seq_fn, blast_dir)
    hit_pos = list(blastGo.blastInfo_dict.keys())[0]
    hitSeqChrom = hit_pos.split(":")[0]
    hitSeqStart = int(hit_pos.split(":")[1].split("-")[0])
    hitSeqEnd = int(hit_pos.split(":")[1].split("-")[1])

    ### step05: search PAM and save as fasta file
    seqDict = seqDeal.readFa(seq_fn)
    pamRes = targetPAM(seqDict, PAM, pam_direct, guide_len)
    rawPamInfo_dict = pamRes.rawPamInfo_dict 
    targetFileName= "{}/{}_target.fa".format(tar_dir, seqName)
    targetFile = open(targetFileName, "w")
    pamInfo_dict = {}
    tarID = 0
    # value is one list
    for k in sorted(list(rawPamInfo_dict.keys())):
        curPos = k
        curInfo_lst = rawPamInfo_dict[k]
        # print(curPos, curInfo_lst)
        # More than one targets for current position k, keep the first two
        if len(curInfo_lst) > 3:
            curInfo_lst2 = curInfo_lst[0:3]
            seqEach =  curInfo_lst2[1]
            # print("seqEach: ", curInfo_lst2[1])
            pamInfo_dict.setdefault(tarID, curInfo_lst2).append(curPos)
            seqCont = ">%s\n%s\n" % (tarID, seqEach)
            targetFile.write(seqCont)
            tarID += 1
            curInfo_lst = curInfo_lst[3:]
        seqEach = curInfo_lst[1]
        pamInfo_dict.setdefault(tarID, curInfo_lst).append(curPos)
        seqCont = ">%s\n%s\n" % (tarID, seqEach)
        targetFile.write(seqCont)
        tarID += 1
    targetFile.close()
    '''
    for k, v in pamInfo_dict.items():
        print(k, v)
    exit(0)
    '''
    ### step06: run the target-sgRNA pairing
    ## pair_dic: key1: tarID, key2: DNA_str, value: [(ref_start, ref_end),(tar_start, tar_end)]
    pair_dic = secondAanlysis_PCE.outPair(pamInfo_dict, sgRNA_dict[seqName], PAM)

    ### step07: run the offtarget prediction
    if genome != "blank":
        batmapOutFileName = "{}/{}.sam".format(sam_dir, seqName)
        genome_file = "{}/{}.fa".format(db_dir, genome)
        batmap_shell = \
            "%s -g %s -q %s -o %s -n 5 -m 20" % (batmap, genome_file, targetFileName, batmapOutFileName)
        print(batmap_shell)
        os.system(batmap_shell)
        # read the batmap result file --sam
        samRes = samReader(genome, hitSeqChrom, hitSeqStart, hitSeqEnd, PAM, guide_len, pamInfo_dict, batmapOutFileName)
        hit_gene = samRes.targetSeqGene
        # print("hit_gene ", hit_gene)
        targetAnnoteInfo_dict = samRes.targetAnnoteInfo_dic
        '''
        print("targetAnnoteInfo_dict")
        for k, v in targetAnnoteInfo_dict.items():
            print(k, v)
        '''
        offTarget_dict = samRes.alignInfo
        '''
        print("alignInfo")
        for k1, v1 in offTarget_dict.items():
            for k2, v2 in v1.items():
                print(k1, k2, v2)
        exit(0)
        '''
        annoteInfo_dict = samRes.annoteInfo_dic
        '''
        print("annoteInfo_dict")
        for k1, v1 in annoteInfo_dict.items():
            for k2, v2 in v1.items():
                print(k1, k2, v2)
        '''
        offDisplay_dict = samRes.offDisplay_dic
        '''
        print("offDisplay_dict")
        for k, v in offDisplay_dict.items():
            print(k, v)
        exit(0)
        '''
    else:
        hit_gene = "NA"
        targetAnnoteInfo_dict = {}
        offTarget_dict = {}
        annoteInfo_dict = {}
        offDisplay_dict = {}

    ### step08: write result table
    tar_result_fn = "{}/{}.result.txt".format(out_dir, seqName)
    tar_result_op = open(tar_result_fn, "w")

    tar_head = "#tarID\tPosition\tStrand\tTarget_seq\tGC\tRegion\tPredicted_offtarget_value\n"
    tar_result_op.write(tar_head)
    del_tarID_lst = []
    for tarID, tar_info_lst in pamInfo_dict.items():
        tarDirect = tar_info_lst[0]
        tarSeq = tar_info_lst[1]
        tarGC = tar_info_lst[2]
        tarPos = tar_info_lst[3]
        annote = targetAnnoteInfo_dict[tarID]
        offDisplay = offDisplay_dict[str(tarID)]
        # pair sgRNA info; pegRNA secondary structure prediction
        cur_pair_info_dict = pair_dic[str(tarID)]
        cur_pair_len_lst = sorted([len(i) for i in cur_pair_info_dict.keys()], reverse=True)
        cur_pair_len = ""
        if len(cur_pair_len_lst)!=0:
            cur_pair_len = ",".join(["{}nt".format(str(i)) for i in cur_pair_len_lst])
        # judge bad target warning:add warning in 3 cases: 1) GC content <=20 or >= 80. 2)protospacer and scaffold+RT >=8 4) max offtarget score >= 0.7
        try:
            warn = seqDeal.judge_bad(tarSeq, tarGC, cur_pair_len_lst, offDisplay) 
        except:
            warn = "" 
        # write tar result
        if '!' in warn:
            del_tarID_lst.append(tarID)
            continue
        cont = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(tarID, tarPos, tarDirect, tarSeq, tarGC, annote, offDisplay)
        tar_result_op.write(cont)
    tar_result_op.close()

    ### step09: Twin peg pairing
    pairs_finder = SequencePairFinder(pamInfo_dict, offDisplay_dict, pair_dic)
    pairs_on_dict = pairs_finder.find_ontarget_pairs()
    pairs_lst = []
    for key, value in pairs_on_dict.items():
        pairs_lst.append(value)
    sorted_pairs_lst = sorted(pairs_lst, key = lambda x:x[14])
    pair_result_fn = "{}/{}.TwinPair.txt".format(out_dir, seqName)
    pair_result_op = open(pair_result_fn, "w")
    pair_head = "Target_seq1\tStart1\tEnd1\tStrand1\tGC1\tMax_offScore1\t2nd_strucLen1\tTarget_seq2\tStart2\tEnd2\tStrand2\tGC2\tMax_offScore2\t2nd_strucLen2\tTwin_dist\n"
    pair_result_op.write(pair_head)

    for pair_item in sorted_pairs_lst:
        pair_result_op.write("\t".join(map(str, pair_item))+'\n')
    pair_result_op.close()
    print(f"Done, paired sequence has been deposited into {pair_result_fn} file.")

    ### End: remove temporary files
    # blastFileName = "{}/{}.blst".format(blast_dir, seqName)
    os.remove(seq_fn)
    os.remove(targetFileName)
    os.remove(batmapOutFileName)
    # os.remove(blastFileName)

