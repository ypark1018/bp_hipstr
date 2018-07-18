#determine which STRs were detected as indels in GATK
#ignore indels of length 1
from collections import defaultdict
import sys

STR_file = sys.argv[1]
indel_file = sys.argv[2]


STR = []
with open(STR_file) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        STR.append(a)

indel = defaultdict(list)
with open(indel_file) as g:
    reflines = g.readlines()
    for line in reflines:
        b = line.strip().split("\t")
        indel[b[0]].append(b)

#print STR

#offset = 10
for STR_locus in STR:
    for indel_locus in indel[STR_locus[0]]:
        offset = 100 #min(len(STR_locus[2]), len(indel_locus[2]))
        if offset <= 1:
            continue
        #elif int(indel_locus[1]) > int(STR_locus[1]) + offset:
        #    break
        elif abs(int(STR_locus[1]) - int(indel_locus[1])) < offset:
            print '\t'.join(STR_locus) + '\t' + "indel:" + '\t'.join(indel_locus)
