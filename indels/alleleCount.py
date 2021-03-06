#count the number of alternative alleles in each loci from STR_indel.txt
from collections import defaultdict
import sys

STR_file = sys.argv[1]

STR = []
with open(STR_file) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        STR.append(a)

for STR_locus in STR:
    slocus = STR_locus[3].split(",")
    slocusCount = len(slocus) + 1
#    if "*" in slocus:
#        slocusCount -= 1
    ilocus = STR_locus[7].split(",")
    ilocusCount = len(ilocus) + 1
    if "*" in ilocus:
        ilocusCount -= 1
#        print str(pos) + " found *"
    difference = slocusCount - ilocusCount
    print '\t'.join(STR_locus[:2]) + '\t' + str(slocusCount) + '\t' +  '\t'.join(STR_locus[4:6]) + '\t' + str(ilocusCount) + '\t' + str(difference)
