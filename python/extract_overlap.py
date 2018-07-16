from collections import defaultdict
import sys

STR_file = sys.argv[1]
indel_file = sys.argv[2]


STR = []
with open(STR_file) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        a = a[:3] + [a[3].split(";")[2][4:]]
        #int(STR_locus[3].split(";")[2][4:])
        STR.append(a)


indel = defaultdict(list)
with open(indel_file) as g:
    reflines = g.readlines()
    for line in reflines:
        b = line.strip().split("\t")
        indel[b[0]].append(b)

overlap_output = sys.argv[2] + ".overlap"
flank_output = syhs.argv[2] + ".flank"

ov_output = open(overlap_output, 'w+')
ov_output.truncate()

fl_output = open(flank_output, 'w+')
fl_output.truncate()
#ov_output.write("#SV\tPOS\tREF\tEND\tSV2\tPOS\tREF\tALT\n")


for STR_locus in STR:
    str_length = int(STR_locus[3]) - int(STR_locus[1])
    for indel_locus in indel[STR_locus[0]]:
        indel_alleles = [indel_locus[2]] + indel_locus[3].split(',')
        allele_end = int(indel_locus[1]) + max(len(allele) for allele in indel_alleles)

        if allele_end >= int(STR_locus[1]) and int(indel_locus[1]) <= int(STR_locus[3]):
            fl_output.write('\t'.join(indel_locus[:2]) + '\n')

        elif allele_end >= (int(STR_locus[1]) - 1000) and int(indel_locus[1]) <= (int(STR_locus[3]) + 1000):
            ov_output.write('\t'.join(indel_locus[:2]) + '\n')
