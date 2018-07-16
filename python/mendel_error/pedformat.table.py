#take the wide format and get table ready for pedfile conversion
import sys

allele_length = sys.argv[1]

alleles = []
with open(allele_length) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        alleles.append(a)

item = 0
new_data = []
mark = []
for y in zip(*alleles):
    STR = list(y)
    if item == 0:
        item += 1
        new_data.append(STR)
        continue
    elif item == 1:
        STR[0] = STR[0].replace("_first_chr", "")
        mark = STR
        item += 1
    else:
        for n,x in enumerate(STR):
            if "DNA" in x:
                continue
            else:
                if mark[n] == "NA":
                    mark[n] = 0
                if x == "NA":
                    x = 0
                mark[n] = str(mark[n]) + " " + str(x)
        new_data.append(mark)
        item -= 1


for chrom in zip(*new_data):
    print '\t'.join(chrom)
