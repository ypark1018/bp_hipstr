import sys

alleles = []
with open(sys.argv[1]) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        alleles.append(a)

for locus in alleles:
    print ':'.join(locus[:2]) + '\t' + locus[2]
