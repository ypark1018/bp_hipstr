#calculate the heterozygosity at each loci

import sys

allele_file = sys.argv[1]

alleles = []
with open(allele_file) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        alleles.append(a)

print "CHR" + "\t" + "N_ALLELE" + "\t" + "HETEROZYGOSITY"
for chrm in alleles:
    het = 1
    for m,freq in enumerate(chrm):
        if m < 4:
            continue
        else:
            alle = freq.split(":")
            het -= float(alle[1])**2
    print ':'.join(chrm[:2]) + '\t' + chrm[2] + '\t' + str(het)
