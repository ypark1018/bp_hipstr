#get the minor allele frequencies at each loci


import sys

allele_file = sys.argv[1]

alleles = []
with open(allele_file) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        alleles.append(a)

indices = []
a1 = []
#a2 will hold the minor allele frequency
a2 = []
for n,chrm in enumerate(alleles):
    if n == 0:
        continue
    else:
        a1 = []
        a2 = []
        for m,freq in enumerate(chrm):
            if m < 4:
                continue
            else:
                alle = freq.split(":")
                if not a1:
                    a1 = alle
                elif alle[1] >= a1[1]:
                    a2 = a1
                    a1 = alle
                elif alle[1] < a1[1]:
                    if not a2:
                        a2 = alle
                    elif alle[1] > a2[1]:
                        a2 = alle
        if a1[1] == 1:
            a2[1] = "NAN"
        print ':'.join(chrm[:2]) + '\t' + a2[1] #'\t'.join(a2)
