#this program parses through the vcf file and filters out STRs with irregular patterns
#currently includes pattern 3

import sys

alleles = []
with open(sys.argv[1]) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        alleles.append(a)

#STR_genotype will contain all genotypes from vcf
STR_genotype = []
typeone = 0
typetwo = 0
for locus in alleles:
    #if one == 0 the genotype is composed of repeat of one nucleotide
    one = 0
    if '#' not in locus[0]:
        #genotype will contain genotypes for individual STR
        genotype = []
        for STRinfo in locus[:4]:
            if STRinfo != locus[2]:
                genotype.append(STRinfo)
        alt = locus[4].split(",")
        for a in alt:
            genotype.append(a)
        for l in genotype[2:]:
            for nuc in l:            
                if nuc != l[0]:
                    one += 1
                    break
        #for now if genotype is only composed of only one nucleotide, break
        if one == 0:
           typeone += 1
        #if one != 0:
        else:
            #if not continue parsing
            #assum the longest number of repeat is units of 5 nucleotides
            repnum = 2
            #if pattern == 1 then the genotype has irregular STR patterns
            pattern = 0
            while repnum <= 5:
                rep = genotype[2][:repnum]
                start = 0
                end = repnum
                while end <= len(genotype[2]):
                    nuc = genotype[2][start:end]
                    if rep == nuc:
                        start += repnum
                        end += repnum
                    else:
                        pattern = 1
                        break
                if pattern == 1 and repnum < 5:
                    repnum += 1
                    pattern = 0
                else:
                    break
            if pattern == 0:
                STR_genotype.append(genotype)
            else:
                typetwo += 1

#print genotype
for STR in STR_genotype:
    print '\t'.join(STR)
#print typeone
#print typetwo
