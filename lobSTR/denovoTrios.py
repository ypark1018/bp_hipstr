#look for sets of complete trios and look for de novo mutations
#chr1.processed.recode.GT.FORMAT
#CR004.key
#CR004_19indivi.ped.txt

import sys

genotypes = []
with open(sys.argv[1]) as f:
    alengths = f.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        genotypes.append(a)

ind_keys = {}
with open(sys.argv[2]) as g:
    alengths = g.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        ind_keys[a[0]] = a[1]

pedigree = []
with open(sys.argv[3]) as h:
    alengths = h.readlines()
    for allele in alengths:
        a = allele.strip().split(" ")
        pedigree.append(a)

#ped_dict = {}
#for individual in pedigree:
#    try:
#        ped_dict[individual[0]].append(individual[1])
#    except KeyError:
#        ped_dict[individual[0]] = [individual[1]]

class Tree(object):
    def __init__(self):
        self.dad = None
        self.mom = None
        self.identity = None
        self.origin = False
        self.checked = False
        self.score = 0

families = {}
for individual in pedigree:
    root = Tree()
    root.identity = individual[1]
    root.dad = individual[2]
    root.mom = individual[3]
    families[root.identity] = root

def get_levels(node):
    node.checked = True
    if node.mom == "0" or node.mom == "0":
        if node.origin == True:
            node.origin = False
        #print "ancestor: " + node.identity
        return 1
    else:
        mom = families[node.mom]
        dad = families[node.dad]
        
        mom_score = 1 + get_levels(mom)
        dad_score = 1 + get_levels(dad)
        if mom_score > dad_score:
            if mom.origin == True:
                mom.origin = False
            #print "path: " + node.mom
            return mom_score
        else:
            if dad.origin == True:
                dad.origin = False
            #print "path: " + node.dad
            return dad_score

level = 0
for individual in families.keys():
    ind_node = families[individual]
    if not ind_node.checked:
        ind_node.origin = True;
        level = get_levels(families[individual])
        ind_node.score = level
        #print "origin: " + ind_node.identity
        #print level

base_node = []
for individual in families.keys():
    ind_node = families[individual]
    if ind_node.origin and ind_node.score > 2:
        base_node.append(ind_node.identity)
        #print ind_node.identity + " " + str(ind_node.score)

for n, location in enumerate(genotypes):
    genotypes[n] = [location[0] + ":" + location[1]] + location[2:]
#print genotypes

#assign index numbers for each individual so we can reference them through index
#go through the list by each STR locux
identity_num = {}
for n,chrom in enumerate(genotypes[0]):
    identity_num[chrom] = n

first = 0
for chrom in genotypes:
    mismatch = []
    parents = []
    loci_output = ""
    loci_output += chrom[0] + "\n"
    if first == 0:
        first = 1
        continue
    for node in base_node:
        denovo = 0
        current = node
        mother = (families[node]).mom
        father = (families[node]).dad
        diverge = []
        diverge.append("0")
        while current != "0":
            #print current + "\t" + father + "\t" + mother
            count = 3
            try:
                genotype = chrom[identity_num[ind_keys[current]]].split('/')
            except KeyError:
                count -= 1
                genotype = [-2, -2]
            try:
                mom_gene = chrom[identity_num[ind_keys[mother]]].split('/')
            except KeyError:
                count -= 1
                mom_gene = [-1, -1]
            try:
                dad_gene = chrom[identity_num[ind_keys[father]]].split('/')
            except KeyError:
                count -= 1
                dad_gene = [-1, -1]
            if count >= 2 and ("." not in genotype) and ("." not in mom_gene) and ("." not in dad_gene):
                #loci_output += current + "\t" + str(genotype) + "\n"
                #loci_output += mother + "\t" + str(mom_gene) + "\t" + father + "\t" + str(dad_gene) + "\n"
                parents.append([mother, str(mom_gene), father, str(dad_gene)])
                combined = mom_gene + dad_gene
                if (genotype[0] not in combined) and (genotype[1] not in combined):
                    #loci_output += "MISMATCH********************************\n"
                    mismatch.append([current, genotype])
                    
            if (families[current]).dad != "0":
                diverge.append((families[current]).dad)
            #print diverge
            if mother == "0" or father == "0":
                current = diverge.pop()
                if current == "0":
                    break
                mother = (families[current]).mom
                father = (families[current]).dad
            else:
                current = (families[current]).mom
                mother = (families[current]).mom
                father = (families[current]).dad
    #if "MISMATCH" in loci_output:
    #    print loci_output
    if mismatch:
       # print chrom[0]
        for a in mismatch:
            #print a[0]
            for b in parents:
                #print b
                if a[0] in b:
                    print chrom[0]
                    print "DENOVO FOUND*******************************"
                    print a
                    print b
