#look for sets of complete trios and look for de novo mutations
#chr1.processed.recode.GT.FORMAT -- file containing only GT info
#FamilyID.conv -- family ID to sample ID conversion in that order
#allfamilies.ped -- pedigree file

import sys
import re
from collections import defaultdict

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
        a = allele.strip().split("\t")
        pedigree.append(a)

qscores = []
with open(sys.argv[4]) as i:
    alengths = i.readlines()
    for allele in alengths:
        a = allele.strip().split("\t")
        qscores.append(a)

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
        self.genotype = None
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
#go through the list by each STR locus
identity_num = {}
for n,chrom in enumerate(genotypes[0]):
    identity_num[chrom] = n





first = 0
denovocount = 0
mendelerror = 0
dtriomendelerror = 0
dtrio = 0
double_parents = defaultdict(list)
doubletrio = 0
ambiguous = 0
dtf = 0
ind_denovo = {}
str_denovo = {}

f_denovo = sys.argv[5] + "_denovo_summary.txt"
denovo_summary = open(f_denovo, 'w+')
denovo_summary.truncate()

denovo_summary.write("CHROM\tDENOVO\tGT\tQSCORE\n")
for chrom in genotypes:
    mismatch = []
    parents = []
    double_parents = defaultdict(list)
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
                genotype = re.split('/|\|',chrom[identity_num[ind_keys[current]]])
                families[current].genotype = genotype
            except KeyError:
                count -= 1
                genotype = [-1, -1]
            try:
                mom_gene = re.split('/|\|',chrom[identity_num[ind_keys[mother]]])
                families[mother].genotype = mom_gene
            except KeyError:
                count -= 1
                mom_gene = [-1, -1]
            
            try:
                dad_gene = re.split('/|\|',chrom[identity_num[ind_keys[father]]])
                families[father].genotype = dad_gene
            except KeyError:
                count -= 1
                dad_gene = [-1, -1]
            
            #if its not a full trio ignore (count == 3)
            if count == 3: #and ("." not in genotype) and ("." not in mom_gene) and ("." not in dad_gene):
                #loci_output += current + "\t" + str(genotype) + "\n"
                #loci_output += mother + "\t" + str(mom_gene) + "\t" + father + "\t" + str(dad_gene) + "\n"
                combined = mom_gene + dad_gene

                #####################################
                #check if both parents are sequenced#
                #####################################
                #if -1 not in combined:
                #parents.append([mother, mom_gene, father, dad_gene, current, genotype])
                allele_flag = 0
                allele_ct = 0
                if genotype[0] in combined:
                    combined.remove(genotype[0])
                    allele_flag += 1
                    allele_ct = 1
                if genotype[1] in combined:
                    #combined.remove(genotype[1])
                    allele_flag += 1
                    allele_ct = 0
                #if (genotype[0] not in combined) and (genotype[1] not in combined) and (genotype != [-1, -1]):

                ################################
                #check if parents are genotyped#
                ################################
                if (genotype != [-1, -1]) and (-1 not in combined) and ("." not in genotype) and ("." not in combined):
                    parents.append([mother, mom_gene, father, dad_gene, current, genotype])

                    if (allele_flag < 2):
                    #loci_output += "MISMATCH********************************\n"
                        mismatch_geno = [current, genotype, genotype[allele_ct]]
                        if mismatch_geno not in mismatch:
                            mismatch.append(mismatch_geno)
                            mendelerror += 1
                    
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

    if dtf == 0:
        temp = 0
        #gen 2
        for par in parents:
            #gen 1
            for parpar in parents:
                if par[0] == parpar[4]:
                    #double_parents = [current0, mom2, dad4, spouse6, kid8]
                    double_parents[par[0]].append(par[:2] + parpar[:4] + par[2:])
                    temp += 1
                    break
                elif par[2] == parpar[4]:
                    #double_parents = [current0, mom2, dad4, spouse6, kid8]
                    double_parents[par[2]].append(par[2:4] + parpar[:4] + par[0:2] + par[4:])
                    temp += 1
                    break
        if temp > doubletrio:
            doubletrio = temp
        #dtf = 1
    #if "MISMATCH" in loci_output:
    if mismatch:
        for m in mismatch:
            denovoflag = 0
            for cur in double_parents.keys():
                if m[0] == cur:
                    dtriomendelerror += 1
                    #one source of denovo may have more than one child
                    for dp in double_parents[cur]:
                        truedn = 0
                        if m[2] in dp[9]:
                            print m
                            print dp

                            denovoflag = 1
                            if dp[1][0] == dp[1][1] or m[2] in dp[7]:
                                category = "AMBIGUOUS"
                            else:
                                category = "TRUE"
                                truedn = 1
                            print "******************************************"
                            print chrom[0] + "\t" + m[0]
                            print "DENOVO:\t" + m[0] + "\t"  + "/".join(m[1]) + "\t" + category
                            print "1st GENERATION:\t" + dp[2] + "\t" + str(dp[3][0]) + "/" + str(dp[3][1]) + "\t" + dp[4] + "\t" + str(dp[5][0]) + "/" + str(dp[5][1])
                            print "2nd GENERATION:\t" + dp[0] + "\t" + str(dp[1][0]) + "/" + str(dp[1][1]) + "\t" + dp[6] + "\t" + str(dp[7][0]) + "/" + str(dp[7][1])
                            print "3rd GENERATION:\t" + dp[8] + "\t" + str(dp[9][0]) + "/" +  str(dp[9][1])
                            person = ind_keys[m[0]]
                            qscore = 0
                            for score in qscores:
                            #print ":".join(score[:2])
                                if ":".join(score[:2]) == chrom[0]:
                                    qscore = score[identity_num[person]]

                                    #we are only counting true denovo mutation not ambiguous
                                    #keep track of now many denovo mutations each person/STR has
                            if category == "TRUE":
                                denovocount += 1
                                ind_denovo[m[0]] = ind_denovo.get(m[0], 0) + 1
                                str_denovo[chrom[0]] = str_denovo.get(chrom[0], 0) + 1
                                s_output = chrom[0] + "\t" + m[0] + "\t" + "/".join(m[1])+ "\t" + str(qscore) + "\n"
                                denovo_summary.write(s_output)
                            else:
                                ambiguous += 1
                            break
                        if truedn:
                            break
                if denovoflag:
                    break

#output file
print "TOTAL NUMBER OF DOUBLE TRIOS: " + str(doubletrio)

f_summary = sys.argv[5] + ".individual.summary.txt"
summary_output = open(f_summary, 'w+')
summary_output.truncate()
s_output = "TOTAL NUMBER OF MENDEL ERRORS: " + str(mendelerror) + "\n"
summary_output.write(s_output)
s_output = "TOTAL NUMBER OF TRUE DEVONO MUTATIONS: " + str(denovocount) + "\n"
summary_output.write(s_output)
s_output = "TOTAL NUMBER OF AMBIGUOUS DEVONO MUTATIONS: " + str(ambiguous) + "\n"
summary_output.write(s_output)
s_output = "TOTAL NUMBER OF DOUBLE TRIOS: " + str(doubletrio) + "\n"
summary_output.write(s_output)
s_output = "TOTAL NUMBER OF MENDEL ERRORS IN DOUBLE TRIOS: " + str(dtriomendelerror) + "\n"
summary_output.write(s_output)

summary_output.write("INDIVIDUAL SUMMARY:\n")
summary_output.write("IND\tCOUNT\n")

if ind_denovo:
    for key in ind_denovo.keys():
        s_output = key + "\t" + str(ind_denovo[key]) + "\n"
        summary_output.write(s_output)

f_summary2 = sys.argv[5] + ".STR.summary.txt"
summary_output2 = open(f_summary2, 'w+')
summary_output2.truncate()
s_output = "TOTAL NUMBER OF DEVONO MUTATIONS: " + str(denovocount) + "\n"
summary_output2.write(s_output)

summary_output2.write("STR SUMMARY:\n")
summary_output2.write("STR\tCOUNT\n")

if str_denovo:
    for key in str_denovo.keys():
        s_output = key + "\t" + str(str_denovo[key]) + "\n"
        summary_output2.write(s_output)


summary_output.close()
denovo_summary.close()
