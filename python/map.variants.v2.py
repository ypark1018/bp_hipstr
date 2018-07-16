#map the vep annotations to known genetic variants

import sys

def main():
    if len(sys.argv) < 4:
        sys.exit("USAGE: python map.variants.py [file path] [variant file] [gene file]\n")
    annot_file = sys.argv[1]
    variant_file = sys.argv[2] #insert variant file here
    gene_file = sys.argv[3]
    variant = []
    genes = []
    with open(variant_file) as f:
        blengths = f.readlines()
        for allele in blengths:
            b = allele.strip().split("\t")
            variant.append(b)
    with open(gene_file) as e:
        blengths = e.readlines()
        for b in blengths:
            genes.append(b.strip())
    #variant_position = 0 #this will track the current position in the variant file
    #for x in range(1,22):
    if 1:
        STR = []
#        STR_file_suffix = STR_file.format(x)
        with open(annot_file) as g:
            alengths = g.readlines()
            for allele in alengths:
                a = allele.strip().split(" ")
                STR.append(a)
        variant_position_s = 0
        variant_position_e = 0
        for locus in STR:
            if "#" in locus[0]:
                continue
            #variant_position_s = 0
            if locus[5] not in genes:
                continue
            while locus[0] != variant[variant_position_s][0]:
                variant_position_s+=1
            variant_position_e = variant_position_s
            while locus[0] == variant[variant_position_e][0]:
                variant_position_e+=1
            for var in variant[variant_position_s:(variant_position_e+1)]:
                if (locus[1] <= var[2] and locus[1] >= var[1]): #genomic position end
                        #print "test" + "\t" + locus[1] + "\t" + locus[4] + "\t" + "\t".join(variant[variant_position][1:4])
                    if locus[5] in var[3].split(", "):
                        print "\t".join(locus) + "\t" +  var[4]
                        break
                    
if __name__ == "__main__":
    main()
