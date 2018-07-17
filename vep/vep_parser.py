#parse vep vcf file for variants at each locus
import sys

HIGH = 3
MODERATE = 2
LOW = 1
MODIFIER = 0
def main():
    if len(sys.argv) < 3:
        print "USAGE: python vep_parser.py [input vep vcf] [output file]"
        exit()
    vep = []
    if ".gz" in sys.argv[1]:
        import gzip
        with gzip.open(sys.argv[1]) as f:
            lines=f.readlines()
            for x in lines:
                a = x.strip().split("\t")[:8]
                if "#" not in a[0]:
                    vep.append(a)
    else:
        with open(sys.argv[1]) as f:
            lines=f.readlines()
            for x in lines:
                a = x.strip().split("\t")[:8]
                if "#" not in a[0]:
                    vep.append(a)
    n_output = sys.argv[2] #output file name
    f_output = open(n_output, 'w+')
    f_output.truncate()
    f_output.write("#CHR POS ALLELE IMPACT VARIANT GENE\n")
    for locus in vep:
        l_output = [] #list holds output per locus
        l_output.extend((locus[0], locus[1]))
        variants = locus[7].split(",")
        genelist = {}
        for variant in variants:
            #print variant
            field = variant.split("|")
            if "CSQ=" in field[0]:
                field[0] = field[0].replace("CSQ=", "")
            if field[3] == '':
                field[3] = "Not_Mapped"

            if field[2] == "HIGH":
                    impact = HIGH
            elif field[2] == "MODERATE":
                    impact = MODERATE
            elif field[2] == "LOW":
                    impact = LOW
            else:
                    impact = MODIFIER

            if field[3] not in genelist.keys():
                genelist[field[3]] = [field[0], field[1], field[2], impact]
            else:
                #compare impact level
                if genelist[field[3]][3] < impact:
                    genelist[field[3]] = [field[0], field[1], field[2], impact]
        for gene in genelist.keys():
            l_output.extend((genelist[gene][0], genelist[gene][1], genelist[gene][2], gene))
        f_output.write(" ".join(l_output)+ '\n')

if __name__ == "__main__":
    main()
