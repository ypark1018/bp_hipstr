#change name of item given table with conversion
import sys

if 

def main():
    if len(sys.argv) < 3:
        print "USAGE: python convert_name.py [conversion table] [input vcf]"
        exit()
    names = {} #dictionary with conversion
    vcf = [] #vcf file for conversion
    with open(sys.argv[1]) as f:
        line = f.readlines()
        for entry in line:
            i = entry.strip().split(" ")
            names[i[1]] = i[2]
    if ".gz" in sys.argv[1]:
        import gzip
        with gzip.open(sys.argv[2]) as g:
            line = g.readlines()
            for entry in line:
                a = entry.strip().split("\t")[:8]
    else:
        with open(sys.argv[1]) as g:
            line = g.readlines()
            for entry in line:
                vcf.append(entry.strip().split("\t"))
            
