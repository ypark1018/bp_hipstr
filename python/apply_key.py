import sys

allele_ref = sys.argv[1]
chr_key = sys.argv[2]
with open(allele_ref) as f:
    chr1 = f.readlines()

with open(chr_key) as g:
    key = g.readlines()

chr1list = []
for x in chr1:
    s = x.strip().split("\t")
    s[1] = s[1].replace("chr1:", "")
    chr1list.append(s)

keylist = []
for y in key:
    s = y.strip().split("\t")
    keylist.append(s)

key_dict = {}
#make dictionary
for x in keylist:
    key_dict[x[0]] = x[1]

#apply key
newlist = []
for n,x in enumerate(chr1list):
    chr1list[n][1] = key_dict[x[1]]

for x in chr1list:
    print '\t'.join(x)
