import sys

chrnum = sys.argv[1]
markers = sys.argv[2]
with open(chrnum) as f:
    chr1 = f.readlines()
with open(markers) as g:
    STR = g.readlines()

#clean
for n,i in enumerate(chr1):
    chr1[n] = i.rstrip()
for n,i in enumerate(STR):
    STR[n] = i.rstrip()

STRlist = []
for x in STR:
    sp = x.split()
    STRlist.append(sp)

for n,x in enumerate(chr1):
    minChr = 251
    STRindex = "NAN"
    for y in STRlist:
        diff = abs(int(y[1]) - int(x))
        if diff < minChr:
            minChr = diff
            STRindex = y[0]
    if minChr < 251:
        chr1[n] = [chr1[n], STRindex, minChr]
        m = n - 5
        while m < n:
            if chr1[m][1] == chr1[n][1]:
                if chr1[m][2] < chr1[n][2]:
                    chr1[n][1] = "NAN"
                else:
                    chr1[m][1] = "NAN"
            m += 1
    else:
        chr1[n] = [chr1[n], "NAN", minChr]

for x in chr1:
   print x[0] + "\t" + x[1]
