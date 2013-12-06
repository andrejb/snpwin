#!/usr/bin/python

import math

# read windows results
f = open('../results/windows-windows.txt')
windows = f.readlines()
f.close()

# read snp names and cromosomes
f = open('/var/tmp/snp-cromo-pos.txt', 'r')
cromopos = f.readlines()
f.close()

snppos = []
snpcromo = []
snpname = []
for l in cromopos:
    cromo, pos, name = l.strip().split('\t')
    snpcromo.append(int(cromo))
    snppos.append(int(pos))
    snpname.append(name)

f = open("../4.figures/tmp/pvalue.txt", "r")
i = 0
d = []
for spval in f:
    pval = float(spval.strip())
    if pval >= 4.0:
        start, end = windows[i].strip().split(' ')
        start = int(start)
        end = int(end)
        if math.isinf(pval):
            pval = 16.0
        d.append( (snpname[start], snpname[end], (end-start+1), (snppos[end]-snppos[start]+1), snpcromo[start], pval, start, end) )
    i += 1
f.close()

for l in sorted(d, key=lambda vals: vals[5], reverse=True):
    if l[5] == 16.0:
      pval = "<10^{-16}"
      formatstr  = "%s %s %d %d %d %s"
    else:
      pval = 10**(-l[5])
      formatstr = "%s %s %d %d %d %e"
    print formatstr  % (l[0], l[1], l[2], l[3], l[4], pval)
