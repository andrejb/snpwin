#!/usr/bin/python

import re

windowsf = '../results/windows-windows.txt'
posf = '/var/tmp/snp-pos.txt'
resultf = '../results/wlen-x-distance.txt'

# read positions from file
position = []
f = open(posf, 'r')
for line in f:
  position.append(int(line))
f.close()

# regular expression to split input
m = re.compile(' ')

# open output file
o = open(resultf, 'w')

# for each window, calculate distances
f = open(windowsf, 'r')
for line in f:
  s, e = m.split(line)
  start = int(s)
  end = int(e)
  diff = end - start + 1
  diff2 = position[end] - position[start] + 1
  o.write("%d %d\n" % (diff, diff2))
f.close()

# close output file
o.close()
