#!/usr/bin/gnuplot

set term postscript
set output './img/wlen-x-snppos.ps'

set title 'Window length versus SNP distance'

set xlabel 'window length'
set ylabel 'SNPs distance'

set size 1,0.5

set yrange [0:1000000]

plot '../results/wlen-x-distance.txt' notitle with points pt 7 lc 0 ps 0.2
