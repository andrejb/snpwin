#!/usr/bin/gnuplot

set term postscript
set output './img/snps-pvalue.ps'

set title 'Chi-square statistics for case and control SNP windows'

set xlabel 'chromosome'
set ylabel '-log10(p-value)'

set size 1,0.5

set xtics ( \
"1" 1887,  \
"2" 6275,  \
"3" 10452, \
"4" 14158, \
"5" 17710, \
"6" 21343, \
"7" 24788, \
"8" 28012, \
"9" 31087, \
"10" 33971, \
"11" 36875, \
"12" 39806, \
"13" 42449, \
"14" 44704, \
"15" 46796, \
"16" 48811, \
"17" 50739, \
"18" 52660, \
"19" 54377, \
"20" 55991, \
"21" 57529, \
"22" 59000 ) font "Helvetica,10"
set ytics font "Helvetica,10"

set xrange [-1000:60000]
set yrange [0:17]

plot 'tmp/cromo-1.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-2.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-3.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-4.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-5.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-6.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-7.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-8.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-9.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-10.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-11.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-12.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-13.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-14.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-15.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-16.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-17.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-18.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-19.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-20.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/cromo-21.txt' notitle with points pt 7 lc 27 ps 0.4, \
    'tmp/cromo-22.txt' notitle with points pt 7 lc 0 ps 0.4, \
    'tmp/special.txt' notitle with points pt 7 lc 0 ps 0.4
