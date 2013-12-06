run: do-figures

do-figures:
	./1.calc-p-value.m
	./2.gen-cromos.pl
	./3.gen-figure.gnuplot
	./4.convert-figures.sh

clean:
	rm -f img/* tmp/*
