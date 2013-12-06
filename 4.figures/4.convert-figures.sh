#!/bin/sh

epstopdf img/snps-pvalue.ps --outfile /tmp/snps-pvalue-rotated.pdf
pdf270 /tmp/snps-pvalue-rotated.pdf --outfile img/snps-pvalue.pdf
convert img/snps-pvalue.pdf img/snps-pvalue.png
