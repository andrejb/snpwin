#!/home/drebs/tmp/octave/bin/octave

x = load("../results/chi-square-chi-square.txt");
y = -log10(1-chi2cdf(x(:,1),x(:,2)));
save "-ascii" "tmp/pvalue.txt" y;
