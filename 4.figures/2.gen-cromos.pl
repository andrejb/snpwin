#!/usr/bin/perl

open(CROMO, "res/cromo-windows.txt"); 
$cromo{0} = 0;
while (<CROMO>) {
  $_ =~ m/([0-9]+) ([0-9]+)/;
  $cromos{ $1 } = $2;
}
close(CROMO);

$cromos{22} = 48698;

$c = 1;
$i = 1;
$ind = 1;

open(SNPS, "tmp/pvalue.txt");

open(CROMO, ">tmp/cromo-$c.txt");

open (SPECIAL, ">tmp/special.txt");

$last = 0;
while (<SNPS>) {
  if ($i > $cromos{$c}) {
    close(CROMO);
    print (($ind+$last)/2);
    print "\n";
    $c++;
    $ind+=500;
    $last = $ind;
    open(CROMO, ">tmp/cromo-$c.txt");
  }
  if ($_ =~ m/([0-9]+\.[0-9]+.*)/) {
    print CROMO "$ind $1\n";
  }
  else {
    print CROMO "$ind -10\n";
    print SPECIAL "$ind 16\n";
  }
#  print "$ind $1\n";
  $i++;
  $ind++;
}

close(CROMO);

close(SPECIAL);

close(SNPS);
