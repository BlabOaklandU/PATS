#!/usr/bin/env perl

use warnings;
use strict;

unless ($ARGV[0]) {
	print STDERR "Usage: formatU.pl FILE_A\n\nExpecting usearch/ublast '-blast6out' format file FILE_A. Removing the description from each gene name (starting with a whitespace character).\n";
	exit;
}

open(FILE,"<$ARGV[0]") || die("Error, could not open file $ARGV[0]: $!");
while(<FILE>) {
	chomp;
	if ($_ =~ /^#/) {next;}
	$_ =~ s/ {1}[^\t]+\t/\t/;
	$_ =~ s/ {1}[^\t]+\t/\t/;
	print STDOUT $_."\n";
}
close(FILE);
