#!/usr/bin/env perl

use warnings;
use strict;

unless (defined($ARGV[0])) {
	print STDERR "Usage: proteinortho_extract_from_graph.pl PROTEINORTHO_TABLE <GRAPH\n\nExtracts an orthology group from a given graph (STDIN)\n";
	exit;
}

my %use_gene;
my @species;
open(PO,"<$ARGV[0]") || die($!);
while (<PO>) {
	chomp;
	my @list = split(/\t+/);
	if ($list[0] eq "# Species") {
		@species = @list;
		shift @species;
		shift @species;
		shift @species;
	}
	else {
		shift @list;
		shift @list;
		shift @list;
		for (my $i = 0; $i < scalar(@list); $i++) {
			foreach my $gene (split(",",$list[$i])) {
				my $id = $gene.' '.$species[$i];
				$use_gene{$id}++;
			}
		}
	}
}
close(PO);

print STDERR scalar(keys %use_gene);
print STDERR " genes found\n";

my $file_a;
my $file_b;
while (<STDIN>) {
	if ($_ =~ /^#/) {
		print $_;
		my @row = split(/\s+/);
		unless ($row[1]) {die();}
		unless ($row[1] eq "file_a" || $row[1] eq "a") {
			$file_a = $row[1];
			$file_b = $row[2];
		}
	}
	else {
		my @row = split(/\s+/);
		unless ($row[1]) {next;}
		my $ida = $row[0].' '.$file_a;
		my $idb = $row[1].' '.$file_b;
		if ($use_gene{$ida} || $use_gene{$idb}) {
			print $_;
		}
	}
}
