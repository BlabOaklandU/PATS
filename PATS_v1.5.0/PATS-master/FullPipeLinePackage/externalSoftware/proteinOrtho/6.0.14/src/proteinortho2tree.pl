#!/usr/bin/env perl

##########################################################################################
#	  This file is part of proteinortho.
#	  (C) 2009 Marcus Lechner
# 
#	  proteinortho is free software; you can redistribute it and/or modify
#	  it under the terms of the GNU General Public License as published
#	  by the Free Software Foundation; either version 2, or (at your
#	  option) any later version.
#
#	  proteinortho is distributed in the hope that it will be useful, but
#	  WITHOUT ANY WARRANTY; without even the implied warranty of
#	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#	  General Public License for more details.
#
#	  You should have received a copy of the GNU General Public License
#	  along with proteinortho; see the file COPYING.  If not, write to the
#	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#	  Boston, MA 02111-1307, USA.	
##########################################################################################

##########################################################################################
# About
##########################################################################################
# orthomatrix2tree
# input ortheset output from proteinortho
# output corresponding UPGMA tree in newick format
#
# branch labels show the number of common genes in the subtrees
# branch lengths represent the number of new common genes since the last node
# 
# @author Marcus Lechner
# @email lechner@staff.uni-marburg.de
# @company Bioinformatics, University of Leipzig
# @version 3.10
# @date 16-03-2016
#
##########################################################################################

##########################################################################################
# Imports
##########################################################################################
use strict;
use warnings "all";
use File::Basename;


##########################################################################################
# Usage, defaults and parameters
##########################################################################################
my $usage= << "JUS";
  usage:   proteinortho2tree.pl [OPTIONS] ORTHOMATRIX(.tsv) >OUTTREE
  input:   output from Proteinortho (version >5) e.g. myproject.proteinortho.tsv (in proteinortho5 id does not have the tsv suffix)
  output:  corresponding UPGMA tree in newick format
  options: -o=[FILE]  prints output to the given file rather than STDOUT  
JUS

# path to the C-part of this program
my @tmppath = fileparse($0);
our $scriptpath = $tmppath[1];

# Output switch
my $o = "";
my $path = "";
my $file = "";

# No parameters
unless (defined($ARGV[0])) {
	print "Error: No input file defined!\n\n$usage";
	exit 1;
}

foreach (@ARGV) {
	if 	($_ =~ /^--?o.*=(.+)/) {$o = $1;}						# output file
	elsif	($_ =~ /^-/) {	print STDERR "Unknown parameter $_\n$usage"; exit 1;}
	elsif	(!-e $_) {	print STDERR "Could not find file $_\n"; exit 1;}
	else 	{$file = $_;}
}

# open the given file, this should be a proteinortho matrix
open(SOURCE,"<$file") || die("Could not open file $file\n$!");
# write temp file
open(TRG,">$file.tmp.matrix2") || die("Could not open tmp-file $file.tmp.matrix2, $!");
my $did = 0;
my @viecher = ();

my @present_in_groups;
# convert the orthoset to a 1/0 matrix
print STDERR "Preparing datafile...\n";
while(<SOURCE>) {
	my @data = split;
	# first line
	if ($data[0] =~ /^\#/) {
		# remove the first elements
		shift(@data);shift(@data);shift(@data);

		# first line
		if ($did == 0) {
			shift(@data);
			# species names to header
			@viecher = @data;
		}
	
	}
	else {
		shift(@data);shift(@data);shift(@data);
		# any futher line
		my $i = 0;
		foreach (@data) {
			# just transform the entries to 0 and 1
			# the gene-names are not important here
			if ($_ eq "*") {
				print TRG 0;
			}
			else {
				# Count number of copies
				my $num = 1;
#				$num += $_ =~ /,/g;		# count number of copies (not yet respected in C)
				print TRG $num;
				$present_in_groups[$i]++; # count group
			}
			print TRG " ";
			$i++;
		}
		print TRG "\n";
	}
}
close(SOURCE);
close(TRG);

open(TMP,">$file.tmp.matrix") || die("Could not open tmp-file $file.tmp.matrix, $!");
	for (my $i = 0; $i < scalar(@viecher); $i++) {
		$viecher[$i] =~ s/\..*?$//; 			# remove point
		$viecher[$i] =~ s/\s/_/g; 			# alphabet
		$viecher[$i] .= "_[$present_in_groups[$i]]";	# add number of groups
	}
	print TMP scalar(@viecher)."\n";
	print TMP join(" ",@viecher)."\n";
close(TMP);
system("cat '$file.tmp.matrix2' >>'$file.tmp.matrix'");

# Run the main algorithm in C
print STDERR "Calculating tree...\n";
my $run = "proteinortho_treeBuilderCore";
my $out = qx($run $ARGV[0].tmp.matrix);
$out =~ s/\[.+?\]//g;
$out =~ s/\s:/:/g;

if ($o ne "") {
	open(OUT,">$o") || die("Could not write output file $!");
	print OUT $out;
	close(OUT);
}
else {
	print $out;
}

# unlink temporary files
unlink("$file.tmp.matrix");
unlink("$file.tmp.matrix2");
