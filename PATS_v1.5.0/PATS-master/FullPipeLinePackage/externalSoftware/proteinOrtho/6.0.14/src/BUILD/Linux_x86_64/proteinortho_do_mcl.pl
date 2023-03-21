#!/usr/bin/env perl
#pk

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
# proteinortho_do_mcl.pl
# input CORES (number of cores for mcl) and BLASTGRAPH1, ... (the blast-graphs)
# output proteinortho 
#
# @author Paul Klemm
# @email klemmp@staff.uni-marburg.de
# @company Bioinformatics, University of Leipzig
# @version 1
# @date 8-04-2019
#
##########################################################################################


use warnings;
use strict;

if ( !$ARGV[0] or !$ARGV[1] ) {
	print STDERR "Usage: perl proteinortho_do_mcl.pl CORES BLASTGRAPH1 (BLASTGRAPH2, ...)\n\nPerformes the mcl algorithm, given a sequence of blastgraphs (tab 6 format)\n";
	exit;
}

my $cores = $ARGV[0];

sub Error {
	print STDERR "\n[ERROR] ".$_[0]."\n";
	exit 0;
}

my $cmd = "mcl --version";
my $out = qx($cmd);
if ($out =~ /mcl (.+)\n/) {
	print STDERR "Detected 'mcl' version $1\n";
}else{
	&Error("Failed to detect 'mcl'! Tried to call 'mcl'.");
	exit;
}

# the output file of mcl (gets deleted in the end)
my $raw_mcl_output = $ARGV[1];
$raw_mcl_output .= 'mcltmp.proteinortho.mcl'; 

# the .proteinortho-graph format (all blast edges that are in a cluster of mcl)
my $outname_proteinorthograph = $ARGV[1];
$outname_proteinorthograph = 'mcl.proteinortho-graph';

# the .proteinortho formatted output of mcl (with *, speciecounter, genecounter, ...)
my $outname_proteinortho_FINAL = $ARGV[1];
$outname_proteinortho_FINAL = 'mcl.proteinortho';

my %speciesmap;

open(my $OF_proteinortho_FINAL, '>:encoding(UTF-8)',$outname_proteinortho_FINAL) or die "Could not open file $outname_proteinortho_FINAL $!";
	print $OF_proteinortho_FINAL "# Species\tGenes\tAlg.-Conn.";

	#
	# 1. parse the input file for mcl (blast tab format -> 3 coloumns with species-gene to species-gene in the first 2 cols and the bitscore in the third)
	# --> L.faa+++Gene1\tG.faa+++GeneX\t3344
	#

	my $tmpfilename= 'tmp.'.rand(10000000).".mclin";

	open(my $OF_tmpfilename, '>:encoding(UTF-8)', $tmpfilename) or die "Could not open temporary file $tmpfilename $!";

		for (my $i = 1 ; $i < scalar(@ARGV) ; $i=$i+1) {
			print STDERR "Preprocessing ".$ARGV[$i]."\n";

			open(my $IF, '<:encoding(UTF-8)', $ARGV[$i])
			  or die "Could not open file $ARGV[$i] $!";
			my $lastSpeciesA='';
			my $lastSpeciesB='';
			my $lastSpeciesmedianBitscore=255;
			my $offset=2;
			while (my $row = <$IF>) {
				$row =~ s/[\r\n]//g;
				my @spl = split("\t",$row);
				if(length($row)==0 || scalar(@spl)==0){ next; }
				if($offset>0){ $offset=$offset-1; next; }
				if(substr($row,0,1) eq "#" && scalar(@spl)==2){ $spl[0] =~ s/\# //g; $lastSpeciesA=$spl[0]; $lastSpeciesB=$spl[1]; next; }
				if(substr($row,0,1) eq "#" && scalar(@spl)==4){ $spl[0] =~ s/\# //g; $lastSpeciesmedianBitscore=($spl[1]+$spl[3])/2; next; }
				if(scalar(@spl)<5){ next; }

				if(!exists $speciesmap{$lastSpeciesA}){
					$speciesmap{$lastSpeciesA}=scalar(keys(%speciesmap));
					print $OF_proteinortho_FINAL "\t".$lastSpeciesA;
				}
				if(!exists $speciesmap{$lastSpeciesB}){
					$speciesmap{$lastSpeciesB}=scalar(keys(%speciesmap));
					print $OF_proteinortho_FINAL "\t".$lastSpeciesB;
				}

				print $OF_tmpfilename $lastSpeciesA."+++".$spl[0]."\t".$lastSpeciesB."+++".$spl[1]."\t".((255/$lastSpeciesmedianBitscore)*$spl[3])."\n";
				# here the species name is encoded in the genename (speciesA+++geneA1)
			}
			close($IF);
		}
		print $OF_proteinortho_FINAL "\n";

	close($OF_tmpfilename);
close($OF_proteinortho_FINAL);


#
# 2. run mcl
#

print STDERR "Running mcl ... \n";

system('mcl '.$tmpfilename.' --abc -te '.$cores.' -o '.$raw_mcl_output.' 2>/dev/null');

unlink($tmpfilename); # remove input tmp file

#
# 3. Format the output of mcl -> .proteinortho and .proteinortho-graph format
#

print STDERR "Postprocessing mcl file ... \n";


my %is_in_mclcluster;

#
# 3.1 proteinortho format
#
open(my $IF_raw_mcl_output, '<:encoding(UTF-8)',$raw_mcl_output) or die "Could not open file $raw_mcl_output $!";
	open($OF_proteinortho_FINAL, '>>:encoding(UTF-8)',$outname_proteinortho_FINAL) or die "Could not open file $outname_proteinortho_FINAL $!";

		# go over the raw output of mcl 
		while (my $row = <$IF_raw_mcl_output>) {
			$row =~ s/[\r\n]//g;
			my @spl = split("\t",$row);
			my @outvector = ("") x (scalar(keys(%speciesmap))+3); # this array contains all genes of a given cluster for each species # +3 because the first 3 cols are speciescounter,genecounter and alg.conn.
			my $genecounter=0;
			my %speciescounter;

			# iterate over all pairs of the given cluster -> the edges for proteinortho-graph file (stored in is_in_mclcluster, lexico ordered key (species and gene id))
			for (my $i = 0 ; $i < scalar(@spl) ; $i=$i+1) {
				for (my $j = $i+1 ; $j < scalar(@spl) ; $j=$j+1) {
					my $key = $spl[$i]."-".$spl[$j];
					if($spl[$j] lt $spl[$i]){
						$key = $spl[$j]."-".$spl[$i];
					}
					# generate a key in the form speciesA+++geneA1-speciesB+++geneB1
					# next: order species+++gene lexicographically such that not both speciesA+++geneA1-speciesB+++geneB1 and speciesB+++geneB1-speciesA+++geneA1 
					# have to be stored (=half the work space)
					# HERE: only the key is generated and the species are allready encoded in the gene name
					
					$is_in_mclcluster{$key}=1;
				}

				#the following is for the proteinortho format (cluster format)
				my @spl2 = split(/\+\+\+/,$spl[$i]);
				if($outvector[$speciesmap{$spl2[0]}+3] ne ""){ $outvector[$speciesmap{$spl2[0]}+3].=","; } 
				$outvector[$speciesmap{$spl2[0]}+3].=$spl2[1]; # +3 because the first 3 cols are speciescounter,genecounter and alg.conn.
				$genecounter++;
				$speciescounter{$spl2[0]}=1;
			}

			for (my $i = 0 ; $i < scalar(@outvector) ; $i=$i+1) {
				if($outvector[$i] eq ""){
					$outvector[$i]="*";
				}
			}

			$outvector[0]=scalar(keys(%speciescounter)); # == the number of species in this row
			$outvector[1]=$genecounter; # == the number of genes in this row
			$outvector[2]="x"; # == no algebraic connectivity
			print $OF_proteinortho_FINAL join("\t",@outvector)."\n"; # join the rest and print
		}
	close($OF_proteinortho_FINAL);
close($IF_raw_mcl_output);

unlink($raw_mcl_output); # remove output tmp file

#
# 3.2 proteinortho-graph format
#
open(my $OF_proteinorthograph, '>:encoding(UTF-8)',$outname_proteinorthograph) or die "Could not open file $outname_proteinorthograph $!";

	print $OF_proteinorthograph "# file_a	file_b\n# a	b	evalue_ab	bitscore_ab	evalue_ba	bitscore_ba\n";

	# now go over the blast graphs and print all rows that are in a cluster (is_in_mclcluster)
	for (my $i = 1 ; $i < scalar(@ARGV) ; $i=$i+1) {
		print STDERR "Postprocessing ".$ARGV[$i]."\n";

		open(my $IF, '<:encoding(UTF-8)', $ARGV[$i])
		  or die "Could not open file $ARGV[$i] $!";
		my $lastSpeciesA='';
		my $lastSpeciesB='';
		my $offset=2;
		while (my $row = <$IF>) {
			$row =~ s/[\r\n]//g; #better than chomp
			my @spl = split("\t",$row);
			if(length($row)==0 || scalar(@spl)==0){ next; }
			if($offset>0){ $offset=$offset-1; next; }
			if(substr($row,0,1) eq "#" && scalar(@spl)==2){ print $OF_proteinorthograph $row."\n"; $spl[0] =~ s/\# //g; $lastSpeciesA=$spl[0]; $lastSpeciesB=$spl[1]; next; }
			if(scalar(@spl)<5){ next; }
			
			my $key = $lastSpeciesA."+++".$spl[0]."-".$lastSpeciesB."+++".$spl[1];
			if($lastSpeciesB."+++".$spl[1] lt $lastSpeciesA."+++".$spl[0] ){
				$key = $lastSpeciesB."+++".$spl[1]."-".$lastSpeciesA."+++".$spl[0];
			}
			# test if the given edge of the blast-graph is part of the cluster (is_in_mclcluster)

			if(exists $is_in_mclcluster{$key}){
				print $OF_proteinorthograph $row."\n";
			}
		}
		
		close($IF);
	}
close($OF_proteinorthograph);
