#!/usr/bin/env perl
#pk 

use strict;
use warnings "all";

if(!defined($ARGV[0]) || $ARGV[0] eq "-h" || $ARGV[0] eq "--help" || $ARGV[0] eq "help" || $ARGV[0] eq "?" || $ARGV[0] eq "h"){
	print STDERR "proteinortho2xml.pl PROTEINORTHOFILE\n";
	print STDERR "Reads Proteinortho file (not proteinortho-graph file!) and produces the OrthoXML format (>stdout).\n\n";
	exit;
}

my $po_file = $ARGV[0];

my $cur_id = 0;

my @species;
my @species_num_prots;
my $headerisset=0;
my %protID2id;
my $orthologygroup="";

open my $fh, $po_file or die "Could not open $po_file: $!";

my $nextIsHeaderLine = 0;
my %species_graph;
my $species_a="";
my $species_b="";

if($po_file=~m/proteinortho-graph/){
	while( my $line = <$fh>)  {  # each line = orthology GRAPH !

		if($line =~ m/(# file_a\tfile_b)|(# a\tb\t)/){
			next
		}elsif($line =~ /^# /){

			chomp($line);
			$line=~s/^# ?//g;
	   		my @linesplt=split(/\t/,$line);
	   		
	   		if(scalar @linesplt != 2){next;}

	   		$species_a = $linesplt[0];
	   		$species_b = $linesplt[1];

			my $NCBITaxId_a=$species_a;
			my $name_a = $species_a;
			if($NCBITaxId_a=~m/UP[0-9]+\_([0-9]+)/){
				$NCBITaxId_a=$1;
			}elsif($NCBITaxId_a=~m/([0-9]+)\_([^.]+).*/){
				$NCBITaxId_a=$1;
				$name_a=$2;
			}else{
				$NCBITaxId_a='-1';
			}
			if($name_a=~m/([^\.]+)\.(f[na]a?|fasta|f)/){
				$name_a=$1;
			}
			my $NCBITaxId_b=$species_b;
			my $name_b = $species_b;
			if($NCBITaxId_b=~m/UP[0-9]+\_([0-9]+)/){
				$NCBITaxId_b=$1;
			}elsif($NCBITaxId_b=~m/([0-9]+)\_([^.]+).*/){
				$NCBITaxId_b=$1;
				$name_b=$2;
			}else{
				$NCBITaxId_b='-1';
			}
			if($name_b=~m/([^\.]+)\.(f[na]a?|fasta|f)/){
				$name_b=$1;
			}

	   		if(!exists $species_graph{$linesplt[0]}){
				$species_graph{$linesplt[0]}{"header"} = '<species name="'.$name_a.'" NCBITaxId="'.$NCBITaxId_a.'">'."\n\t\t".'<database name="'.$name_a.'" version="'.$species_a.'">'."\n\t\t\t<genes>\n"; # save the species
	   		}
	   		if(!exists $species_graph{$linesplt[1]}){
				$species_graph{$linesplt[1]}{"header"} = '<species name="'.$name_b.'" NCBITaxId="'.$NCBITaxId_b.'">'."\n\t\t".'<database name="'.$name_b.'" version="'.$species_b.'">'."\n\t\t\t<genes>\n"; # save the species
	   		}
			
			$nextIsHeaderLine = 0;
		}elsif($line !~ m/^#/){

	   		my @linesplt=split(/\t/,$line); # -> first 3 cols are number of species/genes and the algebraic connectivity of the given group

	   		if(scalar @linesplt < 3){next;}

	   		my $prot_a = $linesplt[0];
	   		my $prot_b = $linesplt[1];

			if($prot_a=~m/^[^\|]+\|([^\|]+)\|[^\|]+/){ # extract the true name (A8JEJ4)
				$prot_a=$1;
			}
			if($prot_b=~m/^[^\|]+\|([^\|]+)\|[^\|]+/){ # extract the true name (A8JEJ4)
				$prot_b=$1;
			}

			if(!exists $species_graph{$species_a}{$linesplt[0]}) { $species_graph{$species_a}{$linesplt[0]} = 1; }
			if(!exists $species_graph{$species_b}{$linesplt[1]}) { $species_graph{$species_b}{$linesplt[1]} = 1; }

			if(!exists($protID2id{$species_a."#".$linesplt[0]})){ # map down to a integer A8JEJ4->0 ...
				$protID2id{$species_a."#".$linesplt[0]}=scalar(keys(%protID2id));
			}
			if(!exists($protID2id{$species_b."#".$linesplt[1]})){ # map down to a integer A8JEJ4->0 ...
				$protID2id{$species_b."#".$linesplt[1]}=scalar(keys(%protID2id));
			}

			my $curorthologygroup="\t\t\t".'<geneRef id="'.$protID2id{$species_a."#".$linesplt[0]}."\"/>\n"."\t\t\t".'<geneRef id="'.$protID2id{$species_b."#".$linesplt[1]}."\"/>\n";

			$orthologygroup.="\t\t<orthologGroup id=\"$cur_id\">\n".$curorthologygroup."\t\t</orthologGroup>\n";
			$cur_id++;
		}

	}
}else{
	while( my $line = <$fh>)  {  # each line = orthology group
	    if(length($line)>0 && substr($line,0,1) ne "#"){ 
	    	chomp($line);

	    	if(!$headerisset){print STDERR "[ERROR] header is missing.\n";die}

	   		my $curorthologygroup = ""; # create template orthology group

	   		my @linesplt=split(/\t/,$line); # -> first 3 cols are number of species/genes and the algebraic connectivity of the given group
	   		my $con="";

			for(my $i = 0 ; $i < scalar(@linesplt) ; $i=$i+1){ # iterate over each column
	   			if( $headerisset == 0 ){

					$curorthologygroup.='<geneRef id="'.$linesplt[$i]."\"/>\n"; # for the .mcl format 

				}elsif( ( $i > 2 ) && $linesplt[$i] ne "*" && $linesplt[$i] ne ""){

	   				my @linesplt2=split(/,/,$linesplt[$i]); # split again for multiple genes of one species (in the current orthology group)

					for(my $j = 0 ; $j < scalar(@linesplt2) ; $j=$j+1){

						if( $linesplt2[$j] eq "*" ){next}
						# $linesplt2[$j] is a protein (name like tr|A8JEJ4|A8JEJ4_CHLRE) of the current ortho group

						if(!exists($protID2id{$i."#".$linesplt2[$j]})){ # map down to a integer A8JEJ4->0 ...
							$protID2id{$i."#".$linesplt2[$j]}=scalar(keys(%protID2id));
							# $i."#".$linesplt2[$j] = the species id (column number i) and the gene name = identifier
						}
						# add to the species and ortho group

						my $protId=$linesplt2[$j];
						$protId=~s/(UniProtKB|Swiss-Prot|TrEMBL|ENSEMBL|\/)+://g; # protein names can be 'UniProtKB/Swiss-Prot:712835' -> '712835'
						$protId=~s/^[sptr]{2}\|//g; # protein names can be 'sp|O08314|TGT_HELPY' -> 'O08314'
						$protId=~s/\|[^|]+$//g; # protein names can be 'sp|O08314|TGT_HELPY' -> 'O08314'

						if(scalar @species_num_prots != scalar @linesplt){ @species_num_prots=(0)x(scalar @linesplt) }
						$species_num_prots[$i-3]++;
						$species[$i-3].="\t\t\t\t".'<gene id="'.$protID2id{$i."#".$linesplt2[$j]}.'" protId="'.$protId."\"/>\n";
						$curorthologygroup.="\t\t\t".'<geneRef id="'.$protID2id{$i."#".$linesplt2[$j]}."\"/>\n";
					}

	   			}elsif($i eq 2){
	   				$con=$linesplt[$i]; # 3rd col is connectivity
	   			}
	   		}
	   		if($orthologygroup ne ""){$orthologygroup.="\n";}
	   		if($con eq "NA"){$con="-1"}
	   		$orthologygroup.="\t\t<orthologGroup id=\"$cur_id\">\n\t\t\t<score id=\"algcon\" value=\"$con\"/>\n".$curorthologygroup."\t\t</orthologGroup>";

			$cur_id=$cur_id+1;

	    }elsif(substr($line,0,1) eq "#"){ # first line (with comment #) -> contains the names of the files/species 

	    	$headerisset=1;

			chomp($line);
	   		my @linesplt=split(/\t/,$line);

			for(my $i = 0 ; $i < scalar(@linesplt) ; $i=$i+1){
				if($i>2 && length($linesplt[$i])>0){
					my $NCBITaxId=$linesplt[$i];
					my $name = $linesplt[$i];
					if($NCBITaxId=~m/UP[0-9]+\_([0-9]+)/){
						$NCBITaxId=$1;
					}elsif($NCBITaxId=~m/([0-9]+)\_([^.]+).*/){
						$NCBITaxId=$1;
						$name=$2;
					}else{
						$NCBITaxId='-1';
					}
					if($name=~m/([^\.]+)\.(f[na]a?|fasta|f)/){
						$name=$1;
					}

					push(@species,"\t".'<species name="'.$name.'" NCBITaxId="'.$NCBITaxId.'">'."\n\t\t".'<database name="'.$name.'" version="'.$linesplt[$i].'">'."\n\t\t\t<genes>\n"); # save the species
				}
			}
	    }
	}
}
close $fh;

print '<?xml version="1.0" encoding="utf-8"?>';
print "\n".'<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="proteinortho" originVersion="6" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd">';
#"\n<orthoXML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"0.3\" origin=\"proteinortho\" originVersion=\"6.0\" xsi:schemaLocation=\"http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd\">";
print "\n\t<notes>\n\t\tProteinortho OrthoXML file.\n\t</notes>\n";

if($po_file=~m/proteinortho-graph/){
	foreach my $cur_species (keys %species_graph){

		if(scalar keys %{$species_graph{$cur_species}} > 0){

			print $species_graph{$cur_species}{"header"};

			foreach my $prot (keys %{$species_graph{$cur_species}}){
				if($prot eq "header"){next;}

				my $prot_id = $prot;

				if($prot=~m/^[^\|]+\|([^\|]+)\|[^\|]+/){ # extract the true name (A8JEJ4)
					$prot_id=$1;
				}

				print "\t\t\t\t".'<gene id="'.$protID2id{$cur_species."#".$prot}.'" protId="'.$prot_id."\"/>\n";
			}

			print "\t\t\t</genes>\n\t\t</database>\n\t</species>\n";
		}
	}
}else{
	for (my $i = 0; $i < scalar @species; $i++) {
		if($species_num_prots[$i]>0){
			print $species[$i]."\t\t\t</genes>\n\t\t</database>\n\t</species>\n";
		}
	}
}

print "\t".'<scores>'."\n\t\t".'<scoreDef id="algcon" desc="proteinortho score of connectivity (normalized algebraic connectivity). The higher the value the better the component is connected."/>'."\n\t".'</scores>'."\n";

print "\t<groups>\n".$orthologygroup."\n\t</groups>\n";

print '</orthoXML>';