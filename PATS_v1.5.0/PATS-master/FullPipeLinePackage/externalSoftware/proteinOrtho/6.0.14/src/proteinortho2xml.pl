#!/usr/bin/env perl

use strict;
use warnings "all";

if(!defined($ARGV[0]) || $ARGV[0] eq "-h" || $ARGV[0] eq "--help" || $ARGV[0] eq "help" || $ARGV[0] eq "?" || $ARGV[0] eq "h"){
	print STDERR "proteinortho2xml.pl PROTEINORTHOFILE\n";
	print STDERR "Reads Proteinortho file (not proteinortho-graph file!) and produces the OrthoXML format (>stdout).\n\n";
	exit;
}

my $po_file = $ARGV[0];

my $curOrthoXML = {
	'orthoXML'=>[
		{
			'version'=>'0.3',
			'xmlns'=>'http://orthoXML.org/2011/',
			'origin'=>'proteinortho',
	  		'originVersion'=>'5',
			'xmlns:xsi'=>'http://www.w3.org/2001/XMLSchema-instance',
			'xsi:schemaLocation'=>'http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd',
			'notes' => [{
				'content' => ' Proteinortho OrthoXML file. '
             }],
			'scores' => [{
	            'scoreDef' => [
					{
						'desc' => 'The algebraic connectivity of the given group',
						'id' => 'connectivity'
					}
                ]
			}],
			'species' => [],
			'groups' => [{
                'orthologGroup' => []
            }]
		}	
	]
};

my $cur_id = 0;

my @species;
my $headerisset=0;
my %protID2id;
my $orthologygroup="";

open my $fh, $po_file or die "Could not open $po_file: $!";
while( my $line = <$fh>)  {  # each line = orthology group
    if(length($line)>0 && substr($line,0,1) ne "#"){ 
    	chomp($line);

   		my $curorthologygroup = ""; # create template orthology group

   		my @linesplt=split(/\t/,$line); # -> first 3 cols are number of species/genes and the algebraic connectivity of the given group
   		my $con="";

		for(my $i = 0 ; $i < scalar(@linesplt) ; $i=$i+1){ # iterate over each column
   			if( $headerisset == 0 ){

				$curorthologygroup.='<geneRef id="'.$linesplt[$i]."\"/>\n"; # for the .mcl format 

			}elsif( ( $i > 2 ) && $linesplt[$i] ne "*" && $linesplt[$i] ne ""){

   				my @linesplt2=split(/,/,$linesplt[$i]); # split again for multiple genes of one species (in the current orthology group)

				for(my $j = 0 ; $j < scalar(@linesplt2) ; $j=$j+1){

					# $linesplt2[$j] is a protein (name like tr|A8JEJ4|A8JEJ4_CHLRE) of the current ortho group

					if($linesplt2[$j]=~m/^[\|]+\|([\|]+)\|[\|]+/){ # extract the true name (A8JEJ4)
						$linesplt2[$j]=$1;
					}
					if(!exists($protID2id{$i."#".$linesplt2[$j]})){ # map down to a integer A8JEJ4->0 ...
						$protID2id{$i."#".$linesplt2[$j]}=scalar(keys(%protID2id));
						# $i."#".$linesplt2[$j] = the species id (column number i) and the gene name = identifier
					}
					# add to the species and ortho group

					my $protId=$linesplt2[$j];
					$protId=~s/(UniProtKB|Swiss-Prot|TrEMBL|ENSEMBL|\/)+://g; # protein names can be 'UniProtKB/Swiss-Prot:712835' -> '712835'
					$protId=~s/^[sptr]{2}\|//g; # protein names can be 'sp|O08314|TGT_HELPY' -> 'O08314'
					$protId=~s/\|[^|]+$//g; # protein names can be 'sp|O08314|TGT_HELPY' -> 'O08314'

					$species[$i-3].='<gene id="'.$protID2id{$i."#".$linesplt2[$j]}.'" protId="'.$protId."\"/>\n";
					$curorthologygroup.='<geneRef id="'.$protID2id{$i."#".$linesplt2[$j]}."\"/>\n";
				}

   			}elsif($i eq 2){
   				$con=$linesplt[$i]; # 3rd col is connectivity
   			}
   		}
   		if($orthologygroup ne ""){$orthologygroup.="\n";}
   		$orthologygroup.="<orthologGroup id=\"$cur_id\">\n<score id=\"bit\" value=\"$con\"/>\n".$curorthologygroup."</orthologGroup>";

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

					push(@species,'<species name="'.$name.'" NCBITaxId="'.$NCBITaxId.'">'."\n".'<database name="'.$name.'" version="'.$linesplt[$i].'">'."\n<genes>\n"); # save the species
				}
			}
    }
}
close $fh;

print '<?xml version="1.0" encoding="utf-8"?>';
print "\n<orthoXML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"0.3\" origin=\"proteinortho\" originVersion=\"6.0\" xsi:schemaLocation=\"http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd\">";
print "\n<notes>Proteinortho OrthoXML file.</notes>\n";

foreach my $cur_species (@species){
	print $cur_species."</genes>\n</database>\n</species>\n";
}

print "<groups>\n".$orthologygroup."\n</groups>\n";

print '</orthoXML>';