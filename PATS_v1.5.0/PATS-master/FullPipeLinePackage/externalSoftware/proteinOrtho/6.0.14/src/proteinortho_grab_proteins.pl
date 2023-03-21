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
# 
# @author Paul Klemm
# @email klemmp@staff.uni-marburg.de
# @company Bioinformatics, University of Leipzig
# @version 2
# @date 11-07-2019
#
##########################################################################################

use POSIX;

my $usage = <<'ENDUSAGE';
proteinortho_grab_proteins.pl        greps all genes/proteins of a given fasta file
 
SYNOPSIS
 
proteinortho_grab_proteins.pl (options) QUERY FASTA1 (FASTA2 ...)

	QUERY	proteinortho.tsv FILE or search STRING or '-' for STDIN:
		a)	proteinortho output file (.tsv). This uses by default the -exact option.
		b)	string of one identifier e.g. 'tr|asd3|asd' OR multiple identifier separated by ',' (-F=)
	FASTA*	fasta file(s) (database)

	(options):
		-tofiles, -t  print everything to files instead of stdout files are called OrthoGroup**.fasta for a proteinortho.tsv file
		-E            enables regex matching otherwise the string is escaped (e.g. | -> \|)
		-exact        search patters are extended with a \b, that indicates end of word.
		-source, -s   adds the filename (FASTA1,...) to the found gene-name
		-F=s          char delimiter for multiple identifier if QUERY is a string input (default: ',')

DESCRIPTION
 
	This script finds and extract all given identifier of a list of fasta files. 
	The identifier can be provided as a simple string 'BDNF1', regex string 'BDNF*' 
	or in form of a proteinortho output file (myproject.proteinortho.tsv).
       
	Example:
 
 	# 1. most simple call:

	perl proteinortho_grab_proteins.pl 'BDNF1' *.faa

		STDOUT:
			>BDNF1 Brain derived neurotrophic factor OS=human(...)
			MNNGGPTEMYYQQHMQSAGQPQQPQTVTSGPMSHYPPAQPPLLQPGQPYSHGAPSPYQYG
			>BDNF15 Brain derived neurotrophic factor OS=human(...)
			MAFPLHFSREPAHAIPSMKAPFSRHEVPFGRSPSMAIPNSETHDDVPPPLPPPRHPPCTN

	    The second hit BDNF15 is reported since it also contains 'BDNF1' as a substring. 
	    To prevent such a behaviour use proteinortho_grab_proteins.pl -E 'BDNF1\b'. 
	    The \b marks the end of a word and -E enables regex expressions.

	    Or simply add -exact: perl proteinortho_grab_proteins.pl -exact 'BDNF1' *.faa

 	# 2. multiple ids:

	perl proteinortho_grab_proteins.pl 'BDNF1,BDNF2,BDNF3' *.faa

 	# 3. more complex regex search:

	perl proteinortho_grab_proteins.pl -E 'B?DNF[0-3]3+' *.faa

		This finds: BDNF13, BDNF23, DNF13, DNF033, ... 

 	# 4. proteinortho tsv file and write output to files:

	proteinortho_grab_proteins.pl -tofiles myproject.proteinortho.tsv test/*.faa

		This will produce the files: OrthoGroup0.fasta, OrthoGroup1.fasta, OrthoGroup2.fasta, ...
		Each fasta file contains all genes of one orthology group (one line in myproject.proteinortho.tsv)
 
ENDUSAGE

my $query;
my $help;
my $tofiles=0;
my $justid;
my $prefix=">";
my $doregex=0;
my $source=0;
my $exact=0;
my $del=',';

my @ARGViddone=[];
my $ARGViddone_counter=scalar(@ARGV);
for(my $v = 0 ; $v < scalar @ARGV ; $v++){
	$ARGViddone[$v]=1;
	if($ARGV[$v] =~ m/--?(help|h)$/){$help=1;}
	elsif($ARGV[$v] =~ m/^--?(tofiles|t)$/){$tofiles=1;}
	elsif($ARGV[$v] =~ m/^--?(source|s)$/){$source=1;}
	elsif($ARGV[$v] =~ m/^--?F=(.*)$/){$del=$1;}
	elsif($ARGV[$v] =~ m/^--?E$/){$doregex=1;}
	elsif($ARGV[$v] =~ m/^--?exact$/){$exact=1;}
	elsif($ARGV[$v] =~ m/^-.+/){print $usage; print STDERR "ERROR: invalid option ".$ARGV[$v]."!\n\n";exit(1);}
	elsif(!defined($query)){$query = $ARGV[$v];}
	else{$ARGViddone[$v]=0;$ARGViddone_counter--;}
}
if ($help){
    print $usage;
    exit(0);
}
my $fail="";
if (scalar(@ARGV) == 0){
    $fail.="ERROR: no arguments found!\n";
}
if (!defined($query)){
    $fail.="ERROR: no QUERY provided!\n";
}
if ( $ARGViddone_counter==scalar(@ARGV) ){
    $fail.="ERROR: no FASTA files provided!\n";
}
if($fail){
	print $usage.$fail;
	exit(1);
}

my %qdata;
# my $qdata_count = {};

my $orthogroupcounter=0;
my $genecounter=0;

unless(open(my $FH,'<',$query)) {
	if($query eq "-"){
		foreach my $line (<>) 
		{
			my $linel="";
			#if(!$doregex){$linel=quotemeta($line);}
			if(length($linel)==0|| $linel eq ""){continue;}
			#if($exact){$linel='\b'.$linel.'\b';}
			$qdata{"___"}{$linel}=$linel;
			$genecounter++;
		}
	}else{
		my @sp = split($del,$query); 
		for(my $v = 0 ; $v < scalar @sp ; $v++){
			#if(!$doregex){$sp[$v]=quotemeta($sp[$v]);}
			if(length($sp[$v])==0 || $sp[$v] eq ""){continue;}
			#if($exact){$sp[$v]='\b'.$sp[$v].'\b';}
			$qdata{"___"}{$sp[$v]}=$sp[$v];
			$genecounter++;
		}
	}
}else{
	if(!$exact){print STDERR "[STDERR] WARNING The -exact option is mandatory if a proteinortho file is given. -exact is now set.\n";$exact=1;}
	if($doregex){print STDERR "[STDERR] WARNING The -E option is not allowed if a proteinortho file is given. -E is now unset.\n";$doregex=0;}

	my $query_basename=$query;
	if($query_basename =~ m/\/([^\/]+)$/){
		$query_basename=$1;
	}

	@filenames;
	while(<$FH>){
		$_=~s/[\r\n]+$//;
		my @sp = split(/\t/,$_);
		if(substr($_,0,1) eq "#"){@filenames=@sp; next;}
		if(scalar(@sp)>3){
			for(my $v = 3 ; $v < scalar @sp ; $v++){
				if($sp[$v] eq "*" || $sp[$v] eq ""){next;}
				my @spp = split(",",$sp[$v]);
				#if($exact){$sp[$v]='\b'.$sp[$v].'\b';}
				for(my $vv = 0 ; $vv < scalar @spp ; $vv++){
					#if(!$doregex && !$exact){$spp[$vv]=quotemeta($spp[$vv]);}
					if($spp[$vv] eq "*" || $spp[$vv] eq ""){next;}
					$qdata{$filenames[$v]}{$spp[$vv]}=$query_basename.".OrthoGroup".$orthogroupcounter;
					$genecounter++;
				}
			}
		}
		$orthogroupcounter++;
	}
	close($FH);
	print STDERR "[STDERR] Done reading the query $query file. Now I know $orthogroupcounter groups with $genecounter genes/proteins in total.\n";
}

if( $tofiles==1 && ($orthogroupcounter > 100) ){
	print STDERR "\n!!!\nWARNING : This call will produce $orthogroupcounter files (one for each orthology group) !\nIn the *.html file you can individually extract single groups by clicking on the front part of a row.\n$NC";
	print STDERR "Press 'strg+c' to prevent me from proceeding or wait 20 seconds to continue...\n!!!\n";
  	sleep 20;
	print STDERR "\nWell then, proceeding...\n\n";
	sleep 1;
}

my $cur_gene="";
my $cur_gene_filename="";
my %cur_gene_firsttime;
my $genecounterfound=0;
my $basename = "";
my $numOfFastas=0;

for(my $v = 0 ; $v < scalar @ARGV ; $v++){
	if($ARGViddone[$v]){next;}
	$numOfFastas++;
}

my $fastai=1;

for(my $v = 0 ; $v < scalar @ARGV ; $v++){

	if($ARGViddone[$v]){next;}

	if($tofiles){print STDERR "[STDERR] ($fastai/$numOfFastas) : ";if($basename ne ""){print STDERR "Done reading $basename. "}print STDERR "Start reading the fasta file ".($ARGV[$v])."\n";}
	$fastai++;

	$basename= $ARGV[$v];
	if($basename =~ m/\/([^\/]+)$/){
		$basename=$1;
	}

	open(my $FH,'<',$ARGV[$v]);
	my $geneprintswitch = 0;
	
	while(<$FH>){
		$_=~s/[\r\n]+$//;

		if($_ eq "" || length $_ < 2 || substr($_,0,1) eq "#"){next;}

		my $curLine=$_;

		if(substr($curLine,0,1) eq $prefix ){
			$geneprintswitch = 0;

			if($cur_gene ne ""){
				$cur_gene_filename=~s/\\b//g;
				$cur_gene_filename=~s/[^a-zA-Z0-9.]//g;
				if($tofiles){ # print to files
					my $writemodus=">>";
					if(!exists $cur_gene_firsttime{$cur_gene_filename}){$writemodus=">";$cur_gene_firsttime{$cur_gene_filename}=1;}
					open($FHOUT,$writemodus,$cur_gene_filename);
					print $FHOUT $cur_gene;
					close($FHOUT);
				}else{
					print $cur_gene;
				}
				$cur_gene="";
			}

			if($exact && exists $qdata{$basename}){

				my $genename=$curLine;
				my @arr=split(" ",$genename);
				if(scalar(@arr)>0){$genename=$arr[0];}
				$genename=~s/^>//;

				if(exists $qdata{$basename}{$genename}){
					
					my $headerstr=$curLine;
					if($source){$headerstr=$headerstr." ".$basename;}

					$cur_gene.=$headerstr."\n";
					$cur_gene_filename=$qdata{$basename}{$genename}.".fasta";
					$geneprintswitch = 1;
					$genecounterfound++;
				}

			}else{ # fallback, if the basename (filename) does not exists, try all 

				foreach my $filename (keys %qdata) { 
					if($geneprintswitch){last;}
					foreach my $key (keys %{$qdata{$filename}}) { 
						my $regexv=$key; 
						
						if(!$doregex){$regexv=quotemeta($regexv);}
						if($exact){$regexv='\b'.$regexv.'\b';}

						if( $curLine =~ $regexv ){
									
							my $headerstr=$curLine;
							if($source){$headerstr=$headerstr." ".$basename;}

							$cur_gene.=$headerstr."\n";	
							$cur_gene_filename=$qdata{$filename}{$key}.".fasta";
							$geneprintswitch = 1;
							$genecounterfound++;

							last;
						}
					}
				}	
			}

		}else{
			if($geneprintswitch){
				$cur_gene.=$curLine."\n";
			}
		}
	}
	close($FH);

}

if($cur_gene ne ""){
	$cur_gene_filename=~s/\\b//g;
	$cur_gene_filename=~s/[^a-zA-Z0-9.]//g;
	if($tofiles){ # print to files
		my $writemodus=">>";
		if(!exists $cur_gene_firsttime{$cur_gene_filename}){$writemodus=">";$cur_gene_firsttime{$cur_gene_filename}=1;}
		open($FHOUT,$writemodus,$cur_gene_filename);
		print $FHOUT $cur_gene;
	}else{
		print $cur_gene;
	}
	close($FHOUT);
}

if($genecounter != $genecounterfound){
	print STDERR "[STDERR] WARNING The input ($query) contains $genecounter queries, but I extracted $genecounterfound entries out of the fasta(s).";
	if(!$exact){print STDERR " If this is not desired, please consider using the -exact option";}elsif($genecounter > $genecounterfound){print STDERR "\n-> This should not have happen, maybe some fasta files are missing as input?\n(If you cannot solve this error, please send a report to incoming+paulklemm-phd-proteinortho-7278443-issue-\@incoming.gitlab.com or visit https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error%20Codes for more help. Further more all mails to lechner\@staff.uni-marburg.de are welcome)\n";}
	print STDERR "\n";
}else{
	print STDERR "[STDERR] All entries of the query are found in the fasta(s).\n";
}