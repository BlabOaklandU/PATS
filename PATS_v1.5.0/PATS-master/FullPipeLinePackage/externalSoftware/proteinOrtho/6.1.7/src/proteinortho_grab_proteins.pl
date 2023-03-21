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
# @company Bioinformatics, University of Marburg
# @version 8
# @date 07-06-2022
#
##########################################################################################

use strict;
use threads;
use threads::shared;
use POSIX;
use File::Basename;
use Thread::Queue;
use Cwd 'abs_path';
our $QUEUE = Thread::Queue->new();    # A new empty queue

my $usage = <<'ENDUSAGE';
proteinortho_grab_proteins.pl        greps all genes/proteins of a given fasta file
 
SYNOPSIS
 
proteinortho_grab_proteins.pl (options) QUERY FASTA1 (FASTA2 ...)

	QUERY	proteinortho.tsv FILE or search STRING or '-' for STDIN:
		a)	proteinortho output file (.tsv). This uses by default the -exact option.
		b)	string of one identifier e.g. 'tr|asd3|asd' OR multiple identifier separated by ',' (-F=)
	FASTA*	fasta file(s) (database)

	(options):
		-tofiles, -t  print everything to files instead of STDOUT. Output files are called OrthoGroupX.fasta with X being the Xth group of the proteinortho.tsv file
		-tofiles=DIR  additionally specifies a output directory for the OrthoGroup files
		-E            enables regex matching otherwise the string is escaped (e.g. | -> \|)
		-exact        search patters are extended with a \b, that indicates end of word.
		-cpus=INT     the number of parallel open files for reading, this is strictly limited by the I/O bandwith (default:1).
		              for fast SSD drives, you can increase this to gain speed.
		-minprot=X    if you give a proteinortho.tsv file, this filters out groups with less than X proteins (default:0).
		-source, -s   adds the filename (FASTA1,...) to the found gene-name
		-F=s          char delimiter for multiple identifier if QUERY is a string input (default: ',')
		-isoform      if you use proteinortho with --isoform option, then you need to set this option here too. 
		-core         similar to -tofiles but the output files are named according to the fasta file names instead of the orthology group.
		              Each outputfile contains all proteins of all orthology groups of that species.
		              Implicitly sets -noco: Only one protein is printed for each species and orthology group
		-noco         Omit co-orthologs: Only outputs first entry for each species of each orthology group, i.e. the protein that has the most connections to the group
		                e.g. a group with 6 proteins across 3 species '6  3  0.5  A,B,C  D,E  F' would return sequences of A,D and F
		-singles      Extract only proteins without any orthologs

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



our $query;
our $help;
our $tofiles="-1";
our $isoform=0;
our $justid;
our $prefix=">";
our $doregex=0;
our $source=0;
our $exact=0;
our $del=',';
our @ARGV_copy=@ARGV;
our $ignoreWarning=0;
our $minprot = 0;
our $singles = 0;
our $cpus=1;
our $core = 0;
our $noco = 0;

our @ARGV_copyiddone=(1) x (scalar(@ARGV_copy));
our $ARGV_copyiddone_counter=scalar(@ARGV_copy);
for(my $v = 0 ; $v < scalar @ARGV_copy ; $v++){
	$ARGV_copyiddone[$v]=1;
	if($ARGV_copy[$v] =~ m/--?(help|h)$/){$help=1;}
	elsif($ARGV_copy[$v] =~ m/^--?(tofiles|t)$/){$tofiles="";}
	elsif($ARGV_copy[$v] =~ m/^--?(tofiles|t)=(.*)$/){$tofiles="$2/";}
	elsif($ARGV_copy[$v] =~ m/^--?(source|s)$/){$source=1;}
	elsif($ARGV_copy[$v] =~ m/^--?F=(.*)$/){$del=$1;}
	elsif($ARGV_copy[$v] =~ m/^--?minprots?=(.*)$/){$minprot=$1;}
	elsif($ARGV_copy[$v] =~ m/^--?cpus?=(.*)$/){$cpus=$1;}
	elsif($ARGV_copy[$v] =~ m/^--?E$/){$doregex=1;}
	elsif($ARGV_copy[$v] =~ m/^--?isoform$/){$isoform=1;}
	elsif($ARGV_copy[$v] =~ m/^--?singles?$/){$singles=1;}
	elsif($ARGV_copy[$v] =~ m/^--?ignoreWarning$/){$ignoreWarning=1;}
	elsif($ARGV_copy[$v] =~ m/^--?core$/){$core=1;} # v6.1.0
	elsif($ARGV_copy[$v] =~ m/^--?noco$/){$noco=1;} # v6.1.0
	elsif($ARGV_copy[$v] =~ m/^--?exact$/){$exact=1;}
	elsif($ARGV_copy[$v] =~ m/^-.+/){print $usage; print STDERR "ERROR: invalid option ".$ARGV_copy[$v]."!\n\n";exit(1);}
	elsif(!defined($query)){$query = $ARGV_copy[$v];}
	else{
		# 6.0.32 replace files if there is a *_clean* present
		if(-e $ARGV_copy[$v]){
			my $file_clean = abs_path $ARGV_copy[$v];
			$file_clean=~s/(\.[^.]+)$/_clean$1/;
			if(-e $file_clean){$ARGV_copy[$v] = $file_clean}
		}
		$ARGV_copyiddone[$v]=0;$ARGV_copyiddone_counter--;
	}
}
if ($help){
    print $usage;
    exit(0);
}
my $fail="";
if (scalar(@ARGV_copy) == 0){
    $fail.="ERROR: no arguments found!\n";
}
if (!defined($query)){
    $fail.="ERROR: no QUERY provided!\n";
}
if ($core && $tofiles ne "-1"){
    $fail.="ERROR: -core and -tofile are not compatible!\n";
}
if ( $ARGV_copyiddone_counter==scalar(@ARGV_copy) ){
    $fail.="ERROR: no FASTA files provided!\n";
}
if($fail){
	print $usage.$fail;
	exit(1);
}

if($core){$noco=1}

our %qdata;

our $orthogroupcounter=0;
our $genecounter=0;
my $numOfFastas=0;
my $line_i = 0;

for(my $v = 0 ; $v < scalar @ARGV_copy ; $v++){
	if($ARGV_copyiddone[$v]){next;}
	$numOfFastas++;
}
our @filenames;
our @thread_return :shared = ("")x$cpus;

my $foundHeader=0;

# print STDERR (scalar keys %qdata)." vs ".($numOfFastas);

sub processLine{
	my $line = shift;
	my $prefix_group = shift;

	$line=~s/[\r\n]+$//; 

	my @sp = split(/\t/,$line);
	if(substr($line,0,1) eq "#"){$foundHeader=1;@filenames=@sp; next;}
	if(scalar(@sp)>3){
	 	if( @sp > 3 && $sp[1] < $minprot ){$orthogroupcounter++; next}
	 	if( $core && $sp[0] < $numOfFastas ){next}
	 	if( $singles && $sp[0] == 1 ){next} # skip the singles so we can use !exists later to find them !

	 	for(my $v = 3 ; $v < scalar @sp ; $v++){
			if($sp[$v] eq "*" || $sp[$v] eq ""){next;}
			my @spp = split(",",$sp[$v]);

			for(my $vv = 0 ; $vv < scalar @spp ; $vv++){

				if($spp[$vv] eq "*" || $spp[$vv] eq ""){next;}
				$spp[$vv]=~s/^\(//;$spp[$vv]=~s/\)$//;

				if(!exists $filenames[$v]){ $filenames[$v]="" }

				$qdata{$filenames[$v]}{$spp[$vv]}=$prefix_group.".OrthoGroup".$orthogroupcounter;
				$genecounter++;

				if( $noco ){last}
			}
		}
		$orthogroupcounter++;
	}
}

unless(open(my $FH,'<',$query)) {
	if($query eq "-"){
		foreach my $line (<STDIN>) { &processLine($line, "STDIN") }
	}else{
		my @sp = split($del,$query); 
		for(my $v = 0 ; $v < scalar @sp ; $v++){
			if(length($sp[$v])==0 || $sp[$v] eq ""){next;}
			$sp[$v]=~s/^\(//;$sp[$v]=~s/\)$//;

			my @arr = split(" ",$sp[$v]); # 6.0.32 , -> ; 
	        $arr[0]=~s/,/;/g;
	        $sp[$v]=join(" ",@arr);

			$qdata{"STDIN"}{$sp[$v]}=$sp[$v];
			$genecounter++;
		}
	}
}else{
	if($isoform && $exact){print STDERR "[proteinortho_grab_proteins.pl] WARNING The -isoform option is not compatible with -exact if a proteinortho file is given. -exact is now unset.\n";$exact=0;}
	elsif(!$isoform && !$exact){print STDERR "[proteinortho_grab_proteins.pl] WARNING The -exact option is mandatory if a proteinortho file is given. -exact is now set.\n";$exact=1;}
	if($doregex){print STDERR "[proteinortho_grab_proteins.pl] WARNING The -E option is not allowed if a proteinortho file is given. -E is now unset.\n";$doregex=0;}

	my $query_basename=$query;
	if($query_basename =~ m/\/([^\/]+)$/){
		$query_basename=$1;
	}

	while(<$FH>){ &processLine($_, $query_basename) }
	close($FH);
	print STDERR "[proteinortho_grab_proteins.pl] Done reading the query $query file. Now I know $orthogroupcounter groups with $genecounter genes/proteins in total.\n";
}

if( $foundHeader==0 && $numOfFastas > 3 && $genecounter > 20){
	print STDERR "\nWARNING : The header of the proteinortho file is missing, this can increase the runtime dramatically. Please include the first line (starting with '#'), to accelerate this program.\n\n";
	sleep 1;
}

if( $tofiles ne "-1" && ($orthogroupcounter > 100) && !$ignoreWarning ){
	print STDERR "\n!!!\nWARNING : This call will produce $orthogroupcounter files (one for each orthology group) !\nIn the *.html file you can individually extract single groups by clicking on the front part of a row.\n";
	print STDERR "Press 'strg+c' to prevent me from proceeding or wait 20 seconds to continue...\n!!!\n";
  	sleep 20;
	print STDERR "\nWell then, proceeding...\n\n";
	sleep 1;
}

for (my $v = 0 ; $v < scalar @ARGV_copy ; $v++){ if($ARGV_copyiddone[$v]){next;} $QUEUE->enqueue($ARGV_copy[$v]) }
for (my $i = 0; $i < $cpus; $i++) { $QUEUE->enqueue(undef) }
for (my $i = 0; $i < $cpus; $i++) { threads->create('worker') } # spawn a thread for each core

my %cache_files;
my %master;
my $genecounterfound = 0;
my %genefound;
foreach my $t (threads->list()) {
	$t->join;
	my @ret = split("$;$;",$thread_return[$t->tid]); # await thread
	my $first=shift @ret;

	$genecounterfound += scalar split("$;",$first);
	map { $genefound{$_}=1 } split("$;",$first);

	for (my $i = 0; $i < scalar @ret; $i++) {
		my ($key,$load)=split("$;",$ret[$i]);
		if(!exists $master{$key}){$master{$key}=""}
		$master{$key}.=$load;
	}
}
foreach my $key (keys %master) {
	if($core || $tofiles ne "-1"){
		open(my $FHOUT,">$key");
		print $FHOUT $master{$key};
		close($FHOUT);
	}else{
		print $master{$key};
	}
}

sub worker {
	local $SIG{KILL} = sub { threads->exit };
	my $tid = threads->tid();

	my @genefound;
	my %cache;

	while(my $job = $QUEUE->dequeue()){
		
		my $cur_gene_filename="";
		my %cur_gene_firsttime;
		my $cur_gene="";
		my $basename = basename($job);

		print STDERR "[proteinortho_grab_proteins.pl] Start processing the fasta file ".($job)." (tid=$tid)\n";

		open(my $FH,'<',$job);
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
				
					if($core){
						if(!exists $cache{"$basename.core"}){$cache{"$basename.core"}=""}
						$cache{"$basename.core"} .= $cur_gene;
					}else{
						if(!exists $cache{($tofiles ne "-1" ? $tofiles : "" ).$cur_gene_filename}){$cache{($tofiles ne "-1" ? $tofiles : "" ).$cur_gene_filename}=""}
						$cache{($tofiles ne "-1" ? $tofiles : "" ).$cur_gene_filename} .= $cur_gene;
					}
				
					$cur_gene="";
				}

				if((!$singles && $exact && exists $qdata{$basename}) || $singles){

					my $genename=$curLine;
					my @arr=split(" ",$genename);
					if(scalar(@arr)>0){$genename=$arr[0];}
					$genename=~s/^>//;

					my @arr = split(" ",$genename); # 6.0.32 , -> ; 
			        $arr[0]=~s/,/;/g;
			        my $genename_no_comma=join(" ",@arr);
			        
					if( ($singles && (!exists $qdata{$basename} || !exists $qdata{$basename}{$genename_no_comma})) || 
						(!$singles && exists $qdata{$basename} && exists $qdata{$basename}{$genename_no_comma}) ){
					
						my $headerstr=$curLine;
						if($source){$headerstr=$headerstr." ".$basename;}

						$cur_gene.=$headerstr."\n";
						$cur_gene_filename=$qdata{$basename}{$genename}.".fasta";
						$geneprintswitch = 1;
						push(@genefound,$genename);
					}

				}else{ # fallback, if the basename (filename) does not exists, try all 

					foreach my $filename (keys %qdata) { 
						if($geneprintswitch){last;}
						foreach my $key (keys %{$qdata{$filename}}) { 

							if(!defined $qdata{$filename}{$key}){next}

							my $regexv = $key; 
					        
							my $curLine_test = $curLine;

							if(!$doregex && !$exact){$regexv=quotemeta($regexv);}

							$regexv=~s/;/[,;]/g; # 6.0.32 , -> ; 

							my $test_match = 0;
							if( !$exact ){
								
								# use regular expression
								$test_match = $curLine_test =~ $regexv;

							}else{
								
								# directly compare starting with first 5 character of fasta entry as offset
								my $offset = 1; # start at 1 -> fasta entries starts with ">"
								while( !( $test_match = substr($curLine,$offset,length $key) eq $key ) && $offset < 5 ){
									$offset++;
								} 
							}

							if( $test_match ){

								if( $qdata{$filename}{$key} eq ""){
									print STDERR "[proteinortho_grab_proteins.pl] WARNING The input ($key) was found multiple times in the fasta files ".(!$exact ? "(maybe try --exact)." : ".")."\n";
								}

								my $headerstr=$curLine;
								if($source){$headerstr=$headerstr." ".$basename;}

								$cur_gene.=$headerstr."\n";	
								$cur_gene_filename = $qdata{$filename}{$key}.".fasta";
								$geneprintswitch = 1;
								push(@genefound,$key);

								last;
							}
						}
					}
				}
			}else{ if($geneprintswitch){ $cur_gene.=$curLine."\n" } }
		}
		close($FH);

		if($cur_gene ne ""){
			$cur_gene_filename=~s/\\b//g;
			$cur_gene_filename=~s/[^a-zA-Z0-9.]//g;

			if($core){
				if(!exists $cache{"$basename.core"}){$cache{"$basename.core"}=""}
				$cache{"$basename.core"} .= $cur_gene;
			}else{
				if(!exists $cache{($tofiles ne "-1" ? $tofiles : "" ).$cur_gene_filename}){$cache{($tofiles ne "-1" ? $tofiles : "" ).$cur_gene_filename}=""}
				$cache{($tofiles ne "-1" ? $tofiles : "" ).$cur_gene_filename} .= $cur_gene;
			}
		}
	}

	my $ret=join("$;",@genefound);
	foreach my $key (keys %cache) { $ret.="$;$;$key$;".$cache{$key} }

	$thread_return[$tid]=$ret;
}


if($genecounter != $genecounterfound){
	print STDERR "[proteinortho_grab_proteins.pl] WARNING The input ($query) contains $genecounter queries, but I extracted $genecounterfound entries out of the fasta(s).";
	if(!$exact){print STDERR " If this is not desired, please consider using the -exact option";}
	#elsif($genecounter > $genecounterfound){print STDERR "\n-> This should not have happen, maybe some fasta files are missing as input?\n(If you cannot solve this error, please send a report to incoming+paulklemm-phd-proteinortho-7278443-issue-\@incoming.gitlab.com or visit https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error%20Codes for more help. Further more all mails to lechner\@staff.uni-marburg.de are welcome)\n";}
	print STDERR "\n\n";

	if($genecounterfound < $genecounter){
		print STDERR "The following ids were not found:\n";
		my $counter=0;
		foreach my $filename (keys %qdata) { 
			foreach my $key (keys %{$qdata{$filename}}) { 
				if(10 < $counter++){last}
				print STDERR $key."\n"; 
			}
			if(10 < $counter){
				print STDERR " ...\n";		

				open(FH,">missing_ids.txt");
				foreach my $filenamee (keys %qdata) { 
					foreach my $keyy (keys %{$qdata{$filenamee}}) { 
						if(!exists $genefound{$keyy}){ print FH "$keyy\n" }
					}
				}
				close(FH);
				print STDERR "I produced a file containing all missing ids in the current working directory (missing_ids.txt)\n";

				last
			}
		}
		print STDERR "\nPlease make sure that those ids are part of the given fasta files (try searching these in the given fasta files) !\n";
	}
	if($genecounterfound == 0){
		print STDERR "\nIf you used the --isoform option in proteinortho, then please set -isoform here too (no need to specify the isoform type, e.g. uniprot)!\n";
	}
	
}else{
	print STDERR "[proteinortho_grab_proteins.pl] All entries of the query are found in the fasta(s).\n";
}
