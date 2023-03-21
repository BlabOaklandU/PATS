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
# @version 1
# @date 11-11-2019
#
##########################################################################################

use POSIX;

my $usage = "
proteinortho_history.pl        reports the history of a (or a pair of) gene/protein(s).
 
SYNOPSIS
 
proteinortho_history.pl (-project=myproject) QUERY (FASTA1 FASTA2 ...)

	QUERY	A string of a single gene/protein or 2 separated by a comma or a whitespace (the input is escaped using quotemeta, use -noquotemeta to avoid this)

	-project=MYPROJECT	The project name (as specified in proteinortho with -project) (default:auto detect in the current directory)
	-step=[123] 		(optional) If specified more optput is printed (to STDOUT) for the given step:
		-step=1 : search for the given fasta sequence in the input fasta files
		-step=2 : search in the *.blast-graph
		-step=3 : search in the *.proteinortho file 
		-step=all : prints everything of above to STDOUT
	FASTA*						(optional) input fasta files 
	-noquotemeta, -E			(optional) If set, then the query will not be escaped.
	-plain, -p, -notableformat	(optional) If -step= is set too, then the tables are not formatted and a plain csv is printed instead. 
	-delim= 					(optional) Defines the delimiter character for spliting the query (if you want to search for 2 genes/proteins)

	NOTE: if you use the -keep option and you have the project_cache_proteinortho/ directory, this program additionally searches for all blast hits.

";

our $maxNumOfCharsInOneLine=`tput cols`;
chomp($maxNumOfCharsInOneLine);
if($maxNumOfCharsInOneLine<10){$maxNumOfCharsInOneLine=160;}
our $split_delim="[:\t]";
our @spl_header;
our @spl;
our $last_isHeaderLine=0;
our $last_isHeaderLine=0;$isHeaderLine=1;
our $noheader=0;

my $query;
my $query_esc;
my $help;
my $project="myproject";
my $step="none";
my $delim="[, ]";
my $do_quotemeta=1;
our $notableformat=0;

my @ARGViddone=[];
my $ARGViddone_counter=scalar(@ARGV);
for(my $v = 0 ; $v < scalar @ARGV ; $v++){
	
	if($ARGV[$v] =~ m/^--?(help|h)$/){$help=1;$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?project=?(.*)$/){$project=$1;chomp($project);$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?step=?([123]|(all))$/){$step=$1;chomp($step);$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?delim=?(.*)$/){$delim=$1;chomp($delim);$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?noquotemeta$/){$do_quotemeta=0;$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?E$/){$do_quotemeta=0;$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?notableformat$/){$notableformat=1;$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?plain$/){$notableformat=1;$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^--?p$/){$notableformat=1;$ARGViddone[$v]=1;}
	elsif($ARGV[$v] =~ m/^-.+/){print $usage; print STDERR "ERROR: invalid option ".$ARGV[$v]."!\n\n";$ARGViddone[$v]=1;exit(1);}
	elsif(!defined($query)){$query = $ARGV[$v];}
	else{$ARGViddone[$v]=0;$ARGViddone_counter++;}

}

my $project_corrected="";

if ($help){
    print $usage;
    exit(0);
}
my $fail="";
if (!defined($query)){
    $fail.="ERROR: no QUERY provided!\n";
}

if (! -e $project.".blast-graph" && ! -e $project.".proteinortho.tsv" && ! -e $project.".proteinortho-graph" ){ # && ! -d "proteinortho_cache_".$project 
	@blastglob=glob("*.blast-graph");
	@proteinorthoglob=glob("*.proteinortho.tsv");
	if (scalar(@blastglob)>0){
		$project_corrected=$project;
		$project=$blastglob[0];$project=~s/.blast-graph$//g;
	}elsif(scalar(@proteinorthoglob)>0){
		$project_corrected=$project;
		$project=$proteinorthoglob[0];$project=~s/.blast-graph$//g;
	}else{
		 $fail.="ERROR: '".$project."' files not found, make sure you are in the directory with the proteinortho output! Specify the project name (prefix) with -project=... \n";
	}
}


if($do_quotemeta){
	$query_esc=quotemeta($query);
	$query_esc=~s/\\([, ])/$1/;
	$query_esc=~s/\\($delim)/$1/;
	if($query ne $query_esc){print STDERR "(escaped the input query '$query' -> '$query_esc')\n";}
}else{
	$query_esc=$query;
}

my @query_split=split($delim,$query_esc);
my @query_split_noesc=split($delim,$query);

if (scalar(@query_split) > 2){
    $fail.="ERROR: the QUERY contains too many delimiters ([, ]), specify the delimiter character with -delim=... !\n";
}

if($fail){
	print $usage.$fail;
	exit(1);
}

print STDERR "\nA short history of '"; if(scalar(@query_split)==2){print STDERR $query_split_noesc[0]."' and '".$query_split_noesc[1];}else{print STDERR $query_split_noesc[0]} print STDERR "' in '$project':\n";

my $evalue=1e-5;
my $cov=50;
my $conn=0.1;
my $identity=25;
my $sim=0.95;
if(-e $project.".info"){
	open(my $FHinfo,"<$project.info") or die "ERROR: file not found $! $?\n";
	while(<$FHinfo>){
		if(/Parameter-vector.*/){
			if(/evalue=([^,]+)/){$evalue=$1}
			if(/coverage=([^,]+)/){$cov=$1}
			if(/connectivity=([^,]+)/){$conn=$1}
			if(/identity=([^,]+)/){$identity=$1}
			if(/sim=([^,]+)/){$sim=$1}
		}
	}
	close($FHinfo);
}

print STDERR "\n:::::::[ input files ]::::::::::::::::\n\n";
#
# input files
#

my $did_found_a_in_input=0;
my %gene_len;
my $did_found_b_in_input=0;

for(my $v = 0 ; $v < scalar @ARGV ; $v++){
	if($ARGViddone[$v] || $ARGV[$v] eq $query){next;}	
	open(my $FHfasta,"<$ARGV[$v]") or die "ERROR: file not found $! $?\n";
	my $cur_gene="";
	my $cur_len=0;
	while(<$FHfasta>){
		chomp;
		if(/^[^>]/){
			$cur_len+=length($_);
		}elsif(/>/){
			if($cur_gene ne ""){$gene_len{$cur_gene}=$cur_len}
			$cur_gene=$_;
			$cur_gene=~s/>| .*//g;
			$cur_len=0;
		}
	}
	close($FHfasta);
	if($cur_gene ne ""){$gene_len{$cur_gene}=$cur_len}
}

# foreach my $key (keys %gene_len){print STDERR "$key -> ".$gene_len{$key}."\n";}die;

if($step eq "1" || $step eq "none" || $step eq "all"){

	if ( $ARGViddone_counter != scalar(@ARGV) ){
		print STDERR "checking the input fasta files...\n";

		for(my $v = 0 ; $v < scalar @ARGV ; $v++){

			if($ARGViddone[$v] || $ARGV[$v] eq $query){next;}	
			if(! -e $ARGV[$v]){print STDERR " ! WARNING ! did not find file '".$ARGV[$v]."', proceeding anyway...\n"}
			elsif(! -r $ARGV[$v]){print STDERR "! WARNING ! cannot read '".$ARGV[$v]."', proceeding anyway...\n"}

			$result=`grep -nH '\\b$query_split[0]\\b' $ARGV[$v] 2>/dev/null`;
			$num_result=scalar split("\n",$result);

			if($num_result>0 && $did_found_a_in_input!=0){
				print STDERR "! WARNING ! found '".$query_split_noesc[0]."' again in ".$ARGV[$v].". (proteinortho can handle this) \n";
			}
			if($num_result==0){
				# did not found the query in this file ...
			}elsif($num_result==1 && $did_found_a_in_input==0){
				$did_found_a_in_input=1;
				print STDERR "found '".$query_split_noesc[0]."' (length=".$gene_len{$query_split_noesc[0]}.") in ".$ARGV[$v].". \n";
				if($step ne "none"){print `proteinortho_grab_proteins.pl -E '\\b$query_split[0]\\b' $ARGV[$v]`."\n";}
			}elsif($num_result>1){
				print STDERR "! ERROR ! found '".$query_split_noesc[0]."' multiple times in ".$ARGV[$v].". This seems dangerous, make sure there are no duplicates present in your fasta files ! \n";
				if($step ne "none"){print `proteinortho_grab_proteins.pl -E '\\b$query_split[0]\\b' $ARGV[$v]`."\n";}
			}

			if(scalar(@query_split)==2){

				$result=`grep -nH '\\b$query_split[1]\\b' $ARGV[$v] 2>/dev/null`;
				$num_result=scalar split("\n",$result);

				if($num_result>0 && $did_found_b_in_input!=0){
					print STDERR "! WARNING ! found '".$query_split_noesc[1]."' again in ".$ARGV[$v].". (proteinortho can handle this) \n";
					if($step ne "none"){print "$result\n";}
				}
				if($num_result==0){
					# did not found the query in this file ...
				}elsif($num_result==1 && $did_found_b_in_input==0){
					print STDERR "found '".$query_split_noesc[1]."' in ".$ARGV[$v].". \n";
					$did_found_b_in_input=1;
					if($step ne "none"){print `proteinortho_grab_proteins.pl -E '\\b$query_split[1]\\b' $ARGV[$v]`."\n";}

				}elsif($num_result>1){
					print STDERR "! ERROR ! found '".$query_split_noesc[1]."' multiple times in ".$ARGV[$v].". This seems dangerous, make sure there are no duplicates present in your fasta files ! \n";
					if($step ne "none"){print `proteinortho_grab_proteins.pl -E '\\b$query_split[1]\\b' $ARGV[$v]`."\n";}
				}
			}
		}

		if(!$did_found_a_in_input){
			print STDERR "! ERROR ! did not found '".$query_split_noesc[0]."' in the provided fasta files...\n";
		}
		if(scalar(@query_split)==2 && !$did_found_b_in_input){
			print STDERR "! ERROR ! did not found '".$query_split_noesc[1]."' in the provided fasta files...\n";
		}
	}else{
		print STDERR "input fasta files are not provided: skipping\n";
	}
if($step eq "none"){print STDERR "\n(use -step=1 for more details)\n";}
}

print STDERR "\n:::::::[ all versus all blast ]:::::::\n\n";
#
# step 2 blast-graph
# TODO the single blast graph ... 

my $proteinortho_cache=0;
my $num_a_all_blast_hits =0;
my $num_b_all_blast_hits =0;
my $num_ab_all_blast_hits =0;
my $num_a_aRBH_hits =0;
my $num_b_aRBH_hits =0;
my $num_ab_aRBH_hits =0;

if($step eq "2" || $step eq "none" || $step eq "all"){

	if(-d "proteinortho_cache_".$project){
		$proteinortho_cache=1;
		print STDERR "checking all blast files in 'proteinortho_cache_".$project."'...\n";

		$result=`grep -nH '\\b$query_split[0]\\b' proteinortho_cache_$project/*.vs.* 2>/dev/null`;
		$num_a_all_blast_hits=scalar split("\n",$result);

		print STDERR "the query '".$query_split_noesc[0]."' got $num_a_all_blast_hits hit(s).\n";

		{
			my $num_a_all_blast_hits_filter=scalar grep {my @a=split("\t",$_); scalar(@a)>10 && $a[10]<$evalue} split("\n",$result);
			print STDERR "$num_a_all_blast_hits_filter/$num_a_all_blast_hits are significant (-e=$evalue).\n";
		}
		{
			#print STDERR join "\n", grep {my @a=split("\t",$_); scalar(@a)>10 && exists $gene_len{$a[1]}} map {$_=~s/^[^:]+:[0-9]+://g;$_} split("\n",$result);
			#die;
			my $num_a_all_blast_hits_filter=scalar grep {my @a=split("\t",$_); scalar(@a)>10 && exists $gene_len{$a[0]} && exists $gene_len{$a[1]} && $a[3]/$gene_len{$a[0]}>$cov/100 && $a[3]/$gene_len{$a[1]}>$cov/100} map {$_=~s/^[^:]+:[0-9]+://g;$_} split("\n",$result);
			print STDERR "$num_a_all_blast_hits_filter/$num_a_all_blast_hits are covered enough (-cov=$cov).\n";
		}
		{
			my $num_a_all_blast_hits_filter=scalar grep {my @a=split("\t",$_); scalar(@a)>10 && $a[10]<$evalue && exists $gene_len{$a[0]} && exists $gene_len{$a[1]} && $a[3]/$gene_len{$a[0]}>$cov/100 && $a[3]/$gene_len{$a[1]}>$cov/100} map {$_=~s/^[^:]+:[0-9]+://g;$_} split("\n",$result);
			print STDERR "$num_a_all_blast_hits_filter/$num_a_all_blast_hits are significant and covered enough (combined).\n";
		}

		if($step ne "none"){$result="# file\tline\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n".$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

		if(scalar @query_split == 2){

			$result=`grep -nH '\\b$query_split[1]\\b' proteinortho_cache_$project/*.vs.* 2>/dev/null`;
			$num_b_all_blast_hits=scalar split("\n",$result);

			print STDERR "the query '".$query_split_noesc[1]."' got $num_b_all_blast_hits hit(s).\n";
			if($step ne "none"){$result="# file\tline\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n".$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

			$result=`grep -nHE '(\\b$query_split[0]\\b\\t\\b$query_split[1]\\b)|(\\b$query_split[1]\\b\\t\\b$query_split[0]\\b)' proteinortho_cache_$project/*.vs.* 2>/dev/null`;
			$num_ab_all_blast_hits=scalar split("\n",$result);

			print STDERR "both queries '".$query_split_noesc[0]."' and '".$query_split_noesc[1]."' have $num_ab_all_blast_hits hit(s) with each other use -step=2 for more details).\n";
			if($step ne "none"){$result="# file\tline\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n".$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

			#print "grep -nH '(\\b$query_split[0]\\b\\t\\b$query_split[1]\\b)|(\\b$query_split[1]\\b\\t\\b$query_split[0]\\b)' proteinortho_cache_$project/*.vs.* 2>/dev/null";
		}

	}elsif(glob("*.vs.*")){
		$proteinortho_cache=1;
		print STDERR "checking all blast files...\n";

		$result=`grep -nH '\\b$query_split[0]\\b' *.vs.* 2>/dev/null`;
		$num_a_all_blast_hits=scalar split("\n",$result);

		print STDERR "the query '".$query_split_noesc[0]."' got $num_a_all_blast_hits hit(s).\n";
		if($step ne "none"){$result="# file\tline\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n".$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

		if(scalar @query_split == 2){

			$result=`grep -nH '\\b$query_split[1]\\b' *.vs.* 2>/dev/null`;
			$num_b_all_blast_hits=scalar split("\n",$result);

			print STDERR "the query '".$query_split_noesc[1]."' got $num_b_all_blast_hits hit(s).\n";
			if($step ne "none"){$result="# file\tline\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n".$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

			$result=`grep -nHE '(\\b$query_split[0]\\b\\t\\b$query_split[1]\\b)|(\\b$query_split[1]\\b\\t\\b$query_split[0]\\b)' *.vs.* 2>/dev/null`;
			$num_ab_all_blast_hits=scalar split("\n",$result);

			print STDERR "both queries '".$query_split_noesc[0]."' and '".$query_split_noesc[1]."' have $num_ab_all_blast_hits hit(s) with each other use -step=2 for more details).\n";
			if($step ne "none"){$result="# file\tline\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n".$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

			#print "grep -nH '(\\b$query_split[0]\\b\\t\\b$query_split[1]\\b)|(\\b$query_split[1]\\b\\t\\b$query_split[0]\\b)' *.vs.* 2>/dev/null";
		}

	}else{
		print STDERR "did not found the temporary files *.vs.* neither in 'proteinortho_cache_".$project."/' nor in the current directory. If you want to analyse all blast files please use '-keep' option with proteinortho or provide the *.vs.* files in the current directory...\n";
	}
print STDERR "\n:::::::[ reciprocal best hit ]::::::::\n\n";

	if(-e $project.".blast-graph"){
		print STDERR "checking the blast-graph '$project.blast-graph' (reciprocal adaptive best hit graph)...\n";

		$result=`grep -nH '\\b$query_split[0]\\b' $project.blast-graph 2>/dev/null`;
		$num_a_aRBH_hits=scalar split("\n",$result);

		print STDERR "the query '".$query_split_noesc[0]."' got $num_a_aRBH_hits reciprocal hit(s) with -sim=$sim.\n";
		if($step ne "none"){$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",`head -n2 $project.blast-graph |tail -n1| perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result)){processLine($_);}print "\n";}

		if(scalar @query_split == 2){

			$result=`grep -nH '\\b$query_split[1]\\b' $project.blast-graph 2>/dev/null`;
			$num_b_aRBH_hits=scalar split("\n",$result);

			print STDERR "the query '".$query_split_noesc[1]."' got $num_b_aRBH_hits reciprocal hit(s).\n";
			if($step ne "none"){$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",`head -n2 $project.blast-graph |tail -n1| perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result)){processLine($_);}print "\n";}

			$result=`grep -nHE '(\\b$query_split[0]\\b\\t\\b$query_split[1]\\b)|(\\b$query_split[1]\\b\\t\\b$query_split[0]\\b)' $project.blast-graph 2>/dev/null`;
			$num_ab_aRBH_hits=scalar split("\n",$result);

			print STDERR "both queries '".$query_split_noesc[0]."' and '".$query_split_noesc[1]."' have $num_ab_aRBH_hits reciprocal hit(s) with each other.\n";
			if($step ne "none"){$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",`head -n2 $project.blast-graph |tail -n1| perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result)){processLine($_);}print "\n";}
		}

	}else{
		print STDERR "did not found the blast-graph '$project.blast-graph': skipping...\n";
	}
if($step eq "none"){print STDERR "\n(use -step=2 for more details)\n";}
}

print STDERR "\n:::::::[ clustering ]:::::::::::::::::\n\n";
#
# step 3 proteinortho
#

my $num_a_cluster_hits =0;
my $num_b_cluster_hits =0;
my $num_ab_cluster_hits =0;
my $num_a_cluster_groups =0;
my $num_b_cluster_groups =0;
my $num_ab_cluster_groups =0;

if($step eq "3" || $step eq "none" || $step eq "all"){

	if(-e $project.".proteinortho-graph"){
		print STDERR "checking the proteinortho-graph '$project.proteinortho-graph' (result of the clustering)...\n";

		$result=`grep -nH '\\b$query_split[0]\\b' $project.proteinortho-graph 2>/dev/null`;
		$num_a_cluster_hits=scalar split("\n",$result);

		print STDERR "the query '".$query_split_noesc[0]."' got $num_a_cluster_hits putative ortholog(s) with -conn=$conn.\n";
		if($step ne "none"){$result=`echo '# file\tline\t'|tr -d '\n'; head -n2 $project.proteinortho-graph| tail -n1 | sed 's/#//g'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

		if(scalar @query_split == 2){

			$result=`grep -nH '\\b$query_split[1]\\b' $project.proteinortho-graph 2>/dev/null`;
			$num_b_cluster_hits=scalar split("\n",$result);

			print STDERR "the query '".$query_split_noesc[1]."' got $num_b_cluster_hits putative ortholog(s).\n";
			if($step ne "none"){$result=`echo '# file\tline\t'|tr -d '\n'; head -n2 $project.proteinortho-graph| tail -n1 | sed 's/#//g'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

			$result=`grep -nHE '(\\b$query_split[0]\\b\\t\\b$query_split[1]\\b)|(\\b$query_split[1]\\b\\t\\b$query_split[0]\\b)' $project.proteinortho-graph 2>/dev/null`;
			$num_ab_cluster_hits=scalar split("\n",$result);

			if($num_ab_cluster_hits!=0){
				print STDERR "both queries '".$query_split_noesc[0]."' and '".$query_split_noesc[1]."' are putative ortholog(s).\n";
				if($step ne "none"){$result=`echo '# file\tline\t'|tr -d '\n'; head -n2 $project.proteinortho-graph| tail -n1 | sed 's/#//g'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}
			}else{
				print STDERR "both queries '".$query_split_noesc[0]."' and '".$query_split_noesc[1]."' are NOT putative ortholog(s).\n";
			}
			
		}

	}else{
		print STDERR "did not found the proteinortho-graph '$project.proteinortho-graph': skipping...\n";
	}
print STDERR "\n:::::::[ clustering groups ]::::::::::\n\n";


	if(-e $project.".proteinortho.tsv"){
		print STDERR "checking the proteinortho.tsv file '$project.proteinortho.tsv' (result of the clustering)...\n";

		$result=`grep -nH '\\b$query_split[0]\\b' $project.proteinortho.tsv 2>/dev/null`;
		$num_a_cluster_groups=scalar split("\n",$result);

		print STDERR "the query '".$query_split_noesc[0]."' is part of $num_a_cluster_groups group(s) of putative orthologs.\n";
		if($step ne "none"){$result=`head -n1 $project.proteinortho.tsv | perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

		if(scalar @query_split == 2){

			$result=`grep -nH '\\b$query_split[1]\\b' $project.proteinortho.tsv 2>/dev/null`;
			$num_b_cluster_groups=scalar split("\n",$result);

			print STDERR "the query '".$query_split_noesc[1]."' is part of $num_b_cluster_groups group(s) of putative orthologs.\n";
			if($step ne "none"){$result=`head -n1 $project.proteinortho.tsv | perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

			$result=`grep -nHE '(\\b$query_split[0]\\b.*\\b$query_split[1]\\b)|(\\b$query_split[1]\\b.*\\b$query_split[0]\\b)' $project.proteinortho.tsv 2>/dev/null`;
			$num_ab_cluster_groups=scalar split("\n",$result);

			print STDERR "both queries '".$query_split_noesc[0]."' and '".$query_split_noesc[1]."' have $num_ab_cluster_groups group(s) in common.\n";
			if($step ne "none"){$result=`head -n1 $project.proteinortho.tsv | perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

		}

	}elsif(-e $project.".proteinortho"){
		print STDERR "checking the proteinortho file '$project.proteinortho' (result of the clustering)...\n";

		$result=`grep -nH '\\b$query_split[0]\\b' $project.proteinortho 2>/dev/null`;
		$num_a_cluster_groups=scalar split("\n",$result);

		print STDERR "the query '".$query_split_noesc[0]."' is part of $num_a_cluster_groups group(s) of putative orthologs.\n";
		if($step ne "none"){$result=`head -n1 $project.proteinortho | perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

		if(scalar @query_split == 2){

			$result=`grep -nH '\\b$query_split[1]\\b' $project.proteinortho 2>/dev/null`;
			$num_b_cluster_groups=scalar split("\n",$result);

			print STDERR "the query '".$query_split_noesc[1]."' is part of $num_b_cluster_groups group(s) of putative orthologs.\n";
			if($step ne "none"){$result=`head -n1 $project.proteinortho | perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

			$result=`grep -nHE '(\\b$query_split[0]\\b.*\\b$query_split[1]\\b)|(\\b$query_split[1]\\b.*\\b$query_split[0]\\b)' $project.proteinortho 2>/dev/null`;
			$num_ab_cluster_groups=scalar split("\n",$result);

			print STDERR "both queries '".$query_split_noesc[0]."' and '".$query_split_noesc[1]."' have $num_ab_cluster_groups group(s) in common.\n";
			if($step ne "none"){$result=`head -n1 $project.proteinortho | perl -lne 'chomp;s/#//g;print "# file\tline\t\$_"'`.$result;$noheader=0;$last_isHeaderLine=0;$isHeaderLine=1;@spl_header=();@spl=();foreach(split("\n",$result)){processLine($_);}print "\n";}

		}

	}else{
		print STDERR "did not found the proteinortho.tsv '$project.proteinortho.tsv': skipping...\n";
	}

if($step eq "none"){print STDERR "\n(use -step=3 for more details)\n";}

}
print STDERR "\n:::::::[ summary ]::::::::::::::::::::\n\n";
if(scalar @query_split == 2){
	if($did_found_b_in_input && $did_found_a_in_input){print STDERR "Found both queries in the input files.\n"}else{print STDERR "Did found the two queries in the input files...\n\n";exit(0);}
	if($proteinortho_cache){
		if($num_ab_all_blast_hits){print STDERR "Both queries hit each other"}else{print STDERR "Sadly, the two queries did not hit each other using the blast algorithm...\n\n";exit(0);}
		if($num_ab_aRBH_hits){print STDERR " AND they are adaptive reciprocal best hits.\n"}else{print STDERR " but they are NOT adaptive reciprocal best hits...\n\n";exit(0);}
	}else{
		if($num_ab_aRBH_hits){print STDERR "Furthermore, Both queries are adaptive reciprocal best hits.\n"}else{print STDERR "Sadly, both queries are NOT adaptive reciprocal best hits...\n\n";exit(0);}	
	}
	if($num_ab_cluster_groups){print STDERR "Finally, both queries are putative orthologs, since they occure in the same group after the clustering step.\n"}else{print STDERR "Sadly, both queries are NOT putative orthologs, since they are not occuring in the same group  after the clustering step....\n\n";exit(0);}	
}else{
	if($did_found_a_in_input){print STDERR "Found the query in the input files.\n"}else{print STDERR "I did NOT find the query in the input files...\n\n";exit(0);}
	if($proteinortho_cache){
		if($num_a_all_blast_hits){print STDERR "The query hit $num_a_all_blast_hits other protein(s)/gene(s) (blast-graph)"}else{print STDERR "the query did not hit anything (blast-graph)...\n\n";exit(0);}
		if($num_a_aRBH_hits){print STDERR " AND $num_a_aRBH_hits of these hits are adaptive reciprocal best hits.\n"}else{print STDERR " BUT none of these hits are adaptive reciprocal best hits...\n\n";exit(0);}
	}else{
		if($num_a_aRBH_hits){print STDERR "The query is part of the adaptive reciprocal best hit graph.\n"}else{print STDERR "The query is NOT part of the adaptive reciprocal best hit graph...\n\n";exit(0);}	
	}
	if($num_a_cluster_groups){print STDERR "Furthermore, the query is a putative ortholog with $num_a_cluster_hits other protein(s)/gene(s)\n"}else{print STDERR "Sadly, the query is NOT a putative ortholog of any other protein(s)/gene(s)....\n\n";exit(0);}	
}

print STDERR "\n";
if($project_corrected ne ""){
	print STDERR "WARNING: The project '".$project_corrected."' was not found, I automatically detected '$project' ! (Specify the project name (prefix) with -project=...)\n";
	print STDERR "\n";
}



sub processLine{
	$_=shift;
	if($notableformat == 1){print "$_\n";next;}
	chomp;
	if(length($_)<1){next;}

	@spl;

	if($split_delim eq ""){
		@spl_t=split("\t",$_);
		@spl_c=split(",",$_);
		@spl_s=split(";",$_);

		if(scalar @spl_t < 2 && scalar @spl_c < 2 && scalar @spl_s < 2){next;}

		if(scalar @spl_t > scalar @spl_c && scalar @spl_t > scalar @spl_s ){ @spl = @spl_t; $split_delim='\t';}
		elsif(scalar @spl_c > scalar @spl_t && scalar @spl_c > scalar @spl_s ){ @spl = @spl_c; $split_delim=",";}
		elsif(scalar @spl_s > scalar @spl_c && scalar @spl_s > scalar @spl_t ){ @spl = @spl_s; $split_delim=";";}

	}else{
		@spl=split($split_delim,$_);
	}

	@spl_backup=@spl;

	if(scalar @spl_header > 0 && scalar @spl != scalar @spl_header){$isHeaderLine=1;}
	if(scalar @spl < 2 ){next;}
	if(substr($spl[0],0,1) eq "#"){$spl[0]=~s/^# ?//g;}
	if(scalar(@spl)*2-1>$maxNumOfCharsInOneLine){$maxNumOfCharsInOneLine= -1+2*scalar @spl;print STDERR "Corrected minimum table width: -w=$maxNumOfCharsInOneLine such that at least 1 character per column is displayed.\n";}

	$sumOfCharsLine=length(join("",@spl));

	if($isHeaderLine){ # is a header row 
		while(($sumOfCharsLine + scalar @spl-1) > $maxNumOfCharsInOneLine){ # shave of chars from widest cell
			$max_l=0;
			@max_l_is;
			for (my $i = 0; $i < scalar @spl; $i++) {
				if($max_l < length $spl[$i]){$max_l=length $spl[$i];@max_l_is=();push(@max_l_is,$i)}elsif($max_l == length $spl[$i]){push(@max_l_is,$i)}
			}
			for (my $i = 0; $i < scalar @max_l_is; $i++) {
				if(length $spl[$max_l_is[$i]] > 8 && substr($spl[$max_l_is[$i]],-3) ne "..." ){
					$spl[$max_l_is[$i]]=substr($spl[$max_l_is[$i]],0,length($spl[$max_l_is[$i]])-3-1)."..."
				}
				else{
					$spl[$max_l_is[$i]]=substr($spl_backup[$max_l_is[$i]],0,length($spl[$max_l_is[$i]])-1)
				}
			}
			$sumOfCharsLine=length(join("",@spl));
		}


		while(($sumOfCharsLine + scalar @spl-1) < $maxNumOfCharsInOneLine ){ # add of chars to smallest cell
			$min_l=$maxNumOfCharsInOneLine*10;
			@min_l_is;
			for (my $i = 0; $i < scalar @spl; $i++) {
				if($min_l > length $spl[$i]){$min_l=length $spl[$i];@min_l_is=();push(@min_l_is,$i)}
			}
			for (my $i = 0; $i < scalar @min_l_is; $i++) {

				$leftPad=0;
				$rightPad=0;
				if($spl[$min_l_is[$i]]=~m/( +)$/){$rightPad=length $1}
				if($spl[$min_l_is[$i]]=~m/^( +)/){$leftPad=length $1}

				if( $leftPad < $rightPad ){
					$spl[$min_l_is[$i]]=" ".$spl[$min_l_is[$i]];
				}else{
					$spl[$min_l_is[$i]]=$spl[$min_l_is[$i]]." ";
				}
				
			}
			$sumOfCharsLine=length(join("",@spl));
		}

		@spl_header=@spl;

	}else{ # is not headerline -> do the same as in headerline
		
		while(scalar @spl > scalar @spl_header){pop @spl;}

		for (my $i = 0; $i < scalar @spl; $i++) {
			while(length $spl[$i]< length $spl_header[$i]){ # add pads
				$leftPad=0;
				$rightPad=0;
				if($spl[$i]=~m/( +)$/){$rightPad=length $1}
				if($spl[$i]=~m/^( +)/){$leftPad=length $1}

				if( $leftPad < $rightPad ){
					$spl[$i]=" ".$spl[$i];
				}else{
					$spl[$i]=$spl[$i]." ";
				}
			}
			while(length $spl[$i]>length $spl_header[$i]){ # trim
				if(length $spl[$i] > 5 && substr($spl[$i],-3) ne "..." ){
					$spl[$i]=substr($spl[$i],0,length($spl[$i])-3-1)."..."
				}
				else{
					$spl[$i]=substr($spl_backup[$i],0,length($spl[$i])-2)."#"
				}
			}
		}
	}

	if($isHeaderLine && !$last_isHeaderLine ){$tmp=join("|",@spl);$tmp=~s/\|/+/g;$tmp=~s/[^+]/-/g; print "$tmp\n";}
	print join("|",@spl);
	if($isHeaderLine ){print "\n";$tmp=join("|",@spl);$tmp=~s/\|/+/g;$tmp=~s/[^+]/-/g; print "$tmp";}
	print "\n";
	$last_isHeaderLine=$isHeaderLine;
	$isHeaderLine=0;


}

