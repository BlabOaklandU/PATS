#!/bin/perl.exe
use strict;
use warnings;

use Bio::TreeIO;
use Bio::Tree::Node;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::AnnotatableNode;
use List::MoreUtils qw(firstidx);
use List::Util qw(first);
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;


my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $mcBaseDir =  $ENV{'MolecularClockBaseDir'};
my $mcScriptsDir = $ENV{'MolecularClockScriptDir'};
my $mcsplitGenesDir = $ENV{'splitGenesDir'};


open my $logMolClockDebug, '>>', "$mcBaseDir/molClockDebugLog.txt" or die "Cannot open molClockDebugLog.txt: $!";

print $logMolClockDebug "Running phyToMeg.pl... \n";

chdir $mcsplitGenesDir or die ("Cannot change directory: $mcsplitGenesDir");

#USER INPUTS
# Location of Gene Folders containing .phy files (C:\...\): ;
my $genes_loc = $mcsplitGenesDir;
chomp($genes_loc);

my $outfile_nameA = $ENV{'outFileTreeA'};
chomp($outfile_nameA);
print $logMolClockDebug "outFileTreeA: $outfile_nameA\n";

my $outfile_nameB = $ENV{'outFileTreeB'};
chomp($outfile_nameB);
print $logMolClockDebug "outfileB: $outfile_nameB\n";

my $output_gene_root = $ENV{'rootOutputDirName'};
chomp($output_gene_root);
print $logMolClockDebug "Name of root output gene folder: $output_gene_root";


my $OGDir = "$baseDir/OGGapFilter";
opendir my $OGdirfh, $OGDir  or die "Can't open $OGDir: $!";
my $num_genes =  grep { -f "$OGDir/$_" } readdir($OGdirfh);
print $logMolClockDebug "Number of Genes: $num_genes";

# "Genes to be skipped, seperated by spaces: ";
my $skip_genes = $ENV{'numberOfGenesSkippedMC'};
chomp($skip_genes);
print $logMolClockDebug "Genes to be skipped: $skip_genes";

# SPLITS MULTIPLE GENES TO BE SKIPPED INTO ARRAY 
my @genes_to_skip = split(/ /, $skip_genes);

#open location of gene folders
my $genes_dir = $genes_loc;
opendir (my $genes_dh, $genes_dir) or die $!;
my $skip_num = 0;
for(my $i = 1; $i <= $num_genes; $i++){
	if($skip_num < scalar(@genes_to_skip)){
		if($i == $genes_to_skip[$skip_num]){
			$skip_num += 1;
			next;
		}
	}

my $iterate = 0;
my $tree_gene_loc = $genes_loc."/".$output_gene_root.$i;
my $A_Tree_dir = $tree_gene_loc."/"."A_tree";
my $B_Tree_dir = $tree_gene_loc."/"."B_tree";

print $logMolClockDebug "A TREE DIR: $A_Tree_dir \n";
print $logMolClockDebug "B_TREE DIR: $B_Tree_dir \n\n";

	my $outfile_name_pathA = $A_Tree_dir."/".$outfile_nameA;
	my $outfile_name_pathB = $B_Tree_dir."/".$outfile_nameB;
print $logMolClockDebug "\n***outfile_name_pathA: $outfile_name_pathA***\n";
print $logMolClockDebug "\n***outfile_name_pathB: $outfile_name_pathB***\n";
	
	opendir (my $tree_genes_dh, $tree_gene_loc) or die $!;
	#A-tree and B-tree (iterates twice, through a-tree and b-tree folder)
	while (my $tree_file_name = readdir($tree_genes_dh)){
		if ($tree_file_name eq ".." or $tree_file_name eq "."){
			next;
		}
		#print "\t", $tree_file_name, "\n"; ##REMOVE
		
		my $phy_file_loc = $tree_gene_loc."/".$tree_file_name;
print $logMolClockDebug "\n***tree_gene_loc: $tree_gene_loc***\n";
print $logMolClockDebug "\n***aaaaphy_file_loc: $phy_file_loc***\n";
		opendir(my $phy_dh, $phy_file_loc) or die $!;
		
		#Each .phy file (only iterates through once, for 1 .phy file)
		while (my $phy_file_name = readdir($phy_dh)){
			if ($phy_file_name eq ".." or $phy_file_name eq "."){
				next;
			}
			elsif(index($phy_file_name, ".txt") != -1){
				next;
			}
			elsif(index($phy_file_name, ".csv") != -1){
				next;
			}
			elsif(index($phy_file_name, ".meg") != -1){
				next;
			}
			elsif(index($phy_file_name, ".ctl") != -1){
				next;
			}
			elsif(index($phy_file_name, ".pl") != -1){
				next;
			}
			
			#print "\t\t", $phy_file_name, "\n"; #REMOVE
			
			#open the .phy file and read all lines into an array
			my $phy_path = $phy_file_loc."/".$phy_file_name;
			open(my $fh_phy, "<", $phy_path) or die "Could not open file '$phy_path'";
			chomp(my @lines = <$fh_phy>);
			my $tot_num_arr = scalar @lines;
			#print $tot_num_arr, "\n"; #REMOVE
			close $fh_phy;
			
			#create a .meg file to write to
			my $meg_file_root = substr $phy_path, 0, -4;
			my $meg_file_name = $meg_file_root.".meg";
			print $logMolClockDebug "$meg_file_name\n";
			open(my $fh_meg, ">", $meg_file_name) or die "Could not open file '$meg_file_name'";
			
			print $fh_meg "#Mega\n";
			print $fh_meg "!Title ", $phy_file_name, ";\n\n";
			
			for (my $i = 1; $i<$tot_num_arr; $i++){
				my @spec_seq = split / /, $lines[$i];
				print $fh_meg "#",$spec_seq[0],"\n";
				print $fh_meg $spec_seq[2],"\n\n";
				
			}
			
			
			close $fh_meg;
			
			
			
		};
		closedir($phy_dh);
	}
	closedir($tree_genes_dh);
	
	
}

#CONVERT .PHY FILE INTO .MEG FILE

#RUN .MEG FILE WITH .MAO FILE

#PERFORM CALCULATIONS AND STORE 

#OUTPUT INTO A TEXT FILE TO APPEAR IN SAME FOLDER

closedir($genes_dh);