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
my $ogindvGeneSequences = $ENV{'indvGenSequencesDir'};


open my $logMolClockDebug, '>>', "$mcBaseDir/molClockDebugLog.txt" or die "Cannot open molClockDebugLog.txt: $!";


print $logMolClockDebug "\nRunning Gene_Splitter.pl...\n\n";

#################################
#USER INPUTS                    #
#################################
#
#
#

#Location of A Tree
my $split_dirA = $ENV{'splitTreeADir'};
chomp($split_dirA);

#Location of B Tree
my $split_dirB = $ENV{'splitTreeBDir'};
chomp($split_dirB);

print $logMolClockDebug "\n";
print $logMolClockDebug "$ENV{'splitTreeADir'}";
print $logMolClockDebug "\n";

#Name of A-tree file
my @a_tree_in = glob ( "$ENV{'splitTreeADir'}/*.nwk" );

#Name of B-tree file
my @b_tree_in = glob ( "$ENV{'splitTreeBDir'}/*.nwk" );

#Location of .fasta files
my $og_locations = $ENV{'orthoGenSequencesDir'};
chomp($og_locations);

#Name of root fasta file
my $root_fasta_name = $ENV{'rootFastaFile'};
chomp($root_fasta_name);

#Location of output gene folders (C:\\...\\)
my $output_location = $mcsplitGenesDir."/";
chomp($output_location);

#Name of root output gene folder
my $output_gene_root = $ENV{'rootOutputDirName'};
chomp($output_gene_root);

#Number of Genes
my $OGDir = "$baseDir/OGGapFilter";
opendir my $OGdirfh, $OGDir  or die "Can't open $OGDir: $!";
my $num_genes =  grep { -f "$OGDir/$_" } readdir($OGdirfh);


#Genes to be skipped, seperated by spaces
my $skip_genes = $ENV{'numberOfGenesSkippedMC'};
chomp($skip_genes);

# SPLITS MULTIPLE GENES TO BE SKIPPED INTO ARRAY 
my @genes_to_skip = split(/ /, $skip_genes);

#
#
#DONE OBTAINING USER INPUTS
##################################

#Copy over CleanedGap Files
	my @ogFiles = glob("$OGDir/*.fasta");
	for my $ogFile (@ogFiles) {
		copy($ogFile, $ogindvGeneSequences) or die "Copy failed: $!";
	}

chdir($mcsplitGenesDir) or die;

#concatenating file location and tree file name
my $a_whole_path = "$a_tree_in[0]";
my $b_whole_path = "$b_tree_in[0]";

print $logMolClockDebug "a_whole_path: $a_whole_path\n";
print $logMolClockDebug "b_whole_path: $b_whole_path\n";
#creating A tree for use
my $treeio  = Bio::TreeIO->new( -format => 'newick', 
								-file   => $a_whole_path );

my $a_tree = $treeio->next_tree;

#creating B tree for use
my $treeio2  = Bio::TreeIO->new( -format => 'newick', 
								-file   => $b_whole_path );

my $b_tree = $treeio2->next_tree;

#array of all leaves in each tree
my @a_leaf_nodes = $a_tree->get_leaf_nodes();
my @b_leaf_nodes = $b_tree->get_leaf_nodes();
#numbers of leaves in each tree
my $num_a_leaves = scalar @a_leaf_nodes;
my $num_b_leaves = scalar @b_leaf_nodes;

my $skip_num = 0;

for(my $i = 1; $i <= $num_genes; $i++){
	if($skip_num < scalar(@genes_to_skip)){
		if($i == $genes_to_skip[$skip_num]){
			$skip_num += 1;
			next;
		}
	}
	
	#Creating a folder for each gene, will contain both a and b fasta contained within a and b folders
	my $newdir = $output_location.$output_gene_root.$i;
	
	print $logMolClockDebug "newDir: $newdir\n\n";
	mkdir $newdir, 0755;
	
	#Creating an A and B folder to go in outside folder. A will contain a gene and same for b
	my $newdir_A = $output_location.$output_gene_root.$i."/A_tree";
	mkdir $newdir_A, 0755;
	my $newdir_B = $output_location.$output_gene_root.$i."/B_tree";
	mkdir $newdir_B, 0755;
	
	#creating and opening a .fasta for A and B within gene folder within respective a and b folders
	my $a_fasta = $newdir_A."/".$root_fasta_name."A_gene".$i.".fasta";
	my $b_fasta = $newdir_B."/".$root_fasta_name."B_gene".$i.".fasta";
	open(my $fh_a, ">", $a_fasta) or die "Could not open file '$a_fasta'";
	open(my $fh_b, ">", $b_fasta) or die "Could not open file '$b_fasta'";
	
	#prints number of species from A and B into respective files
	print $fh_a $num_a_leaves, " ";
	print $fh_b $num_b_leaves, " ";
	
	
	#opening .fasta file from CleanedGap
	my $fasta_open_name = $og_locations.$root_fasta_name.$i.".fasta";
	open(my $fh_fasta, "<", $fasta_open_name) or die "Could not open file '$fasta_open_name'";
	
	#creates an array with all lines from file
	chomp(my @lines = <$fh_fasta>);
	my $tot_num_arr = scalar @lines;
	
	#Append # species + Gene Length to postion [0] of @lines
	
	
	
	#first line is the num species and seq length
	my $row1 = $lines[0];
	my @nums = split / /, $row1;
	print $fh_a $nums[2],"\n";
	print $fh_b $nums[2],"\n";
	
	my $current_a_leaf;
	# index of species name
	my @spec_and_seq;
	#iterating through a-tree leaves
	for(my $a = 0; $a < $num_a_leaves; $a++){
		$current_a_leaf = $a_leaf_nodes[$a]->id();
		chomp($current_a_leaf);
		##two spaces between species and sequence otherwise MCMCtree won't recognize the sequence
		print $fh_a $current_a_leaf,"  ";
		for (my $ind = 0; $ind < $tot_num_arr; $ind++){
			if (index($lines[$ind], $current_a_leaf) != -1){
				@spec_and_seq = split / /, $lines[$ind];
				
				print $fh_a $spec_and_seq[1], "\n";
			}
		}
		
	}
	#iterating through b-tree leaves
	my $current_b_leaf;
	for(my $b = 0; $b < $num_b_leaves; $b++){
		$current_b_leaf = $b_leaf_nodes[$b]->id();
		chomp($current_b_leaf);
		print $fh_b $current_b_leaf,"  ";
		for (my $ind = 0; $ind < $tot_num_arr; $ind++){
			if (index($lines[$ind], $current_b_leaf) != -1){
				@spec_and_seq = split / /, $lines[$ind];
				
				print $fh_b $spec_and_seq[1], "\n";
			}
		}
		
	}
	
	close $fh_fasta;
	
	close $fh_a;
	close $fh_b;
	print $logMolClockDebug "A-tree Gene$i folder created\n";
	print $logMolClockDebug "B-tree Gene$i folder created\n";
}






