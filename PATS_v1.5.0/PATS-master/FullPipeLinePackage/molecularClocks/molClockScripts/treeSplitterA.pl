#!/bin/perl.exe
use strict;
use warnings;

use Bio::TreeIO;
use Bio::Tree::Node;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::AnnotatableNode;
use List::MoreUtils qw(first_index);
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;



my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $mcBaseDir =  $ENV{'MolecularClockBaseDir'};
my $mcScriptsDir = $ENV{'MolecularClockScriptDir'};
my $mcsplitGenesDir = $ENV{'splitGenesDir'};

open my $logMolClockDebug, '>>', "$mcBaseDir/molClockDebugLog.txt" or die "Cannot open molClockDebugLog.txt: $!";

print $logMolClockDebug "\nRunning treeSplitterA.pl... \n\n";

##########################################################
# 	USER INPUTS 
##########################################################

#Name of tree file input
my $file_in = glob "$ENV{'unsplitMainTreeDir'}/*";

my $directorySplitTreeA = $ENV{'splitTreeADir'};
chdir($directorySplitTreeA) or die;

#In format >filename.nwk
#Name of output file to generate for A-tree
my $file_outA = $ENV{'outputFileNameTreeA'};

#Name of Splitter Species
my $splitter_in = $ENV{'splitterSpecies'};

#Name of Overlap Species 1
my $species_1_in = $ENV{'overlapSpecies1'};

#Name of Overlap Species 2
my $species_2_in = $ENV{'overlapSpecies2'};

chomp($file_in);
chomp($file_outA);
chomp($splitter_in);
chomp($species_1_in);
chomp($species_2_in);

print $logMolClockDebug "In file: $file_in\n\n";
print $logMolClockDebug "In file: $file_outA\n\n";


#
# DONE OBTAINING USER INPUTS
#########################################################

############################################
# CREATING A TREE OUTPUT FILE
############################################
#
#
# parse in newick/new hampshire format

my $treeio  = Bio::TreeIO->new( -format => 'newick', 
								-file   => $file_in);

my $treeout = Bio::TreeIO->new( -format => 'newick',
								-file	=> $file_outA);

my $tree = $treeio->next_tree;


# USING USER INPUTS TO CREATE NODE OBJECTS 
my $splitter = $tree->find_node($splitter_in);
my $species_1 = $tree->find_node($species_1_in);
my $species_2 = $tree->find_node($species_2_in);

print $logMolClockDebug "Generating A-Tree...\n\n";


###########################################
# CARVING OUT OUTGROUP SPECIES 
###########################################
#
#

#Name of outgroup species
my $outgroup_in = $ENV{'outgroupSpecies'};
chomp($outgroup_in);

#CREATING OUTGROUP NODE OBJECT FROM USER INPUT
my $outgroup = $tree->find_node($outgroup_in);

#common ancestor between outgroup and splitter
my $lca_outgroup = $tree->get_lca($splitter, $outgroup);

#STORING DATA: TOTAL BRANCH LENGTH OF OUTGROUP SPECIES UP UNTIL LCA BETWEEN OUTGROUP AND SPLITTER
my $node_to_add_out = $outgroup;
my $tot_branch_len_out = 0;
while($node_to_add_out->id() ne $lca_outgroup->id()){
	$tot_branch_len_out += $node_to_add_out->branch_length();
	$node_to_add_out = $node_to_add_out->ancestor();
}
#descendents of lca
my @node_desc_out = $lca_outgroup->each_Descendent();
$node_desc_out[1]->remove_all_Descendents();
$tree->remove_Node($node_desc_out[1]);

#adding back outgroup
my $add_back_outgroup = Bio::Tree::Node->new( -branch_length => $tot_branch_len_out,
											  -id => $outgroup_in);
$lca_outgroup->add_Descendent($add_back_outgroup);
my $lca_outgroupid = $lca_outgroup->id();
print $logMolClockDebug "LCA OUTGROUP: $lca_outgroupid \n\n";
print $logMolClockDebug "Outgroup Species completed. \n\n";
#
#
# OUTGROUP CREATED
###########################################

#the common node ancestor of species_1 and species_2
my $lca_species = $tree->get_lca($species_1, $species_2);

#DATA STORAGE FOR SPECIES_1 and SPECIES_2 FOR ADDING BACK LATER
my $species_1_id = $species_1->id();
my $species_2_id = $species_2->id();
my $lca_species_id = $lca_species->id();
my $splitter_id = $splitter->id();

print $logMolClockDebug "Common ancestor to $species_1_id and $species_2_id: Node $lca_species_id \n\n";

#species 1 and 2 total branch length up to lca_species
my $node_to_add = $species_1;
my $tot_branch_len_spec = 0;
while($node_to_add->id() ne $lca_species->id()){
	$tot_branch_len_spec += $node_to_add->branch_length();
	$node_to_add = $node_to_add->ancestor();
}

#the least common ancestor of lca_species and splitter
my $lca_split = $tree->get_lca($lca_species, $splitter);
my $lca_split_id = $lca_split->id();
print $logMolClockDebug "Common ancestor to Node $lca_species_id and $splitter_id: Node $lca_split_id \n\n"; 

#total branch length from lca_species to lca_split
my $node_to_add_2 = $lca_species;
my $tot_branch_len_split = 0;
while($node_to_add_2->id() ne $lca_split->id()){
	$tot_branch_len_split += $node_to_add_2->branch_length();
	$node_to_add_2 = $node_to_add_2->ancestor();
}


#array of descendents from lca_split [top_des, bot_des]
my @lca_split_desc = $lca_split->each_Descendent();

#remove all descendents from bottom descendent in lca_split_des 
#because top contains the splitter, so everything under except for
#the two species will not be in A-tree
$lca_split_desc[1]->remove_all_Descendents();

#removes bottom descendent 
$tree->remove_Node($lca_split_desc[1]);

#add lca_spec as descendent to lca_split 
my $add_back_lca_spec = Bio::Tree::Node->new( -branch_length => $tot_branch_len_split,
											  -id => $lca_species_id);
$lca_split->add_Descendent($add_back_lca_spec);


my $add_back_spec1 = Bio::Tree::Node->new( -branch_length => $tot_branch_len_spec,
										   -id => $species_1_id);
										   
$add_back_lca_spec->add_Descendent($add_back_spec1);

my $add_back_spec2 = Bio::Tree::Node->new( -branch_length => $tot_branch_len_spec,
										   -id => $species_2_id);
										   
$add_back_lca_spec->add_Descendent($add_back_spec2);

#####################

my $delete_nodes = $ENV{'nodesToBeDeleted'};
chomp($delete_nodes);
my @delete_nodes_arr = split / /, $delete_nodes;
#iterate through array 
my $del_num = scalar @delete_nodes_arr;

my $del_node2;
my $del_ancs2;
my $splicer2;
my $final_leng2;
my @ancs_desc2;

for (my $i = 0; $i < $del_num; $i++){
	my $del_node = $tree->find_node($delete_nodes_arr[$i]);
	my @delete = $del_node->each_Descendent();
	$delete[1]->remove_all_Descendents();
	$del_ancs2 = $delete[1]->ancestor();
	@ancs_desc2 = $del_ancs2->each_Descendent();
	$splicer2 = $ancs_desc2[0];
	$final_leng2 = $del_ancs2->branch_length() + $splicer2->branch_length();
	$tree->remove_Node($delete[1]);
	$tree->splice($del_ancs2);
	$splicer2->branch_length($final_leng2);
}
##########################################################
#fixing bug where splitter is not last value for A-tree
##########################################################
#
#
my $node_look = $splitter;
my $current_node = $splitter->ancestor();
my $branch_len_add = 0;
my $curr_branch_len = 0;
my @node_desc = $current_node->each_Descendent();
my $splice_me;
my $final_len;
while($current_node->id() ne $lca_outgroup->id()){
	@node_desc = $current_node->each_Descendent();
	if( $node_desc[1]->id() eq $node_look->id() or $node_desc[1]->id() eq $lca_species->id()){
		$node_look = $current_node;
		$current_node = $current_node->ancestor();
	}else{
		$node_desc[1]->remove_all_Descendents();
		$tree->remove_Node($node_desc[1]);
		$splice_me = $current_node;
		#fixing branch lengths because preserve lengths not working
		$branch_len_add = $splice_me->branch_length();
		$curr_branch_len = $node_look->branch_length();
		$current_node = $current_node->ancestor();
		$tree->splice($splice_me);
		$final_len = $branch_len_add + $curr_branch_len;
		$node_look->branch_length($final_len);
	}
}
#
# DONE FIXING BUG
############################################################


if($lca_split_id eq $lca_species_id){
	print $logMolClockDebug "WARNING: Common ancestor node to overlap species is also", 
	" the common ancestor to splitter species. Tree is not longer a binary tree. Program will quit.\n\n";
	exit;
}


$treeout->write_tree($tree);
print $logMolClockDebug "A-Tree completed. \n";
#
#
#DONE CREATING A TREE OUTPUT
#######################################################################




