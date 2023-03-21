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

print $logMolClockDebug "\nRunning treeSplitterB.pl... \n\n";

my $directory = $ENV{'splitTreeDir'};
chdir($directory) or die;

########################################
# USER INPUTS
########################################

#Name of tree file input
my $file_in = glob "$ENV{'unsplitMainTreeDir'}/*";

my $directorySplitTreeB = $ENV{'splitTreeBDir'};
chdir($directorySplitTreeB) or die;

#In format >filename.nwk
#Name of output file to generate for B-tree
my $file_outB = $ENV{'outputFileNameTreeB'};

#First Species in B-Tree
my $first_spec_in = $ENV{'firstSpeciesTreeB'};

chomp($file_in);
chomp($file_outB);
chomp($first_spec_in);
#
#
# DONE WITH USER INPUTS
########################################


# parse in newick/new hampshire format
my $treeio  = Bio::TreeIO->new( -format => 'newick', 
								-file   => $file_in);

my $treeout = Bio::TreeIO->new( -format => 'newick',
								-file	=> $file_outB);

my $tree = $treeio->next_tree;

my $first_spec = $tree->find_node($first_spec_in);

#OUTGROUPING
#Name of outgroup species
my $outgroup_in = $ENV{'outgroupSpecies'};
chomp($outgroup_in);

my $outgroup = $tree->find_node($outgroup_in);

#common ancestor between outgroup and first_spec
my $lca_outgroup = $tree->get_lca($first_spec, $outgroup);

#data storage
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

print $logMolClockDebug "Outgroup Species completed. \n\n";

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

print $logMolClockDebug "Making B-Tree\n";
#b tree
my $node_look = $first_spec;
my $current_node = $first_spec->ancestor();
my $branch_len_add = 0;
my $curr_branch_len = 0;
my @node_desc = $current_node->each_Descendent();
my $splice_me;
my $final_len;
while($current_node->id() ne $lca_outgroup->id()){
	@node_desc = $current_node->each_Descendent();
	if( $node_desc[0]->id() eq $node_look->id()){
		$node_look = $current_node;
		$current_node = $current_node->ancestor();
	}else{
		$node_desc[0]->remove_all_Descendents();
		$tree->remove_Node($node_desc[0]);
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

print $logMolClockDebug "B-Tree completed. \n\n";

$treeout->write_tree($tree);