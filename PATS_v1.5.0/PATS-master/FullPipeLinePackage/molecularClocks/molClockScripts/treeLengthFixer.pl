#!/bin/perl.exe
use strict;
use warnings;

use Bio::TreeIO;
use Bio::Tree::Node;

################################
# USER INPUTS				   #
################################
#
#
print "Tree File Name: ";
my $file_in = <STDIN>;
chomp($file_in);

print "Output File Name: ";
my $file_out = <STDIN>;
chomp($file_out);

print "True Length (Root to Tip) of tree: ";
my $true_len_in = <STDIN>;
chomp ($true_len_in);
#
#
#END USER INPUTS
################################

##################################
# CREATE TREE OBJECT			 #
##################################
#
#
# parse in newick/new hampshire format
my $input = new Bio::TreeIO(-file   => $file_in,
                            -format => "newick");
my $tree = $input->next_tree;
#
# TREE OBJECT CREATED
###################################


#####################################
# ANALYZE ROOT TO TIP LENGTHS FOR   #
# ALL LEAVES IN THE TREE, PRINT     #
# LEAVES THAT NEED TO BE ADJUSTED   #
# TO A TXT FILE                     #
#####################################
#
#
#makes an array of all the leaves in the tree
my @leaves = $tree->get_leaf_nodes();

#opens file to write to
my $filename = $file_out;
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";

#prints number of leaves in the tree
print $fh scalar @leaves, " leaves in the tree\n\n";
#stores number of leaves in the tree
my $leaves_num = scalar @leaves;

#what the root to tip lengths SHOULD be
my $true_len = $true_len_in;
print $fh "The True Length is ", $true_len, "\n\n";

# a count of how many errored leaves there are
my $count_wrong = 0;

print $fh "Leaves that do not match true times:\n";
#iterates through all the leaves stored in leaves array 
for(my $i= 0; $i < $leaves_num; $i++){
	#retrieves and stores the root to tip length of leaf
	my $depth = $leaves[$i]->depth;
	
	#prints if the lengths are not the same
	if($depth != $true_len){
		print $fh $leaves[$i]->id,"\t";
		print $fh $depth;
		print $fh "\tDifference: ",($true_len - $depth), "\n";
		$count_wrong++;
	}
}

print $fh "Number of Errors: ", $count_wrong;
close $fh;
print "done\n";
#
#
# DONE ANALYZING AND PRINTING 
#####################################