#!/bin/perl.exe
use strict;
use warnings;

use Parallel::ForkManager;
use Data::Dumper qw(Dumper);
use File::Copy;


use lib "/mnt/gs18/scratch/users/u5g6jsxd/Scratch_GS18_MEGA/zeroBalanced/master_script/";
use MyModule;


#Number of Genes
my $numGenes = $splitarray[22];
chomp($numGenes);

#Genes to be skipped, seperated by spaces
my $skip_genes = $splitarray[23];
chomp($skip_genes);

# SPLITS MULTIPLE GENES TO BE SKIPPED INTO ARRAY 
my @genes_to_skip = split(/ /, $skip_genes);

#Path of MEGA .mao file
my $mao_MEGA_loc = $splitarray[33];
chomp($mao_MEGA_loc);

#Path of calibration files
my $MEGA_cali_loc = $splitarray[34];
chomp($MEGA_cali_loc);

#Path of OutGroup file
my $MEGA_outgroup_loc = $splitarray[35];
chomp($MEGA_outgroup_loc);

#Path of Topology NWK
my $MEGA_topology_loc = $splitarray[36];
chomp($MEGA_topology_loc);

my $splitGeneDir = "/mnt/gs18/scratch/users/u5g6jsxd/Scratch_GS18_MEGA/zeroBalanced/master_script/split_genes";

#generate log
open my $logfh, '>', "/mnt/gs18/scratch/users/u5g6jsxd/Scratch_GS18_MEGA/zeroBalanced/master_script/MegaLog.txt" or die "Cannot open MegaLog.txt: $!";


my $index = 1;

for (1 .. $numGenes) {

if (grep {$_ eq $index} @genes_to_skip){

$index++;

}else{

chdir ("$splitGeneDir/765phyl_Gene_$index/A_tree"  or die "Can't change dir: $!");

my $timestamp1 = localtime(time);
print $logfh "765phyl_Gene_$index Start: $timestamp1\n";


system("$splitarray[30] -a ${mao_MEGA_loc} -d geneAR_765phyl_A_gene$index.meg -t ${MEGA_topology_loc} -g ${MEGA_outgroup_loc} -c ${MEGA_cali_loc} -o outfile_root");

my $timestamp2 = localtime(time);
print $logfh "765phyl_Gene_$index End: $timestamp2\n";

$index++;

}
}


close $logfh;
