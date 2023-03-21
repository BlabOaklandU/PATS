#!/bin/perl.exe
use strict;
use warnings;

use Parallel::ForkManager;
use Data::Dumper qw(Dumper);
use File::Copy;

use lib "/mnt/gs18/scratch/users/u5g6jsxd/Scratch_GS18/zeroBalanced/mcmctreeProject/master_script";
use MyModule;

#Number of Genes
my $numGenes = $splitarray[22];
chomp($numGenes);

#Genes to be skipped, seperated by spaces
my $skip_genes = $splitarray[23];
chomp($skip_genes);

# SPLITS MULTIPLE GENES TO BE SKIPPED INTO ARRAY 
my @genes_to_skip = split(/ /, $skip_genes);


my $splitGeneDir = "/mnt/home/u5g6jsxd/mcmctreeProject/master_script/split_genes";


#generate log
#open my $logfh, '>', "$baseDir/LOG.txt" or die "Cannot open LOG.txt: $!";

 
#Run in Parallel -- pushing analyses across moultiple cores

my $numChildProc = 30;
my $pm = Parallel::ForkManager->new( $numChildProc );

my $index = 1;

for (1 .. $numGenes) {

if (grep {$_ eq $index} @genes_to_skip){

$index++;

}else{

chdir ("$splitGeneDir/765phyl_Gene_$index/A_tree"  or die "Can't change dir: $!");
$index++;

$pm->start and next;
system("mcmctree mcmctree1.ctl");
copy("out.BV","in.BV") or die "Copy failed: $!";

#Move mcmctree first run outputs elsewhere

my $mcmc1Run = "mcmctreeFirstRun";
unless(-d $mcmc1Run)
{
mkdir $mcmc1Run or die "Unable to create $mcmc1Run\n";
}
move("/2base.t","/mcmctreeFirstRun/2base.t") or die "Can't move file: $!:";
move("/aaa","/mcmctreeFirstRun/aaa") or die "Can't move file: $!:";
move("/out.BV","/mcmctreeFirstRun/out.BV") or die "Can't move file: $!:";
move("/rst","/mcmctreeFirstRun/rst") or die "Can't move file: $!:";
move("/rst1","/mcmctreeFirstRun/rst1") or die "Can't move file: $!:";
move("/rub","/mcmctreeFirstRun/rub") or die "Can't move file: $!:";
move("/mcmc.txt","/mcmctreeFirstRun/mcmc.txt") or die "Can't move file: $!:";
move("/rst2","/mcmctreeFirstRun/rst2") or die "Can't move file: $!:";
move("/inf","/mcmctreeFirstRun/inf") or die "Can't move file: $!:";
move("/tmp001.ctl","/mcmctreeFirstRun/tmp001.ctl") or die "Can't move file: $!:";
move("/tmp001.out","/mcmctreeFirstRun/tmp001.out") or die "Can't move file: $!:";
move("/tmp001.trees","/mcmctreeFirstRun/tmp001.trees") or die "Can't move file: $!:";
move("/tmp001.txt","/mcmctreeFirstRun/tmp001.txt") or die "Can't move file: $!:";

system("mcmctree mcmctree2.ctl");
$pm->finish;

}
}

$pm->wait_all_children();
$pm->reap_finished_children();