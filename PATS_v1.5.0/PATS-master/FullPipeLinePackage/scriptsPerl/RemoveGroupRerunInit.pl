#!/usr/bin/perl

use strict;
use warnings;
use File::Find;
use File::Copy;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);

my $baseDir = $ENV{'PoPipe_srcBaseDir'};

my $fastTreeSetting = $ENV{'PoPipe_fastTree'};
my $raxmlTreeSetting = $ENV{'PoPipe_RAxMLTree'};



#Create temp directory to hold moved files

my $directoryRemoveGroupRerunTemp = "$baseDir/fasta/temp";
unless(-d $directoryRemoveGroupRerunTemp)
{
mkdir $directoryRemoveGroupRerunTemp or die "Unable to create $directoryRemoveGroupRerunTemp\n";
}

my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp; 
	 s/\r//g;
    push @removeGroupList,$_;
} 
close $removeGroupFH;

my $numRemoveGroup =  scalar @removeGroupList;



#Move All Remove Group Files to temp folder.

my $index = 0;
foreach (1 .. $numRemoveGroup) {

my $fromDir = "$baseDir/fasta/$removeGroupList[$index]";
my $toDir = "$baseDir/fasta/temp";

fmove ($fromDir, $toDir) or die  "All removeGroup species were not present in directory!\n";

$index ++;
}