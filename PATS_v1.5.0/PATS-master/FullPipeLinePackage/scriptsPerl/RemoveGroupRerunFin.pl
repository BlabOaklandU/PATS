#!/usr/bin/perl

use strict;
use warnings;
use File::Path;
use File::Find;
use File::Copy;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);




my $baseDir = $ENV{'PoPipe_srcBaseDir'};

my $fastTreeSetting = $ENV{'PoPipe_fastTree'};
my $raxmlTreeSetting = $ENV{'PoPipe_RAxMLTree'};

my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp;
	s/\r//g;
    push @removeGroupList,$_;
} 

close $removeGroupFH;

my $numRemoveGroup =  scalar @removeGroupList;


#Move All Remove Group Files to fasta folder.

my $index = 0;
foreach (1 .. $numRemoveGroup) {

my $fromDir = "$baseDir/fasta/temp/$removeGroupList[$index]";
my $toDir = "$baseDir/fasta";

fmove ($fromDir, $toDir) or die  "All removeGroup species were not present in directory!\n";

$index ++;
}
