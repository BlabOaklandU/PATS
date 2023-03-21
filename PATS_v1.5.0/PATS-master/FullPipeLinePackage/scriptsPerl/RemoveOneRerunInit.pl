#!/usr/bin/perl

use strict;
use warnings;
use File::Find;
use File::Copy;
use File::Path;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);

my $baseDir = $ENV{'PoPipe_srcBaseDir'};

my $fastTreeSetting = $ENV{'PoPipe_fastTree'};
my $raxmlTreeSetting = $ENV{'PoPipe_RAxMLTree'};


#Make a backup of removeOne.txt

my $RemoveOneFileBU = "$baseDir/permutations/removeOneBU.txt";
my $RemoveOneFile = "$baseDir/permutations/removeOne.txt";


if (-e $RemoveOneFileBU) {
#Nothing
} else {
fcopy($RemoveOneFile, $RemoveOneFileBU) or die ("Unable to rename $RemoveOneFile to $RemoveOneFileBU.bak: $!\n");
}

#Create temp directory to hold moved files

my $directoryRemoveOneRerunTemp = "$baseDir/fasta/temp";
unless(-d $directoryRemoveOneRerunTemp)
{
mkdir $directoryRemoveOneRerunTemp or die "Unable to create $directoryRemoveOneRerunTemp\n";
}

my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp;
	s/\r//g;
    push @removeOneList,$_;
} 
close $removeOneFH;

#Move the remove One file to temp folder.

my $fromDir = "$baseDir/fasta/$removeOneList[0]";
my $toDir = "$baseDir/fasta/temp";

fmove ($fromDir, $toDir) or die  "All removeOne species were not present in directory!\n";

