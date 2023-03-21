#!/usr/bin/perl

use strict;
use warnings;
use File::Find;
use File::Copy;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);

my $baseDir = $ENV{'PoPipe_srcBaseDir'};

my $fastTreeSetting = $ENV{'PoPipe_fastTree'};
my $raxmlTreeSetting = $ENV{'PoPipe_RAxMLTree'};

#Make a backup of KeepOne .txt

my $KeepOneFileBU = "$baseDir/permutations/keepOneBU.txt";
my $KeepOneFile = "$baseDir/permutations/keepOne.txt";

if (-e $KeepOneFileBU) {
#Do nothing
} else {
fcopy($KeepOneFile, $KeepOneFileBU) or die ("Unable to rename $KeepOneFile to $KeepOneFileBU.bak: $!\n");
}

#Create temp directory to hold moved files

my $directoryKeepOneRerunTemp = "$baseDir/fasta/temp";
unless(-d $directoryKeepOneRerunTemp)
{
mkdir $directoryKeepOneRerunTemp or die "Unable to create $directoryKeepOneRerunTemp\n";
}

my @keepOneListBU;

open my $keepOneFHBU, '<', "$baseDir/permutations/keepOneBU.txt" or die "Cannot open keepOneBU.txt: $!";
while(<$keepOneFHBU>) { 
    chomp;
    s/\r//g;	
    push @keepOneListBU,$_;
} 
close $keepOneFHBU;

my $numKeepOneBU =  scalar @keepOneListBU;

#Move All Keep One Files to temp folder.

my $index = 0;
my $toDir = "$baseDir/fasta/temp";
foreach (1 .. $numKeepOneBU) {

my $fromDir = "$baseDir/fasta/$keepOneListBU[$index]";


fmove ($fromDir, $toDir) or die  "All keepOne species were not present in fasta directory!\n";

undef $fromDir;

$index ++;
}


my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp;
    s/\r//g;	
    push @keepOneList,$_;
} 
close $keepOneFH;

my $numKeepOne =  scalar @keepOneList;

my $indexTwo = 0;

#Move the Keep One back to run with analysis

my $fromDir2 = "$baseDir/fasta/";
my $toDir2 = "$baseDir/fasta/temp/$keepOneList[$indexTwo]";

fmove ($toDir2, $fromDir2) or die "keepOne species for analysis not present in Temp directory!\n";


undef @keepOneList;
undef @keepOneListBU;