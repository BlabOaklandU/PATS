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

my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp;
	s/\r//g;
    push @removeOneList,$_;
} 

close $removeOneFH;


#Move the remove One file back to fasta folder.

my $index = 0;

my $fromDir = "$baseDir/fasta/temp/$removeOneList[$index]";
my $toDir = "$baseDir/fasta";

fmove ($fromDir, $toDir) or die  "All keepOne species were not present in directory!\n";


my $filename = "$baseDir/permutations/removeOne.txt";

if (-z $filename) {

#Nothing

} else {

#Open up KeepOne list and remove index 0

rename($filename, $filename . '.bak') or die ("Unable to rename $filename to $filename.bak: $!\n");

open(INFILE, $filename . '.bak') or die("Can't open $filename.bak for input: $!\n");
open(OUTFILE, '>' . $filename) or die("Can't open $filename for output: $!\n");
  
  my $linecount = 0;
  while (my $line = <INFILE>) { 
    print(OUTFILE $line) if ($linecount);
    $linecount++;
  }
close(INFILE);
close(OUTFILE);

my $filenameDel = "$baseDir/permutations/removeOne.txt.bak";
unlink $filenameDel;

my $directoryFastaTemp = "$baseDir/fasta/temp";
rmtree $directoryFastaTemp;

}











