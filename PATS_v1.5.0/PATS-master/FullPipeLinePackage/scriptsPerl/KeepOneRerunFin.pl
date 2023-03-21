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

my @keepOneListBU;

open my $keepOneFH, '<', "$baseDir/permutations/keepOneBU.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp; 
	s/\r//g;
    push @keepOneListBU,$_;
} 

close $keepOneFH;

my $numKeepOne =  scalar @keepOneListBU;


#Move All Keep One Files to fasta folder.

my $index = 0;
foreach (1 .. $numKeepOne) {

my $fromDir = "$baseDir/fasta/temp/$keepOneListBU[$index]";
my $toDir = "$baseDir/fasta";

fmove ($fromDir, $toDir);

$index ++;
}

my $filename = "$baseDir/permutations/keepOne.txt";

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

my $filenameDel = "$baseDir/permutations/keepOne.txt.bak";
unlink $filenameDel;


my $directoryFastaTemp = "$baseDir/fasta/temp";
rmtree $directoryFastaTemp;

}

undef @keepOneListBU;







