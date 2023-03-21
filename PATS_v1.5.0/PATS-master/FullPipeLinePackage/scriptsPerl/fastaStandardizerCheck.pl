#!/usr/bin/perl
  
use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Copy;


my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $fastaDir = "$baseDir/fasta";



#generate logs
open my $logCheckFH, '>', "$baseDir/DataStandardizerLogs/checkLog.txt" or die "Cannot open checkLog.txt: $!";
open my $goldNumFH, '>', "$baseDir/DataStandardizerLogs/goldNum.txt" or die "Cannot open goldNum.txt: $!";


opendir(DIR,"$fastaDir") or die "Couldn't open $fastaDir\n";

my @fastaFiles = grep(/\.fasta$/,readdir(DIR));
my $numFasta = scalar(@fastaFiles);

closedir(DIR);


my $itter = 0;
my @fastaArrayTemp;
my $fileStart = '>';
my @fastaCountArray;
my @fastaMaxCountArray;
my @fastaMinCountArray;


foreach (1 .. $numFasta) {

open my $fastaFH, '<', "$fastaDir/$fastaFiles[$itter]" or die "Cannot open fasta file: $fastaFiles[$itter]: $!";
while(<$fastaFH>) {
push @fastaArrayTemp, $_;
}

close $fastaFH;


my @countArrayTemp;

foreach my $lineSearch2 (@fastaArrayTemp) {
if ($lineSearch2 =~ /^$fileStart/) {
my $fastaLength = length($lineSearch2);
push @countArrayTemp, "$fastaLength\n";
} else {

}
}

my @sortedCountArrayTemp = sort { $a <=> $b } @countArrayTemp;
my $maxCharCount = $sortedCountArrayTemp[-1];
my $minCharCount = $sortedCountArrayTemp[0];

push @fastaMaxCountArray, $maxCharCount;
push @fastaMinCountArray, $minCharCount;

undef @fastaArrayTemp;

$itter ++;

}

my @sortedMaxArray = sort { $a <=> $b } @fastaMaxCountArray;
my $goldenNumberMax = $sortedMaxArray[-1];

print $logCheckFH "The maximum number of characters in the acession line across all .fasta files is:\n";
print $logCheckFH "$goldenNumberMax\n\n";

my @sortedMinArray = sort { $a <=> $b } @fastaMinCountArray;
my $goldenNumberMin = $sortedMinArray[0];

print $logCheckFH "The minimum number of characters in the acession line across all .fasta files is:\n";
print $logCheckFH "$goldenNumberMin\n";


if ($goldenNumberMax == $goldenNumberMin) {

print $logCheckFH "You are good to run the pipeline!\n";
print "\n";
print "Fasta files are ready for the pipeline!!\n";
print "All acession lines are equal to $goldenNumberMax in length\n";
print "\n";

print $goldNumFH "$goldenNumberMax";

} else {

print $logCheckFH "ERROR: Fasta lengths are NOT equal!!\n";

print "\n";
die "ERROR: Fasta lengths are NOT equal!!\n";
print "\n";

}


