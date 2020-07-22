#!/usr/bin/perl
  
use strict;
use warnings;

use File::Copy;
use File::Path;
use File::Find;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);



my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $fastaDir = "$baseDir/fasta";


#Create directory
my $directoryfBU = "$fastaDir/fastaOriginal";
unless(-d $directoryfBU)
{
mkdir $directoryfBU;
}

my $directoryFileStand = "$baseDir/DataStandardizerLogs";
unless(-d $directoryFileStand)
{
mkdir $directoryFileStand;
}

my $directoryFastaEdit = "$fastaDir/fastaEdit";
unless(-d $directoryFastaEdit)
{
mkdir $directoryFastaEdit;
}

my $directoryFastaEdit2 = "$fastaDir/fastaEdit2";
unless(-d $directoryFastaEdit2)
{
mkdir $directoryFastaEdit2;
}

my $directoryFastaEdit3 = "$fastaDir/fastaChangedBackUp";
unless(-d $directoryFastaEdit3)
{
mkdir $directoryFastaEdit3;
}


my $fastaBUDir = "$fastaDir/fastaOriginal";
my $fastaEditDir = "$fastaDir/fastaEdit";
my $fastaEditDir2 = "$fastaDir/fastaEdit2";
my $fastaEditDir3 = "$fastaDir/fastaChangedBackUp";


#generate logs
open my $logfastaFH, '>', "$baseDir/DataStandardizerLogs/fastaStandardLog.txt" or die "Cannot open fastaStandardLog.txt: $!";
open my $logfastaCountFH, '>', "$baseDir/DataStandardizerLogs/fastaCountLog.txt" or die "Cannot open fastaCountLog.txt: $!";
open my $logfastaDataFH, '>', "$baseDir/DataStandardizerLogs/fastaDataLog.txt" or die "Cannot open fastaDataLog.txt: $!";


opendir(DIR,"$fastaDir") or die "Couldn't open $fastaDir\n";

my @fastaFiles = grep(/\.fasta$/,readdir(DIR));
my $numFasta = scalar(@fastaFiles);




foreach my $filename (@fastaFiles) {
print $logfastaFH "$filename\n";
}

closedir(DIR);


my $itter = 0;
my @fastaArrayTemp;
my $fileStart = '>';
my @fastaCountArray;
my @fastaFinalCountArray;



foreach (1 .. $numFasta) {

open my $fastaFH, '<', "$fastaDir/$fastaFiles[$itter]" or die "Cannot open fasta file: $fastaFiles[$itter]: $!";
while(<$fastaFH>) {
push @fastaArrayTemp, $_;
}

close $fastaFH;

if (-f "$fastaBUDir/$fastaFiles[$itter].bak")
{
print $logfastaFH "The file: $fastaFiles[$itter].bak already exists in this directory\n";
}else{
fmove("$fastaDir/$fastaFiles[$itter].bak","$fastaBUDir/$fastaFiles[$itter].bak");
}


open my $fastaChangeFH, '>', "$fastaEditDir/$fastaFiles[$itter]" or die "Cannot open $fastaFiles[$itter]: $!";


foreach my $lineSearch (@fastaArrayTemp) {
if ($lineSearch =~ /^$fileStart/) {
for ($lineSearch){
s/$fileStart/$fileStart$fastaFiles[$itter]|*|/g;
}
} else {
}
}

print $fastaChangeFH @fastaArrayTemp;
close $fastaChangeFH;


print $logfastaDataFH "Fasta File: $fileStart$fastaFiles[$itter]\n";
print $logfastaDataFH "Accession Sizes (number of characters): \n";

my @countArrayTemp;

foreach my $lineSearch2 (@fastaArrayTemp) {
if ($lineSearch2 =~ /^$fileStart/) {
my $fastaLength = length($lineSearch2);
my $fastaLengthFix = length($lineSearch2) - 1;
push @countArrayTemp, "$fastaLength\n";
print $logfastaDataFH "$fastaLengthFix\n";
} else {

}
}

my @sortedCountArrayTemp = sort { $a <=> $b } @countArrayTemp;
my $maxCharCount = $sortedCountArrayTemp[-1];

my $maxCharCountFix = $maxCharCount - 1;

print $logfastaCountFH "Fasta File: $fastaFiles[$itter]\n";
print $logfastaCountFH "$maxCharCountFix\n";

push @fastaFinalCountArray, $maxCharCount;


undef @fastaArrayTemp;

$itter ++;

}

my @sortedMaxArray = sort { $a <=> $b } @fastaFinalCountArray;
my $goldenNumber = $sortedMaxArray[-1];

my $goldenNumberFix = $goldenNumber - 1;


print $logfastaCountFH "\n\n\n";
print $logfastaCountFH "The longest acession string is:\n";
print $logfastaCountFH $goldenNumberFix;


#Reopen newly created fasta edit files and adjust accession length

my @fastaEditArrayTemp;

my @charPlaceHold = ( 'X' );

opendir(DIR2,"$fastaEditDir") or die "Couldn't open $fastaEditDir\n";

my @fastaEditFiles = grep(/\.fasta/,readdir(DIR2));
my $numEditFasta = scalar(@fastaEditFiles);

closedir(DIR2);

my $itter2 = 0;

foreach (1 .. $numEditFasta) {

open my $fastaFH2, '<', "$fastaEditDir/$fastaEditFiles[$itter2]" or die "Cannot open fasta file: $fastaEditFiles[$itter2]: $!";
while(<$fastaFH2>) {
push @fastaEditArrayTemp, $_;
}
close $fastaFH2;

foreach my $lineEditSearch (@fastaEditArrayTemp) {

if ($lineEditSearch =~ /^$fileStart/) {

for ($lineEditSearch){

my $addNumOfChar = $goldenNumber - length($lineEditSearch);
my $stringSize = join '' => map $charPlaceHold[rand @charPlaceHold], 1 .. $addNumOfChar;

s/\*/$stringSize/g;

}
} else {

}
}


open my $fastaChangeFH2, '>', "$fastaEditDir2/$fastaFiles[$itter2]" or die "Cannot open $fastaFiles[$itter2]: $!";
print $fastaChangeFH2 @fastaEditArrayTemp;
close $fastaChangeFH2;

fcopy("$fastaEditDir2/$fastaFiles[$itter2]","$fastaEditDir3/$fastaFiles[$itter2]");

fmove("$fastaEditDir2/$fastaFiles[$itter2]","$fastaDir/$fastaFiles[$itter2]");

undef @fastaEditArrayTemp;

$itter2 ++;
}



#Clean Up


rmtree $fastaEditDir;
rmtree $fastaEditDir2;









