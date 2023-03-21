#!/usr/bin/perl

use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);

my $baseDir = $ENV{'PoPipe_srcBaseDir'};


my $directory3 = "$baseDir/OGGapFilter";
unless(-d $directory3)
{
mkdir $directory3 or die "Unable to create $directory3\n";
}

#Create base RAxML Tree directory
my $directoryBaseTreeRAxML = "$baseDir/treeRAxML";
unless(-d $directoryBaseTreeRAxML)
{
mkdir $directoryBaseTreeRAxML or die "Unable to create $directoryBaseTreeRAxML\n";
}

#Create base fastTree Tree directory
my $directoryBaseTreefastTree = "$baseDir/fastTree";
unless(-d $directoryBaseTreefastTree)
{
mkdir $directoryBaseTreefastTree or die "Unable to create $directoryBaseTreefastTree\n";
}


#generate log
open my $gapLogFH, '>', "$baseDir/logFiles/GapFilterLog.txt" or die "Cannot open GapFilterLog.txt: $!";

#Generate gap count index log
open my $gapIndexLogFH, '>', "$baseDir/logFiles/GapFilterIndexLog.txt" or die "Cannot open GapFilterIndexLog.txt: $!";

############################################################
#Create .txt file that stores the order based on species ID#
############################################################

my $fastaDir = "$baseDir/fasta/";
opendir my $dh, $fastaDir  or die "Can't open $fastaDir: $!";
my @numSpeicesFiles = grep(/\.fasta$/,readdir($dh));
my $numSpecies =  scalar(@numSpeicesFiles);

my $OGDir = "$baseDir/muscleFasta";
opendir my $OGdirfh, $OGDir  or die "Can't open $OGDir: $!";
my $numOGFiles =  grep { -f "$OGDir/$_" } readdir($OGdirfh);

print $gapLogFH "Number of Species in analysis = $numSpecies\n\n";

#######################
#Read in first OG file#
#######################

my @muscleOutputArrayTemp;
my @muscleOutputArray;
my $fileItter = 0;

#################
#Threshold input#
#################

my $gapThreshholdInput = $ENV{'PoPipe_gapFilterCutoff'};
chomp($gapThreshholdInput);

print $gapLogFH "User Threshold percentage input = $gapThreshholdInput %\n";


my $fileStart = ">";

foreach (1 .. $numOGFiles){

open my $muscleFH, '<', "$baseDir/muscleFasta/OG$fileItter.fasta" or die "Cannot open OG$fileItter.fasta: $!";
chomp(@muscleOutputArrayTemp = <$muscleFH>);
close $muscleFH;

foreach my $muscleLine (@muscleOutputArrayTemp) {
if ($muscleLine =~ m/^$fileStart/){
#do nothing
} else {
$muscleLine =~ s/\n//g;
}
}

my $fileStartReplace = "<>";

foreach my $muscleLine (@muscleOutputArrayTemp) {
if ($muscleLine =~ m/^$fileStart/){
} else {
$muscleLine = $muscleLine.'^';
}
}

foreach my $muscleLine2 (@muscleOutputArrayTemp) {
if ($muscleLine2 =~ m/^$fileStart/){
$muscleLine2 =~ s/$fileStart/$fileStartReplace/;
} else {
}
}

my $largeStringMuscle = join ("^",@muscleOutputArrayTemp);

my @newmuscleOutputArrayTemp = split(/</,$largeStringMuscle);
splice(@newmuscleOutputArrayTemp, 0, 1);

foreach my $muscleLine3 (@newmuscleOutputArrayTemp) {
$muscleLine3 =~ s/\^/ /;
}

foreach my $muscleLine4 (@newmuscleOutputArrayTemp) {
$muscleLine4 =~ s/\^//g;
}

@muscleOutputArray = @newmuscleOutputArrayTemp;



my $numSeq =  scalar @muscleOutputArray;

my $gapThreshold = ($gapThreshholdInput/100) * $numSeq;
print $gapLogFH "Maximum number of gaps allowed per alignment = $gapThreshold\n";

my $count = 0;
my @gapArrayTest;
my @gapArray;

foreach (1 .. $numSeq) {

my $fastaString = $muscleOutputArray[$count];
my $gap = '-';
my $offset = 0;

my $result = index($fastaString, $gap, $offset);

my $seqLengthTest = length($fastaString);

while ($result != -1) {

push @gapArrayTest, "$result\n";

$offset = $result + 1;
$result = index($fastaString, $gap, $offset);
}
$count++;
}

my @gapCountArray;
my %countGap;
my $itt = 0;

foreach my $seqElement( @gapArrayTest ) {
  ++$countGap{$seqElement};
}

print $gapIndexLogFH "\nOG$fileItter:\n";

foreach my $seqElement( keys %countGap ) {
print $gapIndexLogFH "$seqElement = $countGap{$seqElement}\n";
}


foreach my $seqElement( keys %countGap ) {

if ($countGap{$seqElement} > $gapThreshold) {
push @gapCountArray, $seqElement;
}
$itt ++;
}

my @sortedGapCountArray = sort { $b <=> $a } @gapCountArray;

if (@sortedGapCountArray == 0){
print $gapLogFH "OG$fileItter: Location of Gaps to be removed: no gaps were removed\n";
} else {
print $gapLogFH "OG$fileItter: Location of Gaps to be removed:\n @sortedGapCountArray \n";
}


#Delete gaps based on threshold. Delete from end to front to avoid shifting

my @gapFreeArrayFinal;
my $gapRemovalItter = 0;
my $itterSeq = 0;

foreach (0 .. $numSeq-1) {

foreach my $seqElement2(@sortedGapCountArray) {

substr($muscleOutputArray[$itterSeq], $seqElement2, 1, "");

}

push @gapFreeArrayFinal, "$muscleOutputArray[$itterSeq]\n";
#Check this is true!
$itterSeq ++;
}

for (@gapFreeArrayFinal) {
   s/\s+/\n/g;
}


open my $singleFilePushGapFree, '>', "$baseDir/OGGapFilter/OG$fileItter.fasta" or die "Cannot open OG$fileItter.fasta: $!";

print $singleFilePushGapFree @gapFreeArrayFinal;

close $singleFilePushGapFree;

undef @newmuscleOutputArrayTemp;
undef @gapFreeArrayFinal;
undef @sortedGapCountArray;
undef @gapCountArray;
undef @gapArrayTest;
undef @muscleOutputArrayTemp;
undef @muscleOutputArray;

$fileItter++;

}



#Grab acession length

open my $gnFH, '<', "$baseDir/DataStandardizerLogs/checkLog.txt" or die "cannot open input file:$!";
while (my $lineFinder = <$gnFH>) {
if ($. == 2) {
my $GoldenNumber = $lineFinder;
}
}

my $directoryNewGap = "$baseDir/cleanedGap";
unless(-d $directoryNewGap)
{
mkdir $directoryNewGap or die "Unable to create $directoryNewGap\n";
}


#Open newly created OG Files (OGGapFilter) and append missing species with just gaps

my @missingGapsSpeciesArray;
my @dataStoreIndexArray;

my $fileItter2 = 0;

foreach (1 .. $numOGFiles){
my $itterSpecies = 0;
open my $addGapsMissingSpecies, '<', "$baseDir/OGGapFilter/OG$fileItter2.fasta" or die "Cannot open OG$fileItter2.fasta: $!";

while(<$addGapsMissingSpecies>) { 
    push @missingGapsSpeciesArray, $_;
} 
close $addGapsMissingSpecies;

my $speciesOGCount = grep{ />/} @missingGapsSpeciesArray;


#Get num of lines + character counts

my @accesIndex = grep { $missingGapsSpeciesArray[$_] =~ />/ } 0..$#missingGapsSpeciesArray;

@dataStoreIndexArray = @missingGapsSpeciesArray[$accesIndex[0] + 1 .. $accesIndex[1] - 1];
foreach (@dataStoreIndexArray) {
s/[a-z]/-/ig;
}

my @allSpeciesArray;

open my $speciesListHandle, '<', "$baseDir/logFiles/speciesPresent.txt" or die "Cannot open speciesPresent.txt: $!";
while(<$speciesListHandle>) { 
    chomp; 
    push @allSpeciesArray,$_;
} 
close $speciesListHandle;

my $speciesString = $allSpeciesArray[0];
my @allSpeciesOrderArray = split ' ', $speciesString;

my $numSpecies = scalar(grep {defined $_} @allSpeciesOrderArray),;


foreach (1 .. $numSpecies){
	
if ( grep( m/$allSpeciesOrderArray[$itterSpecies]/, @missingGapsSpeciesArray ) ) {
#Do nothing
$itterSpecies++;
} else {
push @missingGapsSpeciesArray, ">$allSpeciesOrderArray[$itterSpecies]|SpecieNotPresentInOrtholog|\n";
push @missingGapsSpeciesArray, @dataStoreIndexArray;
$itterSpecies++;
}
}

open my $newOFFile, '>', "$baseDir/cleanedGap/OG$fileItter2.fasta" or die "Cannot open OG$fileItter2.fasta: $!";

print $newOFFile @missingGapsSpeciesArray;
close $newOFFile;

$fileItter2++;

undef @missingGapsSpeciesArray;
undef @dataStoreIndexArray;
undef $itterSpecies;
}
