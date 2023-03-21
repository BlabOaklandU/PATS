#!/usr/bin/perl
  
use strict;
use warnings;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};


#create og directory
my $directoryOGstore = "$baseDir/og";

unless(-d $directoryOGstore)
{
mkdir $directoryOGstore or die "Unable to create $directoryOGstore\n";
}

#create log directory
my $directoryLogs = "$baseDir/logFiles";

unless(-d $directoryLogs)
{
mkdir $directoryLogs or die "Unable to create $directoryLogs\n";
}


#Read in proteinOrtho output into array
my $myOrthoParaPO = "$baseDir/orthologs/myproject.proteinortho.tsv";

my @orthoParaArray;
open(my $fh, "<", $myOrthoParaPO)
    or die "Failed to open 'myproject.proteinortho.tsv' file: $!\n";
while(<$fh>) { 
    chomp; 
    push @orthoParaArray, $_;
} 
close $fh;

#generate log
open my $logFH, '>', "$baseDir/logFiles/speciesFilterLog.txt" or die "Cannot open Log.txt: $!";

#Threshold
my $threshholdInput = $ENV{'PoPipe_speciesFilterCutoff'};
chomp($threshholdInput);

my $alignConnThresholdInput = $ENV{'PoPipe_proteinOrthoConn'};

print $logFH "Threshold input = $threshholdInput %\n";
print $logFH "Alignment connectivity threshold input = $alignConnThresholdInput %\n";

#Number of species count and push species to .txt
my $numSpeciesString = $orthoParaArray[0];
my @numSpeciesArray = split ' ', $numSpeciesString;
my $numSpecies = scalar(@numSpeciesArray) - 4;
my $numSpecies2 = scalar(@numSpeciesArray);
my $numSpecies3 = $numSpecies2 - 1;

my @speciesList = @numSpeciesArray[4 .. $numSpecies3];

my @speciesList2;
my $speciesString = join('  ',@speciesList);
push @speciesList2, $speciesString;

open my $speciesFH, '>', "$baseDir/logFiles/SpeciesList.txt" or die "Cannot open SpeciesList.txt: $!";
print $speciesFH @speciesList2;

print $logFH "Number of species in analysis = $numSpecies\n";

my $threshold = $numSpecies * ($threshholdInput/100);

print $logFH "Actual threshold number of species = $threshold\n";


#Duplicate species check - Removal
my $arrayLength = scalar(grep {defined $_} @orthoParaArray),;
my $counter = 1;
my @orthoArray;
my @paraArray;
my @coOrthoArray;
my $arrayString;

foreach (1 .. $arrayLength -1) {

my $orthoString = $orthoParaArray[$counter];
my @tempArray = split ' ', $orthoString;

if ($tempArray[0] >= $threshold && $tempArray[0] == $tempArray[1]) 
{ 
my $tempArrayLength = scalar(grep {defined $_} @tempArray);
$arrayString = join(' ',@tempArray);
push @orthoArray, $arrayString;
} elsif ($tempArray[0] < $threshold && $tempArray[0] < $tempArray[1]) {
my $tempArrayLength = scalar(grep {defined $_} @tempArray);
$arrayString = join(' ',@tempArray);
push @paraArray, $arrayString;
} elsif ($tempArray[0] >= $threshold && $tempArray[0] < $tempArray[1]) {
my $tempArrayLength = scalar(grep {defined $_} @tempArray);
$arrayString = join(' ',@tempArray);
push @coOrthoArray, $arrayString;
}else{
}
$counter++;
}

print $logFH "Total number of combined groups = $arrayLength\n";

my $numOG = scalar(grep {defined $_} @orthoArray),;
print $logFH "Number of stored orthologous groups = $numOG\n";
my $numcoOG = scalar(grep {defined $_} @coOrthoArray),;
print $logFH "Number of stored co-orthologous groups = $numcoOG\n";
my $numPG = scalar(grep {defined $_} @paraArray),;
print $logFH "Number of stored paralogous groups = $numPG\n";

my $numRemoved = $arrayLength - ($numOG + $numPG + $numcoOG);
print $logFH "Number of removed groups = $numRemoved\n";


#Standardize .txt files with equal length array elements
my $itterAst = 0;
my $GoldenNumber;

open my $gnFH, '<', "$baseDir/DataStandardizerLogs/checkLog.txt" or die "cannot open input file:$!";

while (my $lineFinder = <$gnFH>) {
    if ($. == 2) {
      $GoldenNumber = $lineFinder;
    }
}

my $adjustedGoldenNumber = $GoldenNumber - 3;

my $replace = "*" x $adjustedGoldenNumber;

foreach (1 .. $numOG) {
$orthoArray[$itterAst] =~ s/\*/$replace/g;
$itterAst++;
}


open $fh, '>', "$baseDir/logFiles/orthoOutput.txt" or die "Cannot open orthoOutput.txt: $!";
foreach (@orthoArray)
{
print $fh "$_\n";
}
close $fh;

open $fh, '>', "$baseDir/logFiles/coOrthoOutput.txt" or die "Cannot open coOrthoOutput.txt: $!";
foreach (@coOrthoArray)
{
print $fh "$_\n";
}
close $fh;

open $fh, '>', "$baseDir/logFiles/paraOutput.txt" or die "Cannot open paraOutput.txt: $!";
foreach (@paraArray)
{
print $fh "$_\n";
}
close $fh;

#Push each ortho group into its own text file

my $itter = 0;

foreach (1 .. $numOG) {

open my $OGfh, '>', "$baseDir/og/OG$itter.txt" or die "Cannot open OG$itter.txt: $!";

print $OGfh $orthoArray[$itter];

close $OGfh;

$itter ++;
}

close $logFH; #close log


#Open OrthoOutput and grab all species present in analysis
open my $speciesPresentFH, '>', "$baseDir/logFiles/speciesPresent.txt" or die "Cannot open speciesPresent.txt: $!";

my @speicesPresentArray;
open my $orthoOutputFH, '<', "$baseDir/logFiles/orthoOutput.txt" or die "Cannot open orthoOutput.txt: $!";

while(my $specieLine = <$orthoOutputFH>){
    chomp $specieLine;
	#s/\*//g;
    my @lineArray = split(" ", $specieLine);
    push(@speicesPresentArray, @lineArray);
}
close $orthoOutputFH;

my @speicesPresentFinalStoreArray;

foreach my $lines (@speicesPresentArray) {
if ($lines =~ m/.fasta/) {
$lines =~ s/\|.*//;
$lines =~ s/a$/a /g;
push @speicesPresentFinalStoreArray, $lines;
} else {
#do nothing
}
}

sub uniqSpecies {
    my %speciesSeen;
    grep !$speciesSeen{$_}++, @_;
}

my @speicesPresentFinalArray = uniqSpecies(@speicesPresentFinalStoreArray);
print $speciesPresentFH @speicesPresentFinalArray;
