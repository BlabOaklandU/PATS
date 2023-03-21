#!/usr/bin/perl
  
use strict;
use warnings;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};

#create directory
my $directory1 = "$baseDir/orthoGroupsFasta";

unless(-d $directory1)
{
mkdir $directory1 or die "Unable to create $directory1\n";
}

my $OGDir = "$baseDir/og";
opendir my $OGdirfh, $OGDir  or die "Can't open $OGDir: $!";

my @curSpecies;

my $fileItter = 0;
my $numOGFiles =  grep { -f "$OGDir/$_" } readdir($OGdirfh);

#Open Species list (ordered)
open my $speciesHandle, '<', "$baseDir/logFiles/SpeciesList.txt" or die "Cannot open SpeciesList.txt: $!";
while(<$speciesHandle>) { 
    chomp; 
    push @curSpecies,$_;
} 
close $speciesHandle;

#Push species to array (ordered)
my $speciesString = $curSpecies[0];
my @speciesArray = split ' ', $speciesString;

my $fastaDir = "$baseDir/fasta";
opendir my $dh, $fastaDir  or die "Can't open $fastaDir: $!";
my @numFastaFilesArray = grep(/\.fasta$/,readdir($dh));
my $numFastaFiles =  scalar(@numFastaFilesArray);

#Initiate loop for each OG file
foreach (0 .. $numOGFiles-1) {

my @orthoGroup;

#Read previously created og file
open my $ogfh, '<', "$baseDir/og/OG$fileItter.txt" or die "Cannot open OG$fileItter.txt: $!";

#Push previously created OG file to array (per file basis)
while(<$ogfh>) { 
    chomp; 
    push @orthoGroup, $_;
} 

close $ogfh;

#organize array by removing first 3 columns
my $orthoGroupString = $orthoGroup[0];
my @arrayOG = split ' ', $orthoGroupString;
splice @arrayOG, 0, 3;

#Get OG file length (ordered)
my $arrayLengthOG = scalar(@arrayOG);

my @fastaFileFinal;
my $count = 0;

foreach (1 .. $arrayLengthOG){

my @fastaArrayTemp;

#Open FASTA file of species
open my $fastaHandle, '<', "$baseDir/fasta/$speciesArray[$count]" or die "Cannot open $speciesArray[$count]: $!";

#Push fasta file information to array
while(<$fastaHandle>) { 
    chomp; 
    push @fastaArrayTemp, $_;
} 
close $fastaHandle;

#join fasta by white space, split by >, store each chunk in array

my $fastaString = join ('', @fastaArrayTemp);

my @fastaArray = grep { /\S/ } split(/[>]/, $fastaString);

#my $arrayLengthOFA = scalar(@fastaArray);
my $searchString = $arrayOG[$count];


my ($searchCrit) = grep( m/\Q$searchString\E/, @fastaArray );
push @fastaFileFinal, ">$searchCrit\n";

$count ++;
}

#Create OG fasta file for storing
open my $FastaFH, '>', "$baseDir/orthoGroupsFasta/OG$fileItter.fasta" or die "Cannot open OrthoGroups$fileItter.fasta: $!";

for (@fastaFileFinal) 
{
  s/\s+/\n/g;
  s/^>\n$//g;
}

foreach (@fastaFileFinal){
print $FastaFH $_;
}
close $FastaFH;

$fileItter ++;

undef @fastaFileFinal;
undef @orthoGroup
}
