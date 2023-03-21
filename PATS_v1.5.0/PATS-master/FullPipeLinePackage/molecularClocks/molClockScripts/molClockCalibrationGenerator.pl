#!/usr/bin/perl

use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $mcBaseDir =  $ENV{'MolecularClockBaseDir'};

#Create directory to store generate calibration MEGA settings files
my $directoryMEGASettingsFile = "$mcBaseDir/userSettings/scenarioSettingsMEGA";
unless(-d $directoryMEGASettingsFile)
{
mkdir $directoryMEGASettingsFile or die "Unable to create $directoryMEGASettingsFile\n";
}

#Create create primary calibration folder
my $directoryMEGAPrimaryFile = "$directoryMEGASettingsFile/primaryCalibrations";
unless(-d $directoryMEGAPrimaryFile)
{
mkdir $directoryMEGAPrimaryFile or die "Unable to create $directoryMEGAPrimaryFile\n";
}


#location of primary calibration user settings/input
my @control;

my $controlfile = "$mcBaseDir/userSettings/userPrimaryCalibrationInput.ctl";
open (FH, '<', $controlfile) or die ("Cannot open $controlfile");
for (<FH>) {
next if /^#/;
push @control,$_;
}

for (@control){
	s/\r//g;
}

my @splitarray;
my @control2;

#Count array length
my $arraySizeControl = scalar(@control);

my $line;
my $i;

#iterate through ctl file and ignore line 0, only reads values after '= '
for ($i=0; $i < $arraySizeControl; $i++) {
$line = $control[$i];
chomp($line);
@control2 = split('= ', $line);
push(@splitarray,"$control2[1]\n");
undef(@control2);
}

my $arraySizeSplitArray = scalar (@splitarray);
my $numPrimCali = $arraySizeSplitArray / 3;

#Start generation of calibration file for MEGA
my $mainCaliTempString;

my @megaFile;
my $speciesA;
my $speciesB;
my $maxTime;
my $minTime;
my $speciesAItter = 0;
my $speciesBItter = 1;
my $trueTimeItter = 2;

############
#0 Balanced#
############
for (1 .. $numPrimCali) {
$mainCaliTempString = "!MRCA='taxonASpecies-taxonBSpecies' TaxonA='taxonASpecies' TaxonB='taxonBSpecies' MinTime=tiMin MaxTime=tiMax calibrationName='taxonASpecies-taxonBSpecies-split';";
push (@megaFile,"$mainCaliTempString\n");

$speciesA = $splitarray[$speciesAItter];
$speciesA =~ s/\s+//;
$speciesB = $splitarray[$speciesBItter];
$speciesB =~ s/\s+//;
$minTime = ($splitarray[$trueTimeItter]) - 1;
$minTime =~ s/\s+//;
$maxTime = ($splitarray[$trueTimeItter]) + 1;
$maxTime =~ s/\s+//;

for (@megaFile) {
s/\btaxonASpecies\b/$speciesA/g;
s/\btaxonBSpecies\b/$speciesB/g;
s/\btiMin\b/$minTime/g;
s/\btiMax\b/$maxTime/g;
}

$speciesAItter = $speciesAItter + 3;
$speciesBItter = $speciesBItter + 3;
$trueTimeItter = $trueTimeItter + 3;

undef $mainCaliTempString;
undef $speciesA;
undef $speciesB;
undef $minTime;
undef $maxTime;
}
for (@megaFile) {
	s/\n\r//g;
}

open my $zeroBalancedSettings, '>', "$directoryMEGAPrimaryFile/0BSettings.txt" or die "Cannot open 0BSettings.txt: $!";
print $zeroBalancedSettings @megaFile;
undef @megaFile;


#############
#10 Balanced#
#############

$speciesAItter = 0;
$speciesBItter = 1;
$trueTimeItter = 2;

for (1 .. $numPrimCali) {
$mainCaliTempString = "!MRCA='taxonASpecies-taxonBSpecies' TaxonA='taxonASpecies' TaxonB='taxonBSpecies' MinTime=tiMin MaxTime=tiMax calibrationName='taxonASpecies-taxonBSpecies-split';";
push (@megaFile,"$mainCaliTempString\n");

$speciesA = $splitarray[$speciesAItter];
$speciesA =~ s/\s+//;
$speciesB = $splitarray[$speciesBItter];
$speciesB =~ s/\s+//;
$minTime = ($splitarray[$trueTimeItter]) * 0.95;
$minTime =~ s/\s+//;
$maxTime = ($splitarray[$trueTimeItter]) * 1.05;
$maxTime =~ s/\s+//;

for (@megaFile) {
s/\btaxonASpecies\b/$speciesA/g;
s/\btaxonBSpecies\b/$speciesB/g;
s/\btiMin\b/$minTime/g;
s/\btiMax\b/$maxTime/g;
}

$speciesAItter = $speciesAItter + 3;
$speciesBItter = $speciesBItter + 3;
$trueTimeItter = $trueTimeItter + 3;

undef $mainCaliTempString;
undef $speciesA;
undef $speciesB;
undef $minTime;
undef $maxTime;
}
for (@megaFile) {
	s/\n\r//g;
}

open my $tenBalancedSettings, '>', "$directoryMEGAPrimaryFile/10BSettings.txt" or die "Cannot open 10BSettings.txt: $!";
print $tenBalancedSettings @megaFile;
undef @megaFile;

#############
#20 Balanced#
#############

$speciesAItter = 0;
$speciesBItter = 1;
$trueTimeItter = 2;

for (1 .. $numPrimCali) {
$mainCaliTempString = "!MRCA='taxonASpecies-taxonBSpecies' TaxonA='taxonASpecies' TaxonB='taxonBSpecies' MinTime=tiMin MaxTime=tiMax calibrationName='taxonASpecies-taxonBSpecies-split';";
push (@megaFile,"$mainCaliTempString\n");

$speciesA = $splitarray[$speciesAItter];
$speciesA =~ s/\s+//;
$speciesB = $splitarray[$speciesBItter];
$speciesB =~ s/\s+//;
$minTime = ($splitarray[$trueTimeItter]) * 0.9;
$minTime =~ s/\s+//;
$maxTime = ($splitarray[$trueTimeItter]) * 1.1;
$maxTime =~ s/\s+//;

for (@megaFile) {
s/\btaxonASpecies\b/$speciesA/g;
s/\btaxonBSpecies\b/$speciesB/g;
s/\btiMin\b/$minTime/g;
s/\btiMax\b/$maxTime/g;
}

$speciesAItter = $speciesAItter + 3;
$speciesBItter = $speciesBItter + 3;
$trueTimeItter = $trueTimeItter + 3;

undef $mainCaliTempString;
undef $speciesA;
undef $speciesB;
undef $minTime;
undef $maxTime;
}
for (@megaFile) {
	s/\n\r//g;
}

open my $twentyBalancedSettings, '>', "$directoryMEGAPrimaryFile/20BSettings.txt" or die "Cannot open 20BSettings.txt: $!";
print $twentyBalancedSettings @megaFile;
undef @megaFile;

########
#10 Low#
########

$speciesAItter = 0;
$speciesBItter = 1;
$trueTimeItter = 2;

for (1 .. $numPrimCali) {
$mainCaliTempString = "!MRCA='taxonASpecies-taxonBSpecies' TaxonA='taxonASpecies' TaxonB='taxonBSpecies' MinTime=tiMin MaxTime=tiMax calibrationName='taxonASpecies-taxonBSpecies-split';";
push (@megaFile,"$mainCaliTempString\n");

$speciesA = $splitarray[$speciesAItter];
$speciesA =~ s/\s+//;
$speciesB = $splitarray[$speciesBItter];
$speciesB =~ s/\s+//;
$minTime = ($splitarray[$trueTimeItter]) * 0.9;
$minTime =~ s/\s+//;
$maxTime = ($splitarray[$trueTimeItter]) + 1;
$maxTime =~ s/\s+//;

for (@megaFile) {
s/\btaxonASpecies\b/$speciesA/g;
s/\btaxonBSpecies\b/$speciesB/g;
s/\btiMin\b/$minTime/g;
s/\btiMax\b/$maxTime/g;
}

$speciesAItter = $speciesAItter + 3;
$speciesBItter = $speciesBItter + 3;
$trueTimeItter = $trueTimeItter + 3;

undef $mainCaliTempString;
undef $speciesA;
undef $speciesB;
undef $minTime;
undef $maxTime;
}
for (@megaFile) {
	s/\n\r//g;
}

open my $tenLowSettings, '>', "$directoryMEGAPrimaryFile/10LSettings.txt" or die "Cannot open 10LSettings.txt: $!";
print $tenLowSettings @megaFile;
undef @megaFile;

########
#20 Low#
########

$speciesAItter = 0;
$speciesBItter = 1;
$trueTimeItter = 2;

for (1 .. $numPrimCali) {
$mainCaliTempString = "!MRCA='taxonASpecies-taxonBSpecies' TaxonA='taxonASpecies' TaxonB='taxonBSpecies' MinTime=tiMin MaxTime=tiMax calibrationName='taxonASpecies-taxonBSpecies-split';";
push (@megaFile,"$mainCaliTempString\n");

$speciesA = $splitarray[$speciesAItter];
$speciesA =~ s/\s+//;
$speciesB = $splitarray[$speciesBItter];
$speciesB =~ s/\s+//;
$minTime = ($splitarray[$trueTimeItter]) * 0.8;
$minTime =~ s/\s+//;
$maxTime = ($splitarray[$trueTimeItter]) + 1;
$maxTime =~ s/\s+//;

for (@megaFile) {
s/\btaxonASpecies\b/$speciesA/g;
s/\btaxonBSpecies\b/$speciesB/g;
s/\btiMin\b/$minTime/g;
s/\btiMax\b/$maxTime/g;
}

$speciesAItter = $speciesAItter + 3;
$speciesBItter = $speciesBItter + 3;
$trueTimeItter = $trueTimeItter + 3;

undef $mainCaliTempString;
undef $speciesA;
undef $speciesB;
undef $minTime;
undef $maxTime;
}
for (@megaFile) {
	s/\n\r//g;
}

open my $twentyLowSettings, '>', "$directoryMEGAPrimaryFile/20LSettings.txt" or die "Cannot open 20LSettings.txt: $!";
print $twentyLowSettings @megaFile;
undef @megaFile;

#########
#10 High#
#########

$speciesAItter = 0;
$speciesBItter = 1;
$trueTimeItter = 2;

for (1 .. $numPrimCali) {
$mainCaliTempString = "!MRCA='taxonASpecies-taxonBSpecies' TaxonA='taxonASpecies' TaxonB='taxonBSpecies' MinTime=tiMin MaxTime=tiMax calibrationName='taxonASpecies-taxonBSpecies-split';";
push (@megaFile,"$mainCaliTempString\n");

$speciesA = $splitarray[$speciesAItter];
$speciesA =~ s/\s+//;
$speciesB = $splitarray[$speciesBItter];
$speciesB =~ s/\s+//;
$minTime = ($splitarray[$trueTimeItter]) - 1;
$minTime =~ s/\s+//;
$maxTime = ($splitarray[$trueTimeItter]) * 1.1;
$maxTime =~ s/\s+//;

for (@megaFile) {
s/\btaxonASpecies\b/$speciesA/g;
s/\btaxonBSpecies\b/$speciesB/g;
s/\btiMin\b/$minTime/g;
s/\btiMax\b/$maxTime/g;
}

$speciesAItter = $speciesAItter + 3;
$speciesBItter = $speciesBItter + 3;
$trueTimeItter = $trueTimeItter + 3;

undef $mainCaliTempString;
undef $speciesA;
undef $speciesB;
undef $minTime;
undef $maxTime;
}
for (@megaFile) {
	s/\n\r//g;
}

open my $tenHighSettings, '>', "$directoryMEGAPrimaryFile/10HSettings.txt" or die "Cannot open 10HSettings.txt: $!";
print $tenHighSettings @megaFile;
undef @megaFile;

#########
#20 High#
#########

$speciesAItter = 0;
$speciesBItter = 1;
$trueTimeItter = 2;

for (1 .. $numPrimCali) {
$mainCaliTempString = "!MRCA='taxonASpecies-taxonBSpecies' TaxonA='taxonASpecies' TaxonB='taxonBSpecies' MinTime=tiMin MaxTime=tiMax calibrationName='taxonASpecies-taxonBSpecies-split';";
push (@megaFile,"$mainCaliTempString\n");

$speciesA = $splitarray[$speciesAItter];
$speciesA =~ s/\s+//;
$speciesB = $splitarray[$speciesBItter];
$speciesB =~ s/\s+//;
$minTime = ($splitarray[$trueTimeItter]) - 1;
$minTime =~ s/\s+//;
$maxTime = ($splitarray[$trueTimeItter]) * 1.2;
$maxTime =~ s/\s+//;

for (@megaFile) {
s/\btaxonASpecies\b/$speciesA/g;
s/\btaxonBSpecies\b/$speciesB/g;
s/\btiMin\b/$minTime/g;
s/\btiMax\b/$maxTime/g;
}

$speciesAItter = $speciesAItter + 3;
$speciesBItter = $speciesBItter + 3;
$trueTimeItter = $trueTimeItter + 3;

undef $mainCaliTempString;
undef $speciesA;
undef $speciesB;
undef $minTime;
undef $maxTime;
}
for (@megaFile) {
	s/\n\r//g;
}

open my $twentyHighSettings, '>', "$directoryMEGAPrimaryFile/20HSettings.txt" or die "Cannot open 20HSettings.txt: $!";
print $twentyHighSettings @megaFile;
undef @megaFile;
