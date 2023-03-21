#!/usr/bin/perl

use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $mcBaseDir =  $ENV{'MolecularClockBaseDir'};
my $ogindvGeneSequences = $ENV{'indvGenSequencesDir'};
my $scenarioMEGADir = "$mcBaseDir/userSettings/scenarioSettingsMEGA";
my $megaLoc = "$ENV{'PoPipe_srcBaseDir'}/externalSoftware/MEGA/11.0.13/RH/usr/lib64/mega/megacc";
my $megaMaoLoc = '$mcBaseDir/userSettings/megaSettings.mao';
my $megaOutgroup = $ENV{'outgroupSpecies'};
my $megaOutputDir = $ENV{'megaOutputDir'};

#Primary calibration locaitons
my $primeCaliLoc0B = "$scenarioMEGADir/primaryCalibrations/0BSettings.txt";
my $primeCaliLoc10B = "$scenarioMEGADir/primaryCalibrations/10BSettings.txt";
my $primeCaliLoc20B = "$scenarioMEGADir/primaryCalibrations/20BSettings.txt";
my $primeCaliLoc10L = "$scenarioMEGADir/primaryCalibrations/10LSettings.txt";
my $primeCaliLoc20L = "$scenarioMEGADir/primaryCalibrations/20LSettings.txt";
my $primeCaliLoc10H = "$scenarioMEGADir/primaryCalibrations/10HSettings.txt";
my $primeCaliLoc20H = "$scenarioMEGADir/primaryCalibrations/20HSettings.txt";

open my $logMolClockDebug, '>>', "$mcBaseDir/molClockDebugLog.txt" or die "Cannot open molClockDebugLog.txt: $!";

my $timestamp0 = localtime(time);

print $logMolClockDebug "\nRunning molClockConcatGene.pl... $timestamp0 \n\n";

if ($ENV{'moleClocksOnly'} == 1 && $ENV{'moleClocksBasicPats'} == 0 && $ENV{'moleClockPermPats'} == 0) {
#Running just Molecular Analysis, no PATS
print $logMolClockDebug "Checking for needed files...\n\n";

my @concatGeneFile = glob ( "$ENV{'concatGenSequencesDir'}/*.meg" );
my $concatArrayCount = scalar @concatGeneFile;
my $concatString = $concatGeneFile[0];
my $concatGeneFileName = "$ENV{'concatGenSequencesDir'}/$concatString";

if ($concatArrayCount == 1) {
print $logMolClockDebug "One .fasta file found.\n";
print $logMolClockDebug "$concatGeneFileName\n";
} else {
print $logMolClockDebug "There is more than one .fasta file present\n";
print $logMolClockDebug "Please check '$ENV{'concatGenSequencesDir'}'\n";
exit;
}

my @nwkTreeFiles = glob ( "$ENV{'unsplitMainTreeDir'}/*.nwk" );
my $nwkArrayCount = scalar @nwkTreeFiles;
my $nwkString = $nwkTreeFiles[0];
my $nwkTreeFileName = "$ENV{'unsplitMainTreeDir'}/$nwkString";

if ($nwkArrayCount == 1) {
print $logMolClockDebug "One .nwk tree file found.\n";
} else {
print $logMolClockDebug "There is more than one .nwk tree file present\n";
print $logMolClockDebug "Please check '$ENV{'unsplitMainTreeDir'}'\n";
exit;
}

undef @nwkTreeFiles;
undef @concatGeneFile;

my $timestamp1 = localtime(time);
print $logMolClockDebug "Starting MEGA v11.0.13 at $timestamp1";
system("${megaLoc} -a ${megaMaoLoc} -d ${concatGeneFileName} -t ${nwkTreeFileName} -g ${megaOutgroup} -c ${primeCaliLoc0B} -o ${megaOutputDir}");
my $timestamp2 = localtime(time);
print $logMolClockDebug "MEGA 11.0.13 has finished at $timestamp2";


}elsif ($ENV{'moleClocksOnly'} == 0 && $ENV{'moleClocksBasicPats'} == 1 && $ENV{'moleClockPermPats'} == 0){
#Running PATS (standard run) + Molecular Analysis








}elsif ($ENV{'moleClocksOnly'} == 0 && $ENV{'moleClocksBasicPats'} == 0 && $ENV{'moleClockPermPats'} == 1){
#Running PATS with permutations then running Molecular clock Analyses	
#This will be added later.
	
}else{
print $logMolClockDebug "Molecular clock settings not set.\n";
print $logMolClockDebug "Please check poPipeSettings.sh\n";
exit;
}

	
	
	
	
	
	
	
	
	
	
	
	