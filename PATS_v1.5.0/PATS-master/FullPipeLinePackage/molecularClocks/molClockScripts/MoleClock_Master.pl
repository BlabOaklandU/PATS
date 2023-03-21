#!/bin/perl.exe
use strict;
use warnings;
use Cwd qw(getcwd);
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;


my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $mcBaseDir =  $ENV{'MolecularClockBaseDir'};
my $mcScriptsDir = $ENV{'MolecularClockScriptDir'};


#generate log
open my $logMolClockDebug, '>', "$mcBaseDir/molClockDebugLog.txt" or die "Cannot open molClockDebugLog.txt: $!";
print $logMolClockDebug "Starting the molecular clock analysis\n";


#Generation of floders will have to be added to initial package

#Generate all Directories for the full MoleClock Analysis
my $mcDirectory1 = "$mcBaseDir/unsplitMainTree";
unless(-d $mcDirectory1)
{
mkdir $mcDirectory1 or die "Unable to create $mcDirectory1\n";
print $logMolClockDebug "Missing directory: $mcDirectory1\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory2 = "$mcBaseDir/nullCheck";
unless(-d $mcDirectory2)
{
mkdir $mcDirectory2 or die "Unable to create $mcDirectory2\n";
print $logMolClockDebug "Missing directory: $mcDirectory2\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory3 = "$mcBaseDir/phyFiles";
unless(-d $mcDirectory3)
{
mkdir $mcDirectory3 or die "Unable to create $mcDirectory3\n";
print $logMolClockDebug "Missing directory: $mcDirectory3\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory4 = "$mcBaseDir/splitTrees";
unless(-d $mcDirectory4)
{
mkdir $mcDirectory4 or die "Unable to create $mcDirectory4\n";
print $logMolClockDebug "Missing directory: $mcDirectory4\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory5 = "$mcBaseDir/splitTrees/treeA";
unless(-d $mcDirectory5)
{
mkdir $mcDirectory5 or die "Unable to create $mcDirectory5\n";
print $logMolClockDebug "Missing directory: $mcDirectory5\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory6 = "$mcBaseDir/splitTrees/treeB";
unless(-d $mcDirectory6)
{
mkdir $mcDirectory6 or die "Unable to create $mcDirectory6\n";
print $logMolClockDebug "Missing directory: $mcDirectory6\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

#Calibration Tree
my $mcDirectory7 = "$mcBaseDir/calibrationTrees";
unless(-d $mcDirectory7)
{
mkdir $mcDirectory7 or die "Unable to create $mcDirectory7\n";
print $logMolClockDebug "Missing directory: $mcDirectory7\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory8 = "$mcBaseDir/calibrationTrees/treeA";
unless(-d $mcDirectory8)
{
mkdir $mcDirectory8 or die "Unable to create $mcDirectory8\n";
print $logMolClockDebug "Missing directory: $mcDirectory8\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory9 = "$mcBaseDir/calibrationTrees/treeB";
unless(-d $mcDirectory9)
{
mkdir $mcDirectory9 or die "Unable to create $mcDirectory9\n";
print $logMolClockDebug "Missing directory: $mcDirectory9\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory10 = "$mcBaseDir/orthoGeneSequences";
unless(-d $mcDirectory10)
{
mkdir $mcDirectory10 or die "Unable to create $mcDirectory10\n";
print $logMolClockDebug "Missing directory: $mcDirectory10\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory11 = "$mcBaseDir/splitGenes";
unless(-d $mcDirectory11)
{
mkdir $mcDirectory11 or die "Unable to create $mcDirectory11\n";
print $logMolClockDebug "Missing directory: $mcDirectory11\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory12 = "$mcBaseDir/groupConcatGenes";
unless(-d $mcDirectory12)
{
mkdir $mcDirectory12 or die "Unable to create $mcDirectory12\n";
print $logMolClockDebug "Missing directory: $mcDirectory12\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}

my $mcDirectory13 = "$mcBaseDir/megaOutput";
unless(-d $mcDirectory13)
{
mkdir $mcDirectory13 or die "Unable to create $mcDirectory13\n";
print $logMolClockDebug "Missing directory: $mcDirectory13\n";
print $logMolClockDebug "Analysis terminated, check directory for issues";
die;
}