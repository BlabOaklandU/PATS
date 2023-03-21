#!/usr/bin/perl
  
use strict;
use warnings;
use File::Path;
use File::Find;
use File::Copy;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);


my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $archBin = $ENV{'PoPipe_createArchive'};

#Tree selection
my $fastTree = $ENV{'PoPipe_fastTree'};
my $raxmlTree = $ENV{'PoPipe_RAxMLTree'};

#Tree versus full run
my $rerunTreeOnly = $ENV{'PoPipe_rerunTreesOnly'};
my $rerunFullPipeline = $ENV{'PoPipe_rerunFullAnalysis'};

#Permutations
my $keepOnePerm = $ENV{'PoPipe_keepOne'};
my $removeOnePerm = $ENV{'PoPipe_removeOne'};
my $removeGroupPerm = $ENV{'poPipe_removeGroup'};


my $permName;
my $treeName;
my $runName;

#Permutation settings
if ($keepOnePerm == 1 && $removeOnePerm == 0 && $removeGroupPerm == 0) {
$permName = 'KeepOne';
}elsif ($keepOnePerm == 0 && $removeOnePerm == 1 && $removeGroupPerm == 0) {
$permName = 'RemoveOne';
}elsif ($keepOnePerm == 0 && $removeOnePerm == 0 && $removeGroupPerm == 1) {
$permName = 'RemoveGroup';
}elsif ($keepOnePerm == 0 && $removeOnePerm == 0 && $removeGroupPerm == 0) {
$permName = 'Standard';
}else {
die ("Error");
}

#Tree settings
if ($fastTree == 1 && $raxmlTree == 0) {
$treeName = 'FastTree';
}elsif ($fastTree == 0 && $raxmlTree == 1) {
$treeName = 'RAxMLTree';
}else{
die ("Error");
}

#Tree settings
if ($rerunTreeOnly == 1 && $rerunFullPipeline == 0) {
$runName = 'RerunTreesOnly';
}elsif ($rerunTreeOnly == 0 && $rerunFullPipeline == 1) {
$runName = 'RerunFullAnalysis';
}elsif ($rerunTreeOnly == 0 && $rerunFullPipeline == 0) {
$runName = 'StandardFullAnalysis';
}else{
die ("Error");
}

my $directoryArchive = "$baseDir/$permName-$treeName-$runName";
rmtree $directoryArchive;
