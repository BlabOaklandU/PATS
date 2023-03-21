#!/usr/bin/perl
  
use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;


if ($ENV{'PoPipe_createArchive'} == 1){

my $baseDir = $ENV{'PoPipe_srcBaseDir'};

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


#Generate Archive Log
open my $logArchiveFH, '>', "$baseDir/archiveLog.txt" or die "Cannot open archiveLog.txt: $!";


#Zip Directory Name

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


#Date/Time
my $dateString = localtime();

print $logArchiveFH "Date of archive: $dateString\n\n";


my @curSpecies1;

#Push all species in analysis to archive log.

print $logArchiveFH "Species in analysis:\n\n";

#Create ordered array with fasta file names
open my $speciesHandle1, '<', "$baseDir/logFiles/SpeciesList.txt" or die "Cannot open SpeciesList.txt: $!";
while(<$speciesHandle1>) { 
    chomp; 
    push @curSpecies1,$_;
} 
close $speciesHandle1;

foreach my $fastaSpeciesList (@curSpecies1) {
print $logArchiveFH "$fastaSpeciesList\n";
}

print $logArchiveFH "\n\n";


#Create a zip directory

my $directoryArchive = "$baseDir/$permName-$treeName-$runName";
unless(-d $directoryArchive)
{
mkdir $directoryArchive or die "Unable to create $directoryArchive\n";
}

#Create Arch

#Create sub directories with in zip directory

#Copy DataStandardizerLogs Directory
my $directoryDSL = "$directoryArchive/DataStandardizerLogs";
unless(-d $directoryDSL)
{
mkdir $directoryDSL or die "Unable to create $directoryDSL\n";
}
my $copyDSL = "$baseDir/DataStandardizerLogs";
dircopy($copyDSL,$directoryDSL) or die $!;

#Copy cleanedGaps Dirctory
my $directoryNewGap = "$directoryArchive/cleanedGap";
unless(-d $directoryNewGap)
{
mkdir $directoryNewGap or die "Unable to create $directoryNewGap\n";
}
my $copyNewGap = "$baseDir/cleanedGap";
dircopy($copyNewGap,$directoryNewGap) or die $!;

#Copy og Directory
my $directoryOG = "$directoryArchive/og";
unless(-d $directoryOG)
{
mkdir $directoryOG or die "Unable to create $directoryOG\n";
}
my $copyOG = "$baseDir/og";
dircopy($copyOG,$directoryOG) or die $!;


#Copy concatFilterGroup Directory
my $directoryconcatFilterGroup = "$directoryArchive/concatFilterGroup";
unless(-d $directoryconcatFilterGroup)
{
mkdir $directoryconcatFilterGroup or die "Unable to create $directoryconcatFilterGroup\n";
}
my $copyconcatFilterGroup = "$baseDir/concatFilterGroup";

dircopy($copyconcatFilterGroup,$directoryconcatFilterGroup) or die $!;


#Copy concatFilterIndv Directory
my $directoryconcatFilterIndv = "$directoryArchive/concatFilterIndv";
unless(-d $directoryconcatFilterIndv)
{
mkdir $directoryconcatFilterIndv or die "Unable to create $directoryconcatFilterIndv\n";
}
my $copyconcatFilterIndv = "$baseDir/concatFilterIndv";
dircopy($copyconcatFilterIndv,$directoryconcatFilterIndv) or die $!;


#Copy fasta Directory
my $directoryfasta = "$directoryArchive/fasta";
unless(-d $directoryfasta)
{
mkdir $directoryfasta or die "Unable to create $directoryfasta\n";
}
my $copyfasta = "$baseDir/fasta";
dircopy($copyfasta,$directoryfasta) or die $!;


#Copy muscleFasta Directory
my $directorymuscleFasta = "$directoryArchive/muscleFasta";
unless(-d $directorymuscleFasta)
{
mkdir $directorymuscleFasta or die "Unable to create $directorymuscleFasta\n";
}
my $copymuscleFasta = "$baseDir/muscleFasta";
dircopy($copymuscleFasta,$directorymuscleFasta) or die $!;


#Copy OGGapFilter Directory
my $directoryOGGapFilter = "$directoryArchive/OGGapFilter";
unless(-d $directoryOGGapFilter)
{
mkdir $directoryOGGapFilter or die "Unable to create $directoryOGGapFilter\n";
}
my $copyOGGapFilter = "$baseDir/OGGapFilter";
dircopy($copyOGGapFilter,$directoryOGGapFilter) or die $!;


#Copy orthoGroupsFasta Directory
my $directoryorthoGroupsFasta = "$directoryArchive/orthoGroupsFasta";
unless(-d $directoryorthoGroupsFasta)
{
mkdir $directoryorthoGroupsFasta or die "Unable to create $directoryorthoGroupsFasta\n";
}
my $copyorthoGroupsFasta = "$baseDir/orthoGroupsFasta";
dircopy($copyorthoGroupsFasta,$directoryorthoGroupsFasta) or die $!;


#Copy orthologs Directory
my $directoryorthologs = "$directoryArchive/orthologs";
unless(-d $directoryorthologs)
{
mkdir $directoryorthologs or die "Unable to create $directoryorthologs\n";
}
my $copyorthologs = "$baseDir/orthologs";
dircopy($copyorthologs,$directoryorthologs) or die $!;


#Copy permutations Directory
my $directorypermutations = "$directoryArchive/permutations";
unless(-d $directorypermutations)
{
mkdir $directorypermutations or die "Unable to create $directorypermutations\n";
}
my $copypermutations = "$baseDir/permutations";
dircopy($copypermutations,$directorypermutations) or die $!;


#Copy poWork Directory
my $directorypoWork = "$directoryArchive/poWork";
unless(-d $directorypoWork)
{
mkdir $directorypoWork or die "Unable to create $directorypoWork\n";
}
my $copypoWork = "$baseDir/poWork";
dircopy($copypoWork,$directorypoWork) or die $!;


#Copy std Directory
my $directorystd = "$directoryArchive/std";
unless(-d $directorystd)
{
mkdir $directorystd or die "Unable to create $directorystd\n";
}
my $copystd = "$baseDir/std";
dircopy($copystd,$directorystd) or die $!;


#Copy treeRAxML Directory
my $directorytreeRAxML = "$directoryArchive/treeRAxML";
unless(-d $directorytreeRAxML)
{
mkdir $directorytreeRAxML or die "Unable to create $directorytreeRAxML\n";
}
my $copytreeRAxML = "$baseDir/treeRAxML";
dircopy($copytreeRAxML,$directorytreeRAxML) or die $!;


#Copy fastTree Directory
my $directoryfastTree = "$directoryArchive/fastTree";
unless(-d $directoryfastTree)
{
mkdir $directoryfastTree or die "Unable to create $directoryfastTree\n";
}
my $copyfastTree = "$baseDir/fastTree";
dircopy($copyfastTree,$directoryfastTree) or die $!;


#Copy log.txt Directory
my $filelogtxt = "$directoryArchive/";
unless(-d $filelogtxt)
{
mkdir $filelogtxt or die "Unable to create $filelogtxt\n";
}
my $copyfilelogtxt = "$baseDir/log.txt";
fcopy($copyfilelogtxt,$filelogtxt) or die $!;


#Copy logFiles Directory
my $directorylogFiles = "$directoryArchive/logFiles";
unless(-d $directorylogFiles)
{
mkdir $directorylogFiles or die "Unable to create $directorylogFiles\n";
}
my $copylogFiles = "$baseDir/logFiles";
dircopy($copylogFiles,$directorylogFiles) or die $!;


my $archName = "$permName-$treeName-$runName";
my $archNameOld = "$permName-$treeName-$runName.7z";
my $archNameNew = "$permName-$treeName-$runName-$dateString.7z";

print $logArchiveFH "Name of archive: $archNameNew\n";

my $fromDir = "$baseDir/archiveLog.txt";
my $toDir = "$baseDir/$archName";

close $logArchiveFH;

fmove ($fromDir, $toDir) or die "error";

#Initialize zip if files exist.

system("7za a $permName-$treeName-$runName.7z $archName");

rename $archNameOld,$archNameNew;

} else {

#Nothing happens
}