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


#Directories
my $logFilesDir = "$baseDir/logFiles";
#my $logFilesTwoDir = "$baseDir/DataStandardizerLogs";
my $muscleFastaDir = "$baseDir/muscleFasta";
my $ogDir = "$baseDir/og";
my $OGGapFilterDir = "$baseDir/OGGapFilter";
my $orthoGroupsFastaDir = "$baseDir/orthoGroupsFasta";
my $orthologsDir = "$baseDir/orthologs";
my $poWorkDir = "$baseDir/poWork";
my $treeRAxMLDir = "$baseDir/treeRAxML";
my $fastTreeDir = "$baseDir/fastTree";
my $concatFilterGroupDir = "$baseDir/concatFilterGroup";
my $concatFilterIndvDir = "$baseDir/concatFilterIndv";
my $directoryNewGap = "$baseDir/cleanedGap";
my $stdDir = "$baseDir/std";

#Files
my $archiveLogFile = "$baseDir/archiveLog.txt";

#Delete Old Folders after Archive

if ($ENV{'PoPipe_createArchive'}  == 1 && $ENV{'PoPipe_folderCleanUp'} == 1) {

rmtree $logFilesDir;
rmtree $muscleFastaDir;
rmtree $ogDir;
rmtree $OGGapFilterDir;
rmtree $orthologsDir;
rmtree $orthoGroupsFastaDir;
rmtree $poWorkDir;
rmtree $treeRAxMLDir;
rmtree $fastTreeDir;
rmtree $concatFilterGroupDir;
rmtree $concatFilterIndvDir;
rmtree $directoryNewGap;
#rmtree $logFilesTwoDir;

unlink $archiveLogFile;

#Delete all files within the std folder

unlink glob "$baseDir/std/*.err";
unlink glob "$baseDir/std/*.out";

} else {

#Nothing Happens

}




#######################
#Keep One Permutations#
#######################

if ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1) {


my $filename = "$baseDir/permutations/keepOne.txt";

if (-z $filename) {

my $filenameOld = "$baseDir/permutations/keepOne.txt";
unlink $filenameOld;


my $filenameBak = "$baseDir/permutations/keepOneBU.txt";
my $filenameNew = "$baseDir/permutations/keepOne.txt";
rename($filenameBak,$filenameNew) or die ("Unable to rename $filenameBak to $filenameNew: $!\n");

my $directoryFastaTemp = "$baseDir/fasta/temp";
rmtree $directoryFastaTemp;


my $directorySTD = "$baseDir/std";
rmtree $directorySTD;

} else {

exec("bash", "poPipeFull.sh") or die ("Unable to run .sh script");

exit;

}


###########################
#Remove Group Permutations#
###########################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 1) {


my $directoryFastaTemp = "$baseDir/fasta/temp";
rmtree $directoryFastaTemp;
my $directorySTD = "$baseDir/std";
rmtree $directorySTD;



#########################
#Remove One Permutations#
#########################

} elsif ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1){


my $filename = "$baseDir/permutations/removeOne.txt";

if (-z $filename) {

my $filenameOld = "$baseDir/permutations/removeOne.txt";
unlink $filenameOld;


my $filenameBak = "$baseDir/permutations/removeOneBU.txt";
my $filenameNew = "$baseDir/permutations/removeOne.txt";
rename($filenameBak,$filenameNew) or die ("Unable to rename $filenameBak to $filenameNew: $!\n");

my $directoryFastaTemp = "$baseDir/fasta/temp";
rmtree $directoryFastaTemp;


my $directorySTD = "$baseDir/std";
rmtree $directorySTD;


}else {

exec("bash", "poPipeFull.sh") or die ("Unable to run .sh script");

exit;

}

#####################################
#Adjusting Cleaner for reruntreeOnly#
#####################################

} elsif ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_createArchive'}  == 1 && $ENV{'PoPipe_folderCleanUp'} == 0){

#Do Nothing

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_createArchive'}  == 1 && $ENV{'PoPipe_folderCleanUp'} == 0){

#Do Nothing

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1 && $ENV{'PoPipe_createArchive'}  == 1 && $ENV{'PoPipe_folderCleanUp'} == 0){

#Do Nothing

#################################
#Adjusting Cleaner for Non-Perms#
#################################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_createArchive'}  == 0 && $ENV{'PoPipe_folderCleanUp'} == 0){

#Do Nothing

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_createArchive'}  == 1 && $ENV{'PoPipe_folderCleanUp'} == 0){

#Do Nothing

} else {

die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");

}
