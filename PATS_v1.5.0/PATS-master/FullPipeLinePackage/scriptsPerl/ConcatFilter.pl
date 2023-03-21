#!/usr/bin/perl

use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};

my $fastTreeSetting = $ENV{'PoPipe_fastTree'};
my $raxmlTreeSetting = $ENV{'PoPipe_RAxMLTree'};


#generate log
open my $logconcatFH, '>', "$baseDir/logFiles/concatFilterLog.txt" or die "Cannot open concatFilterLog.txt: $!";

#Selected permutations
print $logconcatFH "Permutations Selected:\n";
print $logconcatFH "Keep One: $ENV{'PoPipe_keepOne'}\n";
print $logconcatFH "Remove One: $ENV{'PoPipe_removeOne'}\n";
print $logconcatFH "Remove Group: $ENV{'poPipe_removeGroup'}\n";

print $logconcatFH "Rerunning just trees: $ENV{'PoPipe_rerunTreesOnly'}\n";


#Create directory
my $directory2 = "$baseDir/concatFilterIndv";
unless(-d $directory2)
{
mkdir $directory2 or die "Unable to create $directory2\n";
}

my $directorygapFil = "$baseDir/concatFilterGroup";
unless(-d $directorygapFil)
{
mkdir $directorygapFil or die "Unable to create $directorygapFil\n";
}

my $muscleDir = "$baseDir/muscleFasta";
opendir my $muscleDirfh, $muscleDir  or die "Can't open $muscleDir: $!";

my $numMuscleFiles =  grep { -f "$muscleDir/$_" } readdir($muscleDirfh);
my @curSpecies1;

#Create ordered array with fasta file names
open my $speciesHandle1, '<', "$baseDir/logFiles/SpeciesList.txt" or die "Cannot open SpeciesList.txt: $!";
while(<$speciesHandle1>) { 
    chomp;
    push @curSpecies1,$_;
} 
close $speciesHandle1;

my $speciesString1 = $curSpecies1[0];
my @speciesArray = split ' ', $speciesString1;

my $numSpecies1 = scalar(grep {defined $_} @speciesArray),;
my $itter = 0;

#Create individual species files
foreach (1 .. $numSpecies1) {

open my $FH, '>', "$baseDir/concatFilterIndv/$speciesArray[$itter]" or die "Cannot open concatFilter $speciesArray[$itter]: $!";
print $FH ">$speciesArray[$itter]\n";
close $FH;
$itter++;
}

my @muscleArray;
my @muscleArrayTemp;
my @muscleTemp;
my $fileItter = 0;
my $counter = 0;

#Initiate loop for each OG file
foreach (1 .. $numMuscleFiles) {

#Read previously created og file
open my $musclefh, '<', "$baseDir/cleanedGap/OG$fileItter.fasta" or die "Cannot open OG$fileItter.txt: $!";

#Push Muscle file to array (per file basis)
while(<$musclefh>) {  
push @muscleTemp, $_;
} 
close $musclefh;

for (@muscleTemp) 
{
  s/\|.*//;
}

my $muscleString = join ('', @muscleTemp);


my @muscleArrayTemp = grep { /\S/ } split(/[>]/, $muscleString);


#Initiate loop to search for SP name and store in proper file
foreach (1 .. $numSpecies1){

open my $speciesFH, '>>', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

my $searchString = $speciesArray[$counter];

if(grep(m/^$searchString/, @muscleArrayTemp)) {
my ($searchCrit) = grep( m/^$searchString/, @muscleArrayTemp);
push @muscleArray, "$searchCrit";
}else {
#Do Nothing
}

for (@muscleArray) {
   s/.*?\s//;
   s/\s+//g;
}

print $speciesFH @muscleArray;
close $speciesFH;

$counter++;

undef @muscleArray;
}
$fileItter++;
$counter = 0;

undef @muscleTemp;
undef @muscleArrayTemp;
}


#Remove species that do not have any genes

my $itterCheck = 0;
my $numLines;

foreach (1 .. $numSpecies1) {
open my $FH2, '<', "$baseDir/concatFilterIndv/$speciesArray[$itterCheck]" or die "Cannot open concatFilter $speciesArray[$itterCheck]: $!";

$numLines++ while (<$FH2>);
close $FH2;

if ($numLines > 1) {
#do nothing
} else {
unlink ("$baseDir/concatFilterIndv/$speciesArray[$itterCheck]");
}
$itterCheck++;
undef $numLines;
}

my $CFIndvDir = "$baseDir/concatFilterIndv/";
opendir my $CFIndvdirfh, $CFIndvDir  or die "Can't open $CFIndvDir: $!";
my $numCFIndvFiles =  grep { -f "$CFIndvDir/$_" } readdir($CFIndvdirfh);



#########################
#RAXML #######################################################################################################################################################
#########################

if($ENV{'PoPipe_RAxMLTree'} == 1 && $ENV{'PoPipe_fastTree'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 0){

#########################
#Remove One Permutations#
#########################

if ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {

my $directoryRemoveOne = "$baseDir/treeRAxML/removeOne";
unless(-d $directoryRemoveOne)
{
mkdir $directoryRemoveOne or die "Unable to create $directoryRemoveOne\n";
}

#generate log
open my $logRemoveOneFH, '>', "$baseDir/logFiles/removeOneLog.txt" or die "Cannot open removeOneLog.txt: $!";

my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp;
    s/\r//g;	
    push @removeOneList,$_;
} 
close $removeOneFH;


my $numRemoveOne =  scalar @removeOneList;
my $index = 0;

print $logRemoveOneFH "Number of Species to be removed: $numRemoveOne\n\n";


foreach (1 .. $numRemoveOne){

print $logRemoveOneFH "Currently removing $removeOneList[$index] from analysis\n\n";

my $fileRename = "rename.txt";

rename ("$baseDir/concatFilterIndv/$removeOneList[$index]", "$baseDir/concatFilterIndv/$fileRename");


my @finalArray;
my $counterTwo = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]" or die "Cannot open concatFilter $speciesArray[$counterTwo]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counterTwo++;

}else{
#Do Nothing
$counterTwo++;
}
}


for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveOne_GroupFasta_$removeOneList[$index]" or die "Cannot open RemoveOne_GroupFasta_$removeOneList[$index]: $!";
print $singleFilePush @finalArray;
close $singleFilePush;

#Make directoy for each RemoveOne scenario

my $dirRemoveOne = "$baseDir/treeRAxML/removeOne/RemoveOne_$removeOneList[$index]";
unless(-d $dirRemoveOne)
{
mkdir $dirRemoveOne or die "Unable to create $dirRemoveOne\n";
}

rename ("$baseDir/concatFilterIndv/$fileRename", "$baseDir/concatFilterIndv/$removeOneList[$index]");

print $logRemoveOneFH "Group.fasta without $removeOneList[$index] included will be called:\n"; 
print $logRemoveOneFH "RemoveOne_GroupFasta_$removeOneList[$index]\n\n";


undef @finalArray;
undef $counterTwo;

$index ++;

}
 close $logRemoveOneFH;
 

#######################
#Keep One Permutations#
#######################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {


my $directoryKeepOne = "$baseDir/treeRAxML/keepOne";
unless(-d $directoryKeepOne)
{
mkdir $directoryKeepOne or die "Unable to create $directoryKeepOne\n";
}

#Create temp directory to fold moved files

my $directoryKeepOneTemp = "$baseDir/concatFilterIndv/temp";
unless(-d $directoryKeepOneTemp)
{
mkdir $directoryKeepOneTemp or die "Unable to create $directoryKeepOneTemp\n";
}


#Generate KeepOne Log

open my $logKeepOneFH, '>', "$baseDir/logFiles/keepOneLog.txt" or die "Cannot open keepOneLog.txt: $!";

my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp;
    s/\r//g;	
    push @keepOneList,$_;
} 
close $keepOneFH;

my $numKeepOne =  scalar @keepOneList;



#Move All Keep One Files to temp folder.

print $logKeepOneFH "Files being moved for initial setup:\n";

my $index = 0;
foreach (1 .. $numKeepOne) {

my $fromDir = "$baseDir/concatFilterIndv/$keepOneList[$index]";
my $toDir = "$baseDir/concatFilterIndv/temp/";

fmove ($fromDir, $toDir) or die print $logKeepOneFH "All keepOne species were not present in directory!\n";

print $logKeepOneFH "$fromDir was just moved to $toDir\n";

$index ++;
}


print $logKeepOneFH "\n\n";


#Start multiple group.fasta creations
my $indexTwo = 0;

foreach (1 .. $numKeepOne) {

#Move the Keep One back to run with analysis
print $logKeepOneFH "keepOne file being moved back into analysis:\n";

my $fromDir = "$baseDir/concatFilterIndv/";
my $toDir = "$baseDir/concatFilterIndv/temp/$keepOneList[$indexTwo]";

fmove ($toDir, $fromDir) or die print $logKeepOneFH "keepOne species for analysis not present in Temp directory!\n";


print $logKeepOneFH "$toDir was just moved to $fromDir for analysis\n\n";


my @finalArray;
my $counterTwo = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]" or die "Cannot open concatFilter $speciesArray[$counterTwo]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counterTwo++;

}else{
#Do Nothing
$counterTwo++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/KeepOne_GroupFasta_$keepOneList[$indexTwo]" or die "Cannot open KeepOne_GroupFasta_$keepOneList[$indexTwo]: $!";
print $singleFilePush @finalArray;
close $singleFilePush;


#Make directoy for each keepOne scenario

my $dirKeepOne = "$baseDir/treeRAxML/keepOne/KeepOne_$keepOneList[$indexTwo]";
unless(-d $dirKeepOne)
{
mkdir $dirKeepOne or die "Unable to create $dirKeepOne\n";
}

print $logKeepOneFH "Group.fasta with $keepOneList[$indexTwo] included will be called:\n"; 
print $logKeepOneFH "KeepOne_GroupFasta_$keepOneList[$indexTwo]\n\n\n";

undef @finalArray;
undef $counterTwo;

my $fromDir2 = "$baseDir/concatFilterIndv/$keepOneList[$indexTwo]";
my $toDir2 = "$baseDir/concatFilterIndv/temp/";

fmove ($fromDir2, $toDir2) or die print $logKeepOneFH "keepOne species unable to be moved back to Temp directory!\n";


$indexTwo ++;

}

#Move all files back from temp to main folder.

print $logKeepOneFH "Moving files back to fasta directory.\n\n";


my $indexThree = 0;
foreach (1 .. $numKeepOne) {

my $fromDir = "$baseDir/concatFilterIndv/";
my $toDir = "$baseDir/concatFilterIndv/temp/$keepOneList[$indexThree]";

fmove ($toDir, $fromDir) or die print $logKeepOneFH "Cannot move keepOne Species back to fasta folder!\n";

print $logKeepOneFH "$toDir was just moved to $fromDir\n";

$indexThree ++;
}

print $logKeepOneFH "All files were moved back to fasta directory.\n\n\n\n";




if(-e $directoryKeepOneTemp) {
print $logKeepOneFH "Directory '$directoryKeepOneTemp' still exists.\n";
} else {
print $logKeepOneFH "Directory '$directoryKeepOneTemp' has been removed.\n";
}

close $logKeepOneFH;



###########################
#Remove Group Permutations#
###########################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {

my $directoryRemoveGroup = "$baseDir/treeRAxML/removeGroup";
unless(-d $directoryRemoveGroup)
{
mkdir $directoryRemoveGroup or die "Unable to create $directoryRemoveGroup\n";
}


#Create temp directory to hold moved files

my $directoryRemoveGroupTemp = "$baseDir/concatFilterIndv/temp";
unless(-d $directoryRemoveGroupTemp)
{
mkdir $directoryRemoveGroupTemp or die "Unable to create $directoryRemoveGroupTemp\n";
}


#Generate RemoveGroup Log

open my $logRemoveGroupFH, '>', "$baseDir/logFiles/removeGroupLog.txt" or die "Cannot open removeGroupLog.txt: $!";

my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp;
    s/\r//g;	
    push @removeGroupList,$_;
} 
close $removeGroupFH;

my $numRemoveGroup =  scalar @removeGroupList;



#Move All Keep One Files to temp folder.

print $logRemoveGroupFH "Files being moved for initial setup:\n";

my $index = 0;
foreach (1 .. $numRemoveGroup) {

my $fromDir = "$baseDir/concatFilterIndv/$removeGroupList[$index]";
my $toDir = "$baseDir/concatFilterIndv/temp";

fmove ($fromDir, $toDir) or die print $logRemoveGroupFH "All keepOne species were not present in directory!\n";

print $logRemoveGroupFH "$fromDir was just moved to $toDir\n";

$index ++;
}


print $logRemoveGroupFH "\n\n";


#Start Analysis with group removed

my @finalArray;
$counter = 0;

foreach (1 .. $numSpecies1){
	
my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}	

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveGroup_Group.fasta" or die "Cannot open RemoveGroup_Group.fasta: $!";
print $singleFilePush @finalArray;
close $singleFilePush;


#Make directoy for each keepOne scenario

my $dirRemoveGroup = "$baseDir/treeRAxML/removeGroup";
unless(-d $dirRemoveGroup)
{
mkdir $dirRemoveGroup or die "Unable to create $dirRemoveGroup\n";
}


#Move all files back from temp to main folder.

print $logRemoveGroupFH "Moving files back to fasta directory.\n\n";


my $indexTwo = 0;
foreach (1 .. $numRemoveGroup) {

my $fromDir = "$baseDir/concatFilterIndv";
my $toDir = "$baseDir/concatFilterIndv/temp/$removeGroupList[$indexTwo]";

fmove ($toDir, $fromDir) or die print $logRemoveGroupFH "Cannot move keepOne Species back to fasta folder!\n";

print $logRemoveGroupFH "$toDir was just moved to $fromDir\n";

$indexTwo ++;
}

print $logRemoveGroupFH "All files were moved back to fasta directory.\n\n\n\n";


rmdir $directoryRemoveGroupTemp;

if(-e $directoryRemoveGroupTemp) {
print $logRemoveGroupFH "Directory '$directoryRemoveGroupTemp' still exists.\n";
} else {
print $logRemoveGroupFH "Directory '$directoryRemoveGroupTemp' has been removed.\n";
}

close $logRemoveGroupFH;



################################
#Standard Run - No Permutations#
################################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {


my @finalArray;
$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/groupFasta.fasta" or die "Cannot open groupFasta.fasta: $!";

print $singleFilePush @finalArray;


########################################
#Multiple Permutations Selected - Error#
########################################

} else {

print $logconcatFH "You cannot run multiple permutations at once!!\n";
print $logconcatFH "Fix this Error and rerun the pipeline!\n";
die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");

}


###########
#FastTree######################################################################################################################
###########

} elsif ($ENV{'PoPipe_RAxMLTree'} == 0 && $ENV{'PoPipe_fastTree'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {

#########################
#Remove One Permutations#
#########################

if ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0) {

my $directoryRemoveOne = "$baseDir/fastTree/removeOne";
unless(-d $directoryRemoveOne)
{
mkdir $directoryRemoveOne or die "Unable to create $directoryRemoveOne\n";
}

#generate log
open my $logRemoveOneFH, '>', "$baseDir/logFiles/removeOneLog.txt" or die "Cannot open removeOneLog.txt: $!";

my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp;
	s/\r//g;
    push @removeOneList,$_;
} 
close $removeOneFH;


my $numRemoveOne =  scalar @removeOneList;
my $index = 0;

print $logRemoveOneFH "Number of Species to be removed: $numRemoveOne\n\n";


foreach (1 .. $numRemoveOne){

print $logRemoveOneFH "Currently removing $removeOneList[$index] from analysis\n\n";

my $fileRename = "rename.txt";

rename ("$baseDir/concatFilterIndv/$removeOneList[$index]", "$baseDir/concatFilterIndv/$fileRename");


my @finalArray;
my $counterTwo = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]" or die "Cannot open concatFilter $speciesArray[$counterTwo]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counterTwo++;

}else{
#Do Nothing
$counterTwo++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveOne_GroupFasta_$removeOneList[$index]" or die "Cannot open RemoveOne_GroupFasta_$removeOneList[$index]: $!";
print $singleFilePush @finalArray;
close $singleFilePush;

#Make directoy for each RemoveOne scenario

my $dirRemoveOne = "$baseDir/fastTree/removeOne/RemoveOne_$removeOneList[$index]";
unless(-d $dirRemoveOne)
{
mkdir $dirRemoveOne or die "Unable to create $dirRemoveOne\n";
}

rename ("$baseDir/concatFilterIndv/$fileRename", "$baseDir/concatFilterIndv/$removeOneList[$index]");

print $logRemoveOneFH "Group.fasta without $removeOneList[$index] included will be called:\n"; 
print $logRemoveOneFH "RemoveOne_GroupFasta_$removeOneList[$index]\n\n";


undef @finalArray;
undef $counterTwo;

$index ++;

}
 close $logRemoveOneFH;


#######################
#Keep One Permutations#
#######################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {


my $directoryKeepOne = "$baseDir/fastTree/keepOne";
unless(-d $directoryKeepOne)
{
mkdir $directoryKeepOne or die "Unable to create $directoryKeepOne\n";
}

#Create temp directory to fold moved files

my $directoryKeepOneTemp = "$baseDir/concatFilterIndv/temp";
unless(-d $directoryKeepOneTemp)
{
mkdir $directoryKeepOneTemp or die "Unable to create $directoryKeepOneTemp\n";
}


#Generate KeepOne Log

open my $logKeepOneFH, '>', "$baseDir/logFiles/keepOneLog.txt" or die "Cannot open keepOneLog.txt: $!";

my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp;
	s/\r//g;
    push @keepOneList,$_;
} 
close $keepOneFH;

my $numKeepOne =  scalar @keepOneList;



#Move All Keep One Files to temp folder.

print $logKeepOneFH "Files being moved for initial setup:\n";

my $index = 0;
foreach (1 .. $numKeepOne) {
	
my $fromDir = "$baseDir/concatFilterIndv/$keepOneList[$index]";
my $toDir = "$baseDir/concatFilterIndv/temp/";

fmove($fromDir, $toDir) or die "All keepOne species were not present in directory!: $!";

print $logKeepOneFH "$fromDir was just moved to $toDir\n";

$index ++;
}


print $logKeepOneFH "\n\n";


#Start multiple group.fasta creations
my $indexTwo = 0;

foreach (1 .. $numKeepOne) {

#Move the Keep One back to run with analysis
print $logKeepOneFH "keepOne file being moved back into analysis:\n";

my $fromDir = "$baseDir/concatFilterIndv/";
my $toDir = "$baseDir/concatFilterIndv/temp/$keepOneList[$indexTwo]";

fmove ($toDir, $fromDir) or die print $logKeepOneFH "keepOne species for analysis not present in Temp directory!\n";

print $logKeepOneFH "$toDir was just moved to $fromDir for analysis\n\n";


my @finalArray;
my $counterTwo = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counterTwo]" or die "Cannot open concatFilter $speciesArray[$counterTwo]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counterTwo++;

}else{
#Do Nothing
$counterTwo++;

print $logKeepOneFH "\n\n";
print $logKeepOneFH "ERROR";
print $logKeepOneFH "\n\n";


}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/KeepOne_GroupFasta_$keepOneList[$indexTwo]" or die "Cannot open KeepOne_GroupFasta_$keepOneList[$indexTwo]: $!";
print $singleFilePush @finalArray;
close $singleFilePush;


#Make directoy for each keepOne scenario

my $dirKeepOne = "$baseDir/fastTree/keepOne/KeepOne_$keepOneList[$indexTwo]";
unless(-d $dirKeepOne)
{
mkdir $dirKeepOne or die "Unable to create $dirKeepOne\n";
}

print $logKeepOneFH "Group.fasta with $keepOneList[$indexTwo] included will be called:\n"; 
print $logKeepOneFH "KeepOne_GroupFasta_$keepOneList[$indexTwo]\n\n";

undef @finalArray;
undef $counterTwo;

my $fromDir2 = "$baseDir/concatFilterIndv/$keepOneList[$indexTwo]";
my $toDir2 = "$baseDir/concatFilterIndv/temp/";

fmove ($fromDir2, $toDir2) or die print $logKeepOneFH "keepOne species unable to be moved back to Temp directory!\n";


$indexTwo ++;

}

#Move all files back from temp to main folder.

print $logKeepOneFH "Moving files back to fasta directory.\n\n";


my $indexThree = 0;
foreach (1 .. $numKeepOne) {

my $fromDir3 = "$baseDir/concatFilterIndv/";
my $toDir3 = "$baseDir/concatFilterIndv/temp/$keepOneList[$indexThree]";

fmove ($toDir3, $fromDir3) or die print $logKeepOneFH "Cannot move keepOne Species back to fasta folder!\n";

print $logKeepOneFH "$toDir3 was just moved to $fromDir3\n";

$indexThree ++;
}

print $logKeepOneFH "All files were moved back to fasta directory.\n\n\n\n";


my $tempDirKO = "$baseDir/concatFilterIndv/temp";
rmtree $tempDirKO;

if(-e $directoryKeepOneTemp) {
print $logKeepOneFH "Directory '$directoryKeepOneTemp' still exists.\n";
} else {
print $logKeepOneFH "Directory '$directoryKeepOneTemp' has been removed.\n";
}

close $logKeepOneFH;



###########################
#Remove Group Permutations#
###########################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {

my $directoryRemoveGroup = "$baseDir/fastTree/removeGroup";
unless(-d $directoryRemoveGroup)
{
mkdir $directoryRemoveGroup or die "Unable to create $directoryRemoveGroup\n";
}


#Create temp directory to fold moved files

my $directoryRemoveGroupTemp = "$baseDir/concatFilterIndv/temp";
unless(-d $directoryRemoveGroupTemp)
{
mkdir $directoryRemoveGroupTemp or die "Unable to create $directoryRemoveGroupTemp\n";
}


#Generate RemoveGroup Log

open my $logRemoveGroupFH, '>', "$baseDir/logFiles/removeGroupLog.txt" or die "Cannot open removeGroupLog.txt: $!";

my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp;
	s/\r//g;
    push @removeGroupList,$_;
} 
close $removeGroupFH;

my $numRemoveGroup =  scalar @removeGroupList;



#Move All Keep One Files to temp folder.

print $logRemoveGroupFH "Files being moved for initial setup:\n";

my $index = 0;
foreach (1 .. $numRemoveGroup) {

my $fromDir = "$baseDir/concatFilterIndv/$removeGroupList[$index]";
my $toDir = "$baseDir/concatFilterIndv/temp";

fmove ($fromDir, $toDir) or die print $logRemoveGroupFH "All keepOne species were not present in directory!\n";

print $logRemoveGroupFH "$fromDir was just moved to $toDir\n";

$index ++;
}


print $logRemoveGroupFH "\n\n";


#Start Analysis with group removed

my @finalArray;
$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveGroup_Group.fasta" or die "Cannot open RemoveGroup_Group.fasta: $!";
print $singleFilePush @finalArray;
close $singleFilePush;


#Make directoy for each keepOne scenario

my $dirRemoveGroup = "$baseDir/fastTree/removeGroup";
unless(-d $dirRemoveGroup)
{
mkdir $dirRemoveGroup or die "Unable to create $dirRemoveGroup\n";
}


#Move all files back from temp to main folder.

print $logRemoveGroupFH "Moving files back to fasta directory.\n\n";


my $indexTwo = 0;
foreach (1 .. $numRemoveGroup) {

my $fromDir = "$baseDir/concatFilterIndv";
my $toDir = "$baseDir/concatFilterIndv/temp/$removeGroupList[$indexTwo]";

fmove ($toDir, $fromDir) or die print $logRemoveGroupFH "Cannot move keepOne Species back to fasta folder!\n";

print $logRemoveGroupFH "$toDir was just moved to $fromDir\n";

$indexTwo ++;
}

print $logRemoveGroupFH "All files were moved back to fasta directory.\n\n\n\n";


rmdir $directoryRemoveGroupTemp;

if(-e $directoryRemoveGroupTemp) {
print $logRemoveGroupFH "Directory '$directoryRemoveGroupTemp' still exists.\n";
} else {
print $logRemoveGroupFH "Directory '$directoryRemoveGroupTemp' has been removed.\n";
}

close $logRemoveGroupFH;



################################
#Standard Run - No Permutations#
################################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {


my @finalArray;
$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
$counter++;
#Do Nothing
}

}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/groupFasta.fasta" or die "Cannot open groupFasta.fasta: $!";

print $singleFilePush @finalArray;


########################################
#Multiple Permutations Selected - Error#
########################################

} else {

print $logconcatFH "You cannot run multiple permutations at once!!\n";
print $logconcatFH "Fix this Error and rerun the pipeline!\n";
die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");

}

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################



########################################
#Standard Run - Full Rerun Permutations#
########################################


} elsif ($ENV{'PoPipe_RAxMLTree'} == 1 && $ENV{'PoPipe_fastTree'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1){


################
#Keep One - FRR#
################

if ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1) {

my $dirKeepOneFR = "$baseDir/treeRAxML/keepOneFullReRun/";
unless(-d $dirKeepOneFR)
{
mkdir $dirKeepOneFR or die "Unable to create $dirKeepOneFR\n";
}

my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp;
    s/\r//g;	
    push @keepOneList,$_;
} 
close $keepOneFH;

#Move the Keep One back to temp folder

my $fromDir = "$baseDir/fasta/$keepOneList[0]";
my $toDir = "$baseDir/fasta/temp/";

fmove ($fromDir, $toDir) or die "keepOne species for analysis not present in Temp directory!\n";

my $indexTwo = 0;

#Make directoy for each keepOne scenario

my $dirKeepOne = "$baseDir/treeRAxML/keepOneFullReRun/KeepOne_$keepOneList[$indexTwo]";
unless(-d $dirKeepOne)
{
mkdir $dirKeepOne or die "Unable to create $dirKeepOne\n";
}

foreach (1 .. $numSpecies1){


$indexTwo ++;

my @finalArray;

$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}


for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/KeepOne_$keepOneList[0]_Group.fasta" or die "Cannot open groupFasta.fasta: $!";

print $singleFilePush @finalArray;

}


####################
#Remvoe Group - FRR#
####################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 1){


my $dirRemoveGroupFR = "$baseDir/treeRAxML/removeGroupFullReRun/";
unless(-d $dirRemoveGroupFR)
{
mkdir $dirRemoveGroupFR or die "Unable to create $dirRemoveGroupFR\n";
}


my @finalArray;
$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveGroup_FullReRunGroup.fasta" or die "Cannot open RemoveGroup_FullReRunGroup.fasta: $!";

print $singleFilePush @finalArray;


##################
#Remvoe One - FRR#
##################

} elsif ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1){



open my $logRemoveOneFH, '>', "$baseDir/logFiles/removeOneLog.txt" or die "Cannot open removeoneLog.txt: $!";

print $logRemoveOneFH "Running RemoveOne - RAxMLTree\n";

my $dirRemoveOneFR = "$baseDir/treeRAxML/removeOneFullReRun/";
unless(-d $dirRemoveOneFR)
{
mkdir $dirRemoveOneFR or die "Unable to create $dirRemoveOneFR\n";
}

my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp; 
	s/\r//g;
    push @removeOneList,$_;
} 
close $removeOneFH;

#Make directoy for each removeOne scenario

my $dirRemoveOne = "$baseDir/treeRAxML/removeOneFullReRun/RemoveOne_$removeOneList[0]";
unless(-d $dirRemoveOne)
{
mkdir $dirRemoveOne or die "Unable to create $dirRemoveOne\n";
}

foreach (1 .. $numSpecies1){

my @finalArray;

$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveOne_$removeOneList[0]_Group.fasta" or die "Cannot open groupFasta.fasta: $!";

print $singleFilePush @finalArray;

}

} else {

print $logconcatFH "You cannot run multiple permutations at once!!\n";
print $logconcatFH "Fix this Error and rerun the pipeline!\n";
die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");

}








###########
#Fast Tree#
###########


} elsif ($ENV{'PoPipe_RAxMLTree'} == 0 && $ENV{'PoPipe_fastTree'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 1){


################
#Keep One - FRR#
################

if ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1) {


my $dirKeepOneFR = "$baseDir/fastTree/keepOneFullReRun/";
unless(-d $dirKeepOneFR)
{
mkdir $dirKeepOneFR or die "Unable to create $dirKeepOneFR\n";
}

my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp; 
	s/\r//g;
    push @keepOneList,$_;
} 
close $keepOneFH;

#Move the Keep One back to temp folder

my $fromDir = "$baseDir/fasta/$keepOneList[0]";
my $toDir = "$baseDir/fasta/temp/";

fmove ($fromDir, $toDir) or die "keepOne species for analysis not present in Temp directory! ERROR LN1346\n";

my $indexTwo = 0;

#Make directoy for each keepOne scenario

my $dirKeepOne = "$baseDir/fastTree/keepOneFullReRun/KeepOne_$keepOneList[$indexTwo]";
unless(-d $dirKeepOne)
{
mkdir $dirKeepOne or die "Unable to create $dirKeepOne\n";
}

foreach (1 .. $numSpecies1){


$indexTwo ++;

my @finalArray;

$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/KeepOne_$keepOneList[0]_Group.fasta" or die "Cannot open groupFasta.fasta: $!";

print $singleFilePush @finalArray;

}


####################
#Remove Group - FRR#
####################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 1){


my $dirRemoveGroupFR = "$baseDir/fastTree/removeGroupFullReRun/";
unless(-d $dirRemoveGroupFR)
{
mkdir $dirRemoveGroupFR or die "Unable to create $dirRemoveGroupFR\n";
}


my @finalArray;
$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveGroup_FullReRunGroup.fasta" or die "Cannot open RemoveGroup_FullReRunGroup.fasta: $!";

print $singleFilePush @finalArray;



##################
#Remvoe One - FRR#
##################

} elsif ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1){



open my $logRemoveOneFH, '>', "$baseDir/logFiles/removeOneLog.txt" or die "Cannot open removeOneLog.txt: $!";

print $logRemoveOneFH "Running RemoveOne - fastTree\n";

my $dirRemoveOneFR = "$baseDir/fastTree/removeOneFullReRun/";
unless(-d $dirRemoveOneFR)
{
mkdir $dirRemoveOneFR or die "Unable to create $dirRemoveOneFR\n";
}

my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp;
    s/\r//g;	
    push @removeOneList,$_;
} 
close $removeOneFH;

#Make directoy for each removeOne scenario

my $dirRemoveOne = "$baseDir/fastTree/removeOneFullReRun/RemoveOne_$removeOneList[0]";
unless(-d $dirRemoveOne)
{
mkdir $dirRemoveOne or die "Unable to create $dirRemoveOne\n";
}

foreach (1 .. $numSpecies1){

my @finalArray;

$counter = 0;

foreach (1 .. $numSpecies1){

my $fileExistsPath = "$baseDir/concatFilterIndv/$speciesArray[$counter]";
	
if (-e $fileExistsPath){ 

open my $finalSpeciesFH, '<', "$baseDir/concatFilterIndv/$speciesArray[$counter]" or die "Cannot open concatFilter $speciesArray[$counter]: $!";

while(<$finalSpeciesFH>) { 
    chomp;
    push @finalArray, "$_\n";
} 
close $finalSpeciesFH;

$counter++;

}else{
#Do Nothing
$counter++;
}
}

for (@finalArray) {
s/^\s+//;
}

open my $singleFilePush, '>', "$baseDir/concatFilterGroup/RemoveOne_$removeOneList[0]_Group.fasta" or die "Cannot open RemoveOne_$removeOneList[0]_Group.fasta: $!";

print $singleFilePush @finalArray;

}


} else {

print $logconcatFH "You cannot run multiple permutations at once!!\n";
print $logconcatFH "Fix this Error and rerun the pipeline!\n";
die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");

}


}else {

die ("Error");
}