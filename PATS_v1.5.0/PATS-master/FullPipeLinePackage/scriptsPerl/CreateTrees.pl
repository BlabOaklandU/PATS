#!usr/bin/perl
use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $numThreads = $ENV{'cpusPerTask'};

my $myGammaChoice = $ENV{'gammaChoice'};
my $myfastTreeModel = $ENV{'fastTreeModel'};

my $trueGamma;
if ($myGammaChoice == 1) {$trueGamma = '-gamma' }else{
#do nothing
}


if($ENV{'PoPipe_RAxMLTree'} == 1 && $ENV{'PoPipe_fastTree'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {


#############################
#Running RAxML Tree Creation#
#############################

#generate RAxML log
open my $logRAxMLFH, '>', "$baseDir/logFiles/RAxMLLog.txt" or die "Cannot open RAxMLLog.txt: $!";

###################################################################################################################

if ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp; 
    push @removeOneList,$_;
} 
close $removeOneFH;

my $numRemoveOne =  scalar @removeOneList;

print $logRAxMLFH "You are running a Remove One Permutation\n\n";


foreach (1 .. $numRemoveOne){
	
my $currenSpecies = $removeOneList[$index];
	
my $filename = "$baseDir/concatFilterGroup/RemoveOne_GroupFasta_${currenSpecies}";
my $fileOutputNWK = "$baseDir/treeRAxML/removeOne/RemoveOne_${currenSpecies}.nwk";
my $newName = "RemoveOne_${currenSpecies}";

print $logRAxMLFH "Starting on RemoveOne_GroupFasta_${currenSpecies} tree creation\n";
my $timestamp0 = localtime(time);
print $logRAxMLFH "Start Time: $timestamp0\n";

system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");			
print $logRAxMLFH "Finished RemoveOne_GroupFasta_${currenSpecies} tree creation\n";
my $timestamp1 = localtime(time);
print $logRAxMLFH "End Time: $timestamp1\n\n";

$index ++;
}

#################################################################################################################

}elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp; 
    push @keepOneList,$_;
} 
close $keepOneFH;

my $numKeepOne =  scalar @keepOneList;

print $logRAxMLFH "You are running a Keep One Permutation\n\n";

my $currentSpecies;

foreach (1 .. $numKeepOne){
	
my $currenSpecies = $keepOneList[$index];

my $filename = "$baseDir/concatFilterGroup/KeepOne_GroupFasta_${currenSpecies}";
my $fileOutputNWK = "$baseDir/treeRAxML/keepOne/KeepOne_${currenSpecies}.nwk";
my $newName = "KeepOne_${currenSpecies}";

print $logRAxMLFH "Starting on KeepOne_GroupFasta_${currenSpecies} tree creation\n";
my $timestamp0 = localtime(time);
print $logRAxMLFH "Start Time: $timestamp0\n";

system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");		

print $logRAxMLFH "Finished KeepOne_GroupFasta_${currenSpecies} tree creation\n";
my $timestamp1 = localtime(time);
print $logRAxMLFH "End Time: $timestamp1\n\n";

undef $currenSpecies;

$index ++;
}


#################################################################################################################

}elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1) {

my $index = 0;
my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp; 
    push @removeGroupList,$_;
} 
close $removeGroupFH;

my $numremoveGroup =  scalar @removeGroupList;

my $filename = "$baseDir/concatFilterGroup/RemoveGroup_Group.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/removeGroup/RemoveGroup_Group.nwk";
my $newName = "RemoveGroup_Group.fasta";

print $logRAxMLFH "You are running a Remove Group Permutation\n";
my $timestamp0 = localtime(time);
print $logRAxMLFH "Start Time: $timestamp0\n";

system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");

my $timestamp1 = localtime(time);
print $logRAxMLFH "End Time: $timestamp1\n\n";			

$index ++;

print $logRAxMLFH "Finished RemoveGroup_Group.fasta tree creation\n\n";


#################################################################################################################


} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0) {

print $logRAxMLFH "You are running the standard tree creation, no permutations selected.\n";

#Create base directory for standard [no permutations] run
my $directoryBaseNWK = "$baseDir/treeRAxML/standard";
unless(-d $directoryBaseNWK)
{
mkdir $directoryBaseNWK or die "Unable to create $directoryBaseNWK\n";
}


my $filename = "$baseDir/concatFilterGroup/groupFasta.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/standard.nwk";
my @newName;

print $logRAxMLFH "Starting groupFasta.fasta tree creation.\n";
my $timestamp0 = localtime(time);
print $logRAxMLFH "Start Time: $timestamp0\n";

system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n @newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");	

print $logRAxMLFH "Finished groupFasta.fasta tree creation.\n\n";
my $timestamp1 = localtime(time);
print $logRAxMLFH "End Time: $timestamp1\n\n";	


###################################################################################################################

 } else {

die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");
}


###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################


}elsif ($ENV{'PoPipe_RAxMLTree'} == 0 && $ENV{'PoPipe_fastTree'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 0) {

################################
#Running fastTree Tree Creation#
################################

#generate fastTree log
open my $logfastTreeFH, '>', "$baseDir/logFiles/fastTreeLog.txt" or die "Cannot open fastTreeLog.txt: $!";

###################################################################################################################

if ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp;
	s/\r//g;
    push @removeOneList,$_;
} 
close $removeOneFH;

my $numRemoveOne =  scalar @removeOneList;

print $logfastTreeFH "You are running a Remove One Permutation\n\n";

foreach (1 .. $numRemoveOne){
	
print $logfastTreeFH "Starting on RemoveOne_GroupFasta_$removeOneList[$index] tree creation\n\n";

my $filename = "$baseDir/concatFilterGroup/RemoveOne_GroupFasta_$removeOneList[$index]";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree";

system("$fastaTreeMPPath/FastTreeMP $trueGamma -$myfastTreeModel -log fastTree/removeOne/RemoveOne_$removeOneList[$index]/fastTreeLog_$removeOneList[$index].txt $filename > $baseDir/fastTree/removeOne/RemoveOne_$removeOneList[$index]/Tree_$removeOneList[$index].nwk");

print $logfastTreeFH "Finished RemoveOne_GroupFasta_$removeOneList[$index] tree creation\n\n";

$index ++;
}

#################################################################################################################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp;
	s/\r//g;
    push @keepOneList,$_;
} 
close $keepOneFH;

my $numKeepOne =  scalar @keepOneList;

print $logfastTreeFH "You are running a Keep One Permutation\n\n";

foreach (1 .. $numKeepOne){
	
print $logfastTreeFH "Starting on KeepOne_GroupFasta_$keepOneList[$index] tree creation\n";

my $filename = "$baseDir/concatFilterGroup/KeepOne_GroupFasta_$keepOneList[$index]";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree";

system("$fastaTreeMPPath/FastTreeMP $trueGamma -$myfastTreeModel -log fastTree/keepOne/KeepOne_$keepOneList[$index]/fastTreeLog_$keepOneList[$index].txt $filename > $baseDir/fastTree/keepOne/KeepOne_$keepOneList[$index]/Tree_$keepOneList[$index].nwk");

print $logfastTreeFH "Finished KeepOne_GroupFasta_$keepOneList[$index] tree creation\n\n";

$index ++;
}


#################################################################################################################

}elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1) {

my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp;
	s/\r//g;
    push @removeGroupList,$_;
} 
close $removeGroupFH;

my $numremoveGroup =  scalar @removeGroupList;

print $logfastTreeFH "You are running a Remove Group Permutation\n";

my $filename = "$baseDir/concatFilterGroup/RemoveGroup_Group.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree";


system("$fastaTreeMPPath/FastTreeMP $trueGamma -$myfastTreeModel -log fastTree/removeGroup/fastTreeLog.txt $filename > $baseDir/fastTree/removeGroup/Tree_removeGroup.nwk");			

print $logfastTreeFH "Finished RemoveGroup_Group.fasta tree creation\n\n";


#################################################################################################################


} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0) {

print $logfastTreeFH "You are running the standard tree creation, no permutations selected.\n";

#Create base directory for standard [no permutations] run
my $directoryBaseNWK = "$baseDir/fastTree/standard";
unless(-d $directoryBaseNWK)
{
mkdir $directoryBaseNWK or die "Unable to create $directoryBaseNWK\n";
}


my $filename = "$baseDir/concatFilterGroup/groupFasta.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree";

print $logfastTreeFH "Starting groupFasta.fasta tree creation.\n";

#chdir $fastaTreeMPPath;

system("$fastaTreeMPPath/FastTreeMP $trueGamma -$myfastTreeModel -log fastTree/standard/fastTreeLog.txt $filename > $baseDir/fastTree/standard/Tree_standard.nwk");

print $logfastTreeFH "Finished groupFasta.fasta tree creation.\n\n";		


###################################################################################################################

} else {

die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");
}

} else {

die ("You cannot run fastTree and RAxML at the same time. Double check the poPipeSettings file.");

}