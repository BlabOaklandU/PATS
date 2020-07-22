#!usr/bin/perl
use strict;
use warnings;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $numThreads = $ENV{'PoPipe_threadsToUse'};

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

print $logRAxMLFH "You are running a Remove One Permutation\n";

print $logRAxMLFH "Starting on RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n";

foreach (1 .. $numRemoveOne){

my $filename = "$baseDir/concatFilterGroup/RemoveOne_$removeOneList[$index]_Group.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/removeOne/RemoveOne_$removeOneList[$index]";
my $newName = "RemoveOne_$removeOneList[$index]";


system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");			

$index ++;

print $logRAxMLFH "Finished RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n\n";

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

print $logRAxMLFH "You are running a Keep One Permutation\n";

print $logRAxMLFH "Starting on KeepOne_$keepOneList[$index]_Group.fasta tree creation\n";

foreach (1 .. $numKeepOne){

my $filename = "$baseDir/concatFilterGroup/KeepOne_$keepOneList[$index]_Group.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/keepOne/KeepOne_$keepOneList[$index]";
my $newName = "KeepOne_$keepOneList[$index]";


system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");			

$index ++;

print $logRAxMLFH "Finished KeepOne_$keepOneList[$index]_Group.fasta tree creation\n\n";
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

print $logRAxMLFH "You are running a Remove Group Permutation\n";

my $filename = "$baseDir/concatFilterGroup/RemoveGroup_Group.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/removeGroup/RemoveGroup_Group.fasta";
my $newName = "RemoveGroup_Group.fasta";


system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");			

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
my $fileOutputNWK = "$baseDir/treeRAxML/standard";
my @newName;

print $logRAxMLFH "Starting groupFasta.fasta tree creation.\n";

system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n @newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");	

print $logRAxMLFH "Finished groupFasta.fasta tree creation.\n\n";		


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
    push @removeOneList,$_;
} 
close $removeOneFH;

my $numRemoveOne =  scalar @removeOneList;

print $logfastTreeFH "You are running a Remove One Permutation\n";

print $logfastTreeFH "Starting on RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n";

foreach (1 .. $numRemoveOne){

my $filename = "$baseDir/concatFilterGroup/RemoveOne_$removeOneList[$index]_Group.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree";


system("$fastaTreeMPPath/FastTreeMP -lg -log fastTree/removeOne/RemoveOne_$removeOneList[$index]/fastTreeLog_$removeOneList[$index].txt $filename");

$index ++;

print $logfastTreeFH "Finished RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n\n";

}

#################################################################################################################

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp; 
    push @keepOneList,$_;
} 
close $keepOneFH;

my $numKeepOne =  scalar @keepOneList;

print $logfastTreeFH "You are running a Keep One Permutation\n";

print $logfastTreeFH "Starting on KeepOne_$keepOneList[$index]_Group.fasta tree creation\n";

foreach (1 .. $numKeepOne){

my $filename = "$baseDir/concatFilterGroup/KeepOne_$keepOneList[$index]_Group.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree";


system("$fastaTreeMPPath/FastTreeMP -lg -log fastTree/keepOne/KeepOne_$keepOneList[$index]/fastTreeLog_$keepOneList[$index].txt $filename");


$index ++;

print $logfastTreeFH "Finished KeepOne_$keepOneList[$index]_Group.fasta tree creation\n\n";
}


#################################################################################################################

}elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1) {

my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp; 
    push @removeGroupList,$_;
} 
close $removeGroupFH;

my $numremoveGroup =  scalar @removeGroupList;

print $logfastTreeFH "You are running a Remove Group Permutation\n";

my $filename = "$baseDir/concatFilterGroup/RemoveGroup_Group.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree";


system("$fastaTreeMPPath/FastTreeMP -lg -log fastTree/removeGroup/fastTreeLog.txt $filename");			

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

system("$fastaTreeMPPath/FastTreeMP -lg -log fastTree/standard/fastTreeLog.txt $filename");

print $logfastTreeFH "Finished groupFasta.fasta tree creation.\n\n";		


###################################################################################################################

} else {

die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");
}




} else {

die ("You cannot run fastTree and RAxML at the same time. Double check the poPipeSettings file.");

}