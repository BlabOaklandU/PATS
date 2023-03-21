#!usr/bin/perl
use strict;
use warnings;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $numThreads = $ENV{'PoPipe_threadsToUse'};


my $myGammaChoice = $ENV{'gammaChoice'};
my $myfastTreeModel = $ENV{'fastTreeModel'};

my $trueGamma;
if ($myGammaChoice == 1) {$trueGamma = '-gamma' }else{
#do nothing
}

###################################################################################################
#
#Full Rerun Permutations
#
###################################################################################################


############
#RAxML Tree#
############

if ($ENV{'PoPipe_RAxMLTree'} == 1 && $ENV{'PoPipe_fastTree'} == 0 && $ENV{'PoPipe_rerunFullAnalysis'} == 1) {


#generate RAxML log
open my $logRAxMLFH, '>', "$baseDir/logFiles/RAxMLLog.txt" or die "Cannot open RAxMLLog.txt: $!";



############
#Remove One#
############

if ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp; 
    push @removeOneList,$_;
} 
close $removeOneFH;

print $logRAxMLFH "You are running a Remove One Permutation\n";

print $logRAxMLFH "Starting on RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n";

my $filename = "$baseDir/concatFilterGroup/RemoveOne_$removeOneList[$index]_Group.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/removeOneFullReRun/RemoveOne_$removeOneList[$index]";
my $newName = "RemoveOne_$removeOneList[$index]";


system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");			

print $logRAxMLFH "Finished RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n\n";




##########
#Keep One#
##########

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 1 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @keepOneList;

open my $keepOneFH, '<', "$baseDir/permutations/keepOne.txt" or die "Cannot open keepOne.txt: $!";
while(<$keepOneFH>) { 
    chomp; 
    push @keepOneList,$_;
} 
close $keepOneFH;


print $logRAxMLFH "You are running a Keep One Permutation\n";

print $logRAxMLFH "Starting on KeepOne_$keepOneList[$index]_Group.fasta tree creation\n";


my $filename = "$baseDir/concatFilterGroup/KeepOne_$keepOneList[$index]_Group.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/keepOneFullReRun/KeepOne_$keepOneList[$index]";
my $newName = "KeepOne_$keepOneList[$index]";


system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");			


print $logRAxMLFH "Finished KeepOne_$keepOneList[$index]_Group.fasta tree creation\n\n";


##############
#Remove Group#
##############

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1) {

my $index = 0;
my @removeGroupList;

open my $removeGroupFH, '<', "$baseDir/permutations/removeGroup.txt" or die "Cannot open removeGroup.txt: $!";
while(<$removeGroupFH>) { 
    chomp; 
    push @removeGroupList,$_;
} 
close $removeGroupFH;

print $logRAxMLFH "You are running a Remove Group Permutation\n";

print $logRAxMLFH "Starting RemoveGroup_FullReRunGroup.fasta tree creation\n";



my $filename = "$baseDir/concatFilterGroup/RemoveGroup_FullReRunGroup.fasta";
my $fileOutputNWK = "$baseDir/treeRAxML/removeGroupFullReRun";
my $newName = "RemoveGroup_FullReRunGroup";


system("$baseDir/externalSoftware/RAxML/raxmlHPC-PTHREADS-AVX -T $numThreads -s $filename -n $newName.nwk -w $fileOutputNWK -k -f a -p 95137 -x 84621 -# 100 -m PROTGAMMAAUTO");			



print $logRAxMLFH "Finished RemoveGroup_FullReRunGroup.fasta tree creation\n\n";



} else {

die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");

}





##########
#FastTree#
##########

} elsif($ENV{'PoPipe_RAxMLTree'} == 0 && $ENV{'PoPipe_fastTree'} == 1 && $ENV{'PoPipe_rerunFullAnalysis'} == 1) {


#generate fastTree log
open my $logfastTreeFH, '>', "$baseDir/logFiles/fastTreeLog.txt" or die "Cannot open fastTreeLog.txt: $!";


############
#Remove One#
############

if ($ENV{'PoPipe_removeOne'} == 1 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 0) {

my $index = 0;
my @removeOneList;

open my $removeOneFH, '<', "$baseDir/permutations/removeOne.txt" or die "Cannot open removeOne.txt: $!";
while(<$removeOneFH>) { 
    chomp; 
    push @removeOneList,$_;
} 
close $removeOneFH;

print $logfastTreeFH "You are running a Remove One Permutation\n";

print $logfastTreeFH "Starting on RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n";


my $filename = "$baseDir/concatFilterGroup/RemoveOne_$removeOneList[$index]_Group.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree/";


system("$fastaTreeMPPath/FastTreeMP $trueGamma -$myfastTreeModel -log $baseDir/fastTree/removeOneFullReRun/RemoveOne_$removeOneList[$index]/fastTreeLog_$removeOneList[$index].txt $filename > $baseDir/fastTree/removeOneFullReRun/RemoveOne_$removeOneList[$index]/Tree_$removeOneList[$index].nwk");

print $logfastTreeFH "Finished RemoveOne_$removeOneList[$index]_Group.fasta tree creation\n\n";





##########
#Keep One#
##########

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



my $filename = "$baseDir/concatFilterGroup/KeepOne_$keepOneList[$index]_Group.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree/";


system("$fastaTreeMPPath/FastTreeMP $trueGamma -$myfastTreeModel -log $baseDir/fastTree/keepOneFullReRun/KeepOne_$keepOneList[$index]/fastTreeLog_$keepOneList[$index].txt $filename > $baseDir/fastTree/keepOneFullReRun/KeepOne_$keepOneList[$index]/Tree_$keepOneList[$index].nwk");



print $logfastTreeFH "Finished KeepOne_$keepOneList[$index]_Group.fasta tree creation\n\n";


##############
#Remove Group#
##############

} elsif ($ENV{'PoPipe_removeOne'} == 0 && $ENV{'PoPipe_keepOne'} == 0 && $ENV{'poPipe_removeGroup'} == 1) {


print $logfastTreeFH "You are running a Remove Group Permutation\n";

print $logfastTreeFH "Starting on RemoveGroup_FullReRunGroup.fasta tree creation\n";



my $filename = "$baseDir/concatFilterGroup/RemoveGroup_FullReRunGroup.fasta";
my $fastaTreeMPPath = "$baseDir/externalSoftware/fastTree/";


system("$fastaTreeMPPath/FastTreeMP $trueGamma -$myfastTreeModel -log $baseDir/fastTree/removeGroupFullReRun/fastTreeLog_RemoveGroup.txt $filename > $baseDir/fastTree/removeGroupFullReRun/Tree_RemoveGroup.nwk");



print $logfastTreeFH "RemoveGroup_FullReRunGroup.fasta tree creation\n\n";


} else {

die ("You cannot run multiple permutations at once!! Double check the poPipeSettings file.");

}


} else {

die ("You cannot run fastTree and RAxML at the same time. Double check the poPipeSettings file.");

}
