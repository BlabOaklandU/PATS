#!/bin/sh
export PoPipe_baseDir=$PWD
export PoPipe_stdDir=$PWD/std
if [[ -d $PoPipe_stdDir ]]; then if [[ $(ls -A $PoPipe_stdDir) ]]; then echo "PoPipe_stdDir already exists: $PoPipe_stdDir"; exit 10; fi; fi
mkdir $PoPipe_stdDir &> /dev/null;
if [[ ! -d $PoPipe_stdDir ]]; then echo "could not create $PoPipe_stdDir"; exit 10; fi
export PoPipe_gapFilterDirName=gapFilter
export PoPipe_groupName=groups.fasta
export PoPipe_poGroupsName=myproject.proteinortho.tsv
export PoPipe_fastaDirName=fasta
export PoPipe_speciesFilterDirName=speciesFilter
export PoPipe_ogDirName=og
export PoPipe_orthologsDirName=orthologs
export PoPipe_raxmlDirName=treeRAxML
export PoPipe_fasttreeDirName=fasttree
export PoPipe_muscleDirName=muscleFasta
export PoPipe_protOrthoDirName=poWork
export PoPipe_srcBaseDir=/mnt/gs18/scratch/users/u5g6jsxd/PATS_v0.91_28June2020/PATS-master/FullPipeLinePackage
export PoPipe_fastaSrcDir=$PoPipe_srcBaseDir/fasta
export PoPipe_multiInfCount=100
export PoPipe_speciesFilterCutoff=60
export PoPipe_gapFilterCutoff=25
export PoPipe_proteinOrthoConn=0.1
export PoPipe_cpusReserved=0
export PoPipe_fastaDir=$PoPipe_baseDir/$PoPipe_fastaDirName
export PoPipe_protOrthoDir=$PoPipe_baseDir/$PoPipe_protOrthoDirName
export PoPipe_settingsFile=$PoPipe_srcBaseDir/poPipeSettings.sh

#Threads to Use Settings
export OMP_NUM_THREADS=$(( $(cat /proc/cpuinfo |grep -c "^processor") - 2 ));
export PoPipe_threadsToUse=$(( $(cat /proc/cpuinfo |grep -c "^processor") - 2 ));

##################
# User Settings: #
##################

#Create archive of completed pipe? (compressed .7z)
#This will compress all populated files via the pipeline
export PoPipe_createArchive=1


#Delete populated folders after archive is created. 
#ONLY USE IF YOU ARE ARCHIVING!!!!
export PoPipe_folderCleanUp=1


#Tree Selection Protocol
#Only select one tree creation protocol (1 = TRUE; 0 = FALSE)
export PoPipe_fastTree=1
export PoPipe_RAxMLTree=0


#Run only the tree creation (can be combined with permutations)
#This will ONLY rerun the trees based off the orthologs previously created from a full run
export PoPipe_rerunTreesOnly=0


#Run full analysis from start to finish (can be combined with permutations)
#This will permutate then generate the orthologs. 
export PoPipe_rerunFullAnalysis=1


#Permutation Settings (1 = TRUE; 0 = FALSE)
#Only one permutation allowed per run. No permutations = standard run
export PoPipe_keepOne=0
export PoPipe_removeOne=1
export poPipe_removeGroup=0




function polog 
{
echo -e $(date +%Y.%m.%d.%H.%M.%S) '\t' $1 >> $PoPipe_baseDir/log.txt;
}
export -f polog
