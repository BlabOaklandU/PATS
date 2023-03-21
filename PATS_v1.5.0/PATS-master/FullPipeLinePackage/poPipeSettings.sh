#!/bin/sh
export PoPipe_baseDir=$PWD
export PoPipe_stdDir=$PWD/std
export MolClock_stdDir=$PWD/molecularClocks/std
if [[ -d $PoPipe_stdDir ]]; then if [[ $(ls -A $PoPipe_stdDir) ]]; then echo "PoPipe_stdDir already exists: $PoPipe_stdDir"; exit 10; fi; fi;
mkdir $PoPipe_stdDir &> /dev/null;
if [[ ! -d $PoPipe_stdDir ]]; then echo "could not create $PoPipe_stdDir"; exit 10; fi;

if [[ -d $MolClock_stdDir ]]; then if [[ $(ls -A $MolClock_stdDir) ]]; then echo "MolClock_stdDir already exists: $PoPipe_stdDir"; exit 10; fi; fi;
mkdir $MolClock_stdDir &> /dev/null;
if [[ ! -d $MolClock_stdDir ]]; then echo "could not create $MolClock_stdDir"; exit 10; fi;

###########################################################################################
#                       Set Work Directory                                                #
###########################################################################################
export PoPipe_srcBaseDir=/scratch/projects/blab/cpowell/PATS_v1.3.0_MolecularClocks/PATS-master/FullPipeLinePackage
###########################################################################################


######################################
#       General Settings             # 
######################################
export PoPipe_multiInfCount=100
export PoPipe_speciesFilterCutoff=75
export PoPipe_gapFilterCutoff=25
export PoPipe_proteinOrthoConn=0.1
######################################


#########################################################################
#                  Taxon Sampling Settings - Do NOT edit                #
#########################################################################
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
export PoPipe_cpusReserved=0
export PoPipe_fastaSrcDir=$PoPipe_srcBaseDir/fasta
export PoPipe_fastaDir=$PoPipe_baseDir/$PoPipe_fastaDirName
export PoPipe_protOrthoDir=$PoPipe_baseDir/$PoPipe_protOrthoDirName
export PoPipe_settingsFile=$PoPipe_srcBaseDir/poPipeSettings.sh

#########################################################################
#                  Molecular Clock Analysis Settings - Do NOT edit      #
#########################################################################
export MolecularClockBaseDir=$PoPipe_baseDir/molecularClocks
export MolecularClockScriptDir=$PoPipe_baseDir/molecularClocks/molClockScripts
export nullDir=$MolecularClockBaseDir/nullCheck
export unsplitMainTreeDir=$MolecularClockBaseDir/unsplitMainTree
#export phyFileDir=$MolecularClockBaseDir/phyFiles
export splitTreeDir=$MolecularClockBaseDir/splitTrees
export splitTreeADir=$MolecularClockBaseDir/splitTrees/treeA
export splitTreeBDir=$MolecularClockBaseDir/splitTrees/treeB
export calibrationTreeDir=$MolecularClockBaseDir/calibrationTrees
export calibrationATreeDir=$MolecularClockBaseDir/calibrationTrees/treeA
export calibrationBTreeDir=$MolecularClockBaseDir/calibrationTrees/treeB
export indvGenSequencesDir=$MolecularClockBaseDir/indvGeneSequences
export concatGenSequencesDir=$MolecularClockBaseDir/groupConcatGenes
export splitGenesDir=$MolecularClockBaseDir/splitGenes
export megaOutputDir=$MolecularClockBaseDir/megaOutput
export rootFastaFile=OG
########################################################################

##################
# User Settings: #
##################

#ProteinOrtho and MUSCLE Version Choice
#ProteinOrtho versions 6.0.14 (GALAXY/STABLE); 6.1.5 and 6.1.7 (NEWEST/EXPERIMENTAL)
#Type version number to set the correct version
export PoPipe_proteinOrthoVersion="6.1.7"

#MUSCLE versions 3.8.31 and 5.1.0
#Type version number to set the correct version
#Default: version 5.1.0
export PoPipe_muscleVersion="5.1.0"

#MUSCLE alignment algorithm
#Type either align (ppp algorithm) or super5 (Super5 algorithm)
#align - smaller datasets ~<1k sequences; super5 - larger datasers ~>1k sequences
export PoPipe_muscleAlgorithm=align

#Create archive of completed pipe? (compressed .7z)
#This will compress all populated files via the pipeline
export PoPipe_createArchive=0


#Delete populated folders after archive is created. 
#ONLY USE IF YOU ARE ARCHIVING!!!!
export PoPipe_folderCleanUp=0


#Choice of BLAST+ or Diamond
#Diamond: Translated DNA and Proteins
#BLAST+: Nucleotides, Proteins, Translated DNA
#PoPipe_querySoftware: Diamond=0; autoblast=1; Blastn=2; Blastp=3; Blastx=4
export PoPipe_querySoftware=0


#Tree Selection Protocol
#Only select one tree creation protocol (1 = TRUE; 0 = FALSE)
export PoPipe_fastTree=1
export PoPipe_RAxMLTree=0


#Run only the tree creation (can be combined with permutations)
#This will ONLY rerun the trees based off the orthologs previously created from a full run
export PoPipe_rerunTreesOnly=0


#Run full analysis from start to finish (can be combined with permutations)
#This will permutate then generate the orthologs. 
export PoPipe_rerunFullAnalysis=0


#Permutation Settings (1 = TRUE; 0 = FALSE)
#Only one permutation allowed per run. No permutations = standard run
export PoPipe_keepOne=0
export PoPipe_removeOne=0
export poPipe_removeGroup=0

####################
#Fast Tree Settings#
####################
#FastTree model options: (1) lg (2) wag (3) gtr
export fastTreeModel=lg

#Use gamma likelihood (1 = TRUE; 0 = FALSE)
export gammaChoice=1

################
#RAxML Settings#
################
#RAxML model options: (1) lg (2) wag (3) gtr
#export $RAxMLTreeModel=lg

#Use gamma likelihood (1 = TRUE; 0 = FALSE)
#export $gammaChoice=1


################################
#ProteinOrtho parallel Settings#
################################
#ONLY USE IF YOU UNDERSTAND (still in testing)
#Breaks up step 2 of PO and pushes each chunk to a different node
#Node running pipeline does NOT run PO step 2 chunks if parallel is used
#N/M; N = step, M = Number of nodes being used
#ProteinOrtho Large Job Split(1 = True; 0 = False)
export POLargeJobSplit=0
#This is for safety of the pipeline
#ProteinOrtho Large Job Split Varification (1) = YES; 
export POLargeJobSplitVarification=

#Set number of avalibale nodes (either '2', '5', '10', '20'; default ='1' for non-parallel)
export PONumNodes=5

#SLURM pass off settings (should match SBATCH)
export cpusPerTask=40

#Should match SLURM --mem settings in MB
export memUsage100=102400






#####################################
# Molecular Clock Analysis Settings #
#####################################

#Add a moleucular clock anlysis to end of Taxon sampling analysis? (1 = TRUE; 0 = FALSE)
export molecularClockAnalysis=1
#Molecular clock Varification Type: YES (all caps)
export molecularClockVarification=YES


#Run just Molecular clock Analysis without PATS
#If running this setting, plaese make sure you haev populated the proper directories
#You will need to do the following:
#Add full tree.nwk to unsplitMainTree dir
#Add full concatenated ortholog file (concatFilterGroup from PATS) to groupConcatGenes dir
export moleClocksOnly=1

#Run molecualr clock analysis after basic run of PATS (no permutations or file cleanup)
#Files from PATS have to be populated for this setting to work
export moleClocksBasicPats=0

#Run molecular clock analyses after PATS permutations (permutations + file archive on)
#This setting will utilize the .zips created from each permuation run
export moleClockPermPats=0

#Pick the data type for the gene files
#Utilize concatenated genes = 1
#Utilize individual genes = 2
export geneFileType=1


#Root Output gene folder name
#Example: 123Genes
export rootOutputDirName=123Genes

#Number of genes to be skipped
#Enther these one by one with a space between (ex: 1 56 99)
export numberOfGenesSkippedMC="1 56 99"

#nodesToBeDeleted
#Enter these one by one with a space between (ex: 3 14)
export nodesToBeDeleted="3 14"

#Tree A output file name
#Example: Dibamidae_Hele_Dicro_A_Tree.nwk
export outputFileNameTreeA=">Dibamidae_Hele_Dicro_A_Tree.nwk"

#Tree B output file name
#Example: Dibamidae_Hele_Dicro_B_Tree.nwk
export outputFileNameTreeB=">Dibamidae_Hele_Dicro_B_Tree.nwk"

#Species to split tree on (splitter species)
#Example: Dibamidae
export splitterSpecies="Dibamidae"

#Overlap Species 1
#Example: Heleophryn
export overlapSpecies1="Heleophryn"

#Overlap Species 2
#Example: Dicrogloss
export overlapSpecies2="Dicrogloss"

#Outgroup speices
#Example: Polyodonti
export outgroupSpecies="Polyodonti"

#First species in B Tree
#Example: Cryptobran
export firstSpeciesTreeB="Cryptobran"


#MCMCTree Settings
#Outfile names for both Tree A and Tree B
export outFileTreeA=aaa
export outFileTreeA=bbb

#Ingroup to tip disatnce A Tree
export ingroupTipDistanceTreeA=0.3612

#Ingroup to tip disatnce B Tree
export ingroupTipDistanceTreeB=0.2940


function polog 
{
echo -e $(date +%Y.%m.%d.%H.%M.%S) '\t' $1 >> $PoPipe_baseDir/log.txt;
}
export -f polog
