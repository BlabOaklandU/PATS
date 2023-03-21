#!/bin/sh
fullname=$(readlink -f $0)
myname=$(basename $fullname)

if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;


if [ "$molecularClockAnalysis" -eq 1 ] && [ "$molecularClockVarification" == 'YES' ] && [ "$moleClocksOnly" -eq 1 ] && [ "$moleClocksBasicPats" -eq 0 ] && [ "$moleClockPermPats" -eq 0 ]; then

polog " ";
polog "======================================================================= ";
polog "Molecular Clock Analysis Only";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";

#molClockPipeFull.sh
step=molClockPipeFull.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

polog "$myname end";

return 0

else


# Three different itterations of the pipeline:
# 1) standard run from start to end
# 2) tree only with permutations
# 3) full analysis with permutations

# Option #1: This will ALWAYS rerun the FULL analysis. This is recommended for first run



if [ "$PoPipe_rerunTreesOnly" -eq 0 ] && [ "$PoPipe_rerunFullAnalysis" -eq 0 ]; then
now=$(date)

polog " ";
polog "======================================================================= ";
polog "Full Run - Standard";
polog "No Permutations";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";


#Run FastaStandardizer
step=fastaStandardizer.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeProteinOrtho.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeSpeciesFilter_Muscle_gapFilter.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeTreeCreation.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeFileCleanUp.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

polog "$myname end";

return 0

# Option #2: This will rerun ONLY the tree creation. This is recommended for permutation additions

elif [ "$PoPipe_rerunTreesOnly" -eq 1 ] && [ "$PoPipe_rerunFullAnalysis" -eq 0 ]; then
now=$(date)

if [ "$PoPipe_removeOne" -eq 1 ] && [ "$poPipe_removeGroup" -eq 0 ] && [ "$PoPipe_keepOne" -eq 0 ]; then rrperm="Permutation: Remove One";

elif [ "$PoPipe_removeOne" -eq 0 ] && [ "$poPipe_removeGroup" -eq 1 ] && [ "$PoPipe_keepOne" -eq 0 ]; then rrperm="Permutation: Remove Group";

elif [ "$PoPipe_removeOne" -eq 0 ] && [ "$poPipe_removeGroup" -eq 0 ] && [ "$PoPipe_keepOne" -eq 1 ]; then rrperm="Permutation: Keep One";

else rrperm="Permutation: Standard Run"
fi;

polog " ";
polog "======================================================================= ";
polog "Re-run Tree Creation";
polog "$rrperm";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";

step=fastaStandardizer.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeTreeCreation.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeFileCleanUp.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

polog "$myname end";
return 0


# Option #3: This will rerun the full analysis allowing for permutation additions
# This option has three subsets; 1) KeepOne, 2) RemoveOne, 3) RemoveGroup

elif [ "$PoPipe_rerunTreesOnly" -eq 0 ] && [ "$PoPipe_rerunFullAnalysis" -eq 1 ]; then

# Option #3 Subset 1: KeepOne Permutation for a full rerun.

if [ "$PoPipe_removeOne" -eq 0 ] && [ "$poPipe_removeGroup" -eq 0 ] && [ "$PoPipe_keepOne" -eq 1 ] && [ "$PoPipe_rerunFullAnalysis" -eq 1 ]; then
now=$(date)

polog " ";
polog "======================================================================= ";
polog "Multiple Full Run - Permutations";
polog "Keep One Permutation Selected";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";

step=fastaStandardizer.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeKeepOneInit.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeProteinOrtho.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeSpeciesFilter_Muscle_gapFilter.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeTreeCreationFullRun.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeKeepOneFileCleanUp.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

polog "$myname end";
return 0


# Option #3 Subset 2: RemoveOne Permutation for a full rerun.

elif [ "$PoPipe_removeOne" -eq 1 ] && [ "$poPipe_removeGroup" -eq 0 ] && [ "$PoPipe_keepOne" -eq 0 ] && [ "$PoPipe_rerunFullAnalysis" -eq 1 ]; then
now=$(date)

polog " ";
polog "======================================================================= ";
polog "Multiple Full Run - Permutations";
polog "Remove One Permutation Selected";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";

#Run FastaStandardizer
step=fastaStandardizer.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeRemoveOneInit.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeProteinOrtho.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeSpeciesFilter_Muscle_gapFilter.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeTreeCreationFullRun.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeRemoveOneFileCleanUp.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

polog "$myname end";
return 0


# Option #3 Subset 2: RemoveGroup Permutation for a full rerun.

elif [ "$PoPipe_removeOne" -eq 0 ] && [ "$poPipe_removeGroup" -eq 1 ] && [ "$PoPipe_keepOne" -eq 0 ] && [ "$PoPipe_rerunFullAnalysis" -eq 1 ]; then
now=$(date)

polog " ";
polog "======================================================================= ";
polog "Multiple Full Run - Permutations";
polog "Remove Group Permutation Selected";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";

#Run FastaStandardizer
step=fastaStandardizer.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeRemoveGroupInit.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeProteinOrtho.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeSpeciesFilter_Muscle_gapFilter.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeTreeCreationFullRun.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

step=poPipeRemoveGroupFileCleanUp.sh
polog "$step called";
source $step
res=$?; if [ $res -ne 0 ]; then polog "$step error $res"; exit $res; fi; 
polog "$step returned";

polog "$myname end";


else

polog "There seems to be an error in the poPipeSettings.sh file";

fi;


else

polog "There seems to be an error in the poPipeSettings.sh file";

fi;


fi;