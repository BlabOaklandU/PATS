#!/bin/sh
fullname=$(readlink -f $0)
myname=$(basename $fullname)

if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;
if [ -z "$MolecularClockBaseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;


sname=MoleClock_Master
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/MoleClock_Master.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=MoleClock_CalibrationGenerator
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/molClockCalibrationGenerator.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=MoleClock_Concat
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/molClockConcatGene.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";




#Run Molecular Clock Analysis Only
if [ "$molecularClockAnalysis" -eq 1 ] && [ "$molecularClockVarification" == 'YES' ] && [ "$moleClocksOnly" -eq 1 ] && [ "$moleClocksBasicPats" -eq 0 ] && [ "$moleClockPermPats" -eq 0 ]; then


if [ "$geneFileType" -eq 1 ]; then

sname=molClockConcatGene
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/molClockConcatGene.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";


elif [ "$geneFileType" -eq 2 ]; then

sname=molClockSplitGene
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/molClockSplitGene.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";


else
polog "Unsure of the type of data to use";
polog "please make sure geneFileType is set";
fi;




#Run Molecular Clock after PATS standard run
elif [ "$molecularClockAnalysis" -eq 1 ] && [ "$molecularClockVarification" == 'YES' ] && [ "$moleClocksOnly" -eq 0 ] && [ "$moleClocksBasicPats" -eq 1 ] && [ "$moleClockPermPats" -eq 0 ]; then

polog " ";
polog "======================================================================= ";
polog "Molecular Clock Analysis";
polog "Based on PATS standard run output";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";



#Run Molecular Clock after PATS permutation run
elif [ "$molecularClockAnalysis" -eq 1 ] && [ "$molecularClockVarification" == 'YES' ] && [ "$moleClocksOnly" -eq 0 ] && [ "$moleClocksBasicPats" -eq 0 ] && [ "$moleClockPermPats" -eq 1 ]; then

polog " ";
polog "======================================================================= ";
polog "Molecular Clock Analyses";
polog "Based on PATS permutation .zip files";
polog "Timestamp: $now";
polog "======================================================================= ";
polog " ";
polog "$myname begin";


else

polog "There seems to be an error in the poPipeSettings.sh file";
polog "Check to make sure molecular clock settings are correct";

fi;








#if [ -n "$(ls -A $MolecularClockBaseDir/splitTreeDirectory/treeA/ 2>/dev/null)" ] || [ -n #3"$(ls -A $MolecularClockBaseDir/splitTreeDirectory/treeB/ 2>/dev/null)" ]; then

#sname=phyToMeg
#polog "$sname begin";
#std=$MolClock_stdDir/$sname
#perl $MolecularClockBaseDir/molClockScripts/phyToMeg.pl 1>$std.out 2>$std.err;
#res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
#polog "$sname end";

#sname=gammaCalc
#polog "$sname begin";
#std=$MolClock_stdDir/$sname
#perl $MolecularClockBaseDir/molClockScripts/gammaCalc.pl 1>$std.out 2>$std.err;
#res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
#polog "$sname end";


#else

#sname=TreeSplitterA
##polog "$sname begin";
#std=$MolClock_stdDir/$sname
###perl $MolecularClockBaseDir/molClockScripts/treeSplitterA.pl 1>$std.out 2>$std.err;
#res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
#polog "$sname end";

#sname=TreeSplitterB
#polog "$sname begin";
#std=$MolClock_stdDir/$sname
#perl $MolecularClockBaseDir/molClockScripts/treeSplitterB.pl 1>$std.out 2>$std.err;
#res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
#polog "$sname end";

#sname=nullChecker
#polog "$sname begin";
###std=$MolClock_stdDir/$sname
#perl $MolecularClockBaseDir/molClockScripts/nullChecker.pl 1>$std.out 2>$std.err;
#res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
#polog "$sname end";

#sname=geneSplitter
#polog "$sname begin";
#std=$MolClock_stdDir/$sname
#perl $MolecularClockBaseDir/molClockScripts/geneSplitter.pl 1>$std.out 2>$std.err;
#res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
#polog "$sname end";

#fi;