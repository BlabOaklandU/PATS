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


if [ -n "$(ls -A $MolecularClockBaseDir/splitTreeDirectory/treeA/ 2>/dev/null)" ] || [ -n "$(ls -A $MolecularClockBaseDir/splitTreeDirectory/treeB/ 2>/dev/null)" ]; then

sname=phyToMeg
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/phyToMeg.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=gammaCalc
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/gammaCalc.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";


else

sname=TreeSplitterA
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/treeSplitterA.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=TreeSplitterB
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/treeSplitterB.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=nullChecker
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/nullChecker.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=geneSplitter
polog "$sname begin";
std=$MolClock_stdDir/$sname
perl $MolecularClockBaseDir/molClockScripts/geneSplitter.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

fi;