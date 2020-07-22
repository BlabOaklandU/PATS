#!/bin/sh
sname=Concat-TreeCreate-Archive
polog "$sname begin";
if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;

polog "$myname begin";

sname=concatFilter
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/ConcatFilter.pl 2>$std.err 1>$std.out;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=TreeCreation
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/CreateTrees.pl 2>$std.err 1>$std.out;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=Archive
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/Archive.pl 2>$std.err 1>$std.out;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=ArchiveCleanUp
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/ArchiveCleanUp.pl 2>$std.err 1>$std.out;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

polog "$myname end";

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then return 0; else exit 0; fi