#!/bin/sh
fullname=$(readlink -f $0)
myname=$(basename $fullname)

if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;

polog "$myname begin";

sname=KeepOneRerunFinal
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/KeepOneRerunFin.pl 2>$std.err 1>$std.out;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=FileCleanUp
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/FileCleanUp.pl 2>$std.err 1>$std.out;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

polog "$myname end";
return 0