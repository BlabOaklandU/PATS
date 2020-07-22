#!/bin/sh

fullname=$(readlink -f $0)
myname=$(basename $fullname)
if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;
if [ -z "$PoPipe_fastaDir" ]; then echo "$myname PoPipe_fastaDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;
if [ -z "$PoPipe_protOrthoDir" ]; then echo "$myname PoPipe_protOrthoDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;

polog "$myname begin";
mkdir $PoPipe_protOrthoDir;
if [[ ! -d $PoPipe_protOrthoDir ]]; then echo "could not create $PWD/$PoPipe_protOrthoDir"; exit 10; fi
pushd $PoPipe_protOrthoDir >/dev/null

sname=indices
polog "$sname begin";
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/Blast/bin -p=autoblast -step=1 -verbose=1 -checkfasta -keep -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>$std.out 2>$std.err ;
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi
polog "$sname end";

sname=blaGraph
polog "$sname begin";
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/Blast/bin -p=autoblast -step=2 -verbose=1 -keep -cleanblast -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta -clean 1>$std.out 2>$std.err;
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi
polog "$sname end";

sname=cluster
polog "$sname begin";
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/proteinortho6.pl -step=3 -debug -verbose=1 -keep -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta -clean 1>$std.out 2>$std.err;
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi
polog "$sname end";

popd >/dev/null

polog "$myname end";
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then return 0; else exit 0; fi