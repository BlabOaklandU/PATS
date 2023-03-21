#!/bin/sh
fullname=$(readlink -f $0)
myname=$(basename $fullname)
if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;
if [ -z "$PoPipe_fastaDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;
if [ -z "$PoPipe_protOrthoDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;

polog "$myname begin";

sname=speciesFilter
polog "$sname begin";
PoPipe_osfDir=$PoPipe_baseDir/$PoPipe_speciesFilterDirName; export PoPipe_osfDir; 
PoPipe_ogDir=$PoPipe_baseDir/$PoPipe_ogDirName; export PoPipe_ogDir; 
rm -rf $PoPipe_ogDir;
PoPipe_orthoDir=$PoPipe_baseDir/$PoPipe_orthologsDirName; export PoPipe_orthoDir;
rm -rf $PoPipe_orthoDir; mkdir $PoPipe_orthoDir; 
ln -s $PoPipe_protOrthoDir/$PoPipe_poGroupsName $PoPipe_orthoDir; 
rm -rf $PoPipe_osfDir; 
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/SpeciesFilter.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=OGFilter
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/OGFilter.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

sname=muscle
polog "$sname begin";
PoPipe_muscleDir=$PoPipe_baseDir/$PoPipe_muscleDirName; export PoPipe_muscleDir; 
rm -rf $PoPipe_muscleDir; 
mkdir $PoPipe_muscleDir; 
if [[ ! -d $PoPipe_muscleDir ]]; then echo "could not create $PWD/$PoPipe_muscleDir"; exit 10; fi
find ./orthoGroupsFasta -name "*.fasta" -printf "%f\n"| \
while read f; do \
std=$PoPipe_stdDir/$sname

if [ "$PoPipe_muscleVersion" == "5.1.0" ]
then
$PoPipe_srcBaseDir/externalSoftware/muscle/$PoPipe_muscleVersion/muscle5.1.linux_intel64 -$PoPipe_muscleAlgorithm $PoPipe_srcBaseDir/orthoGroupsFasta/$f -output $PoPipe_srcBaseDir/muscleFasta/$f 1>>$std.out 2>>$std.err;
elif [ "$PoPipe_muscleVersion" == "3.8.31" ]
then
$PoPipe_srcBaseDir/externalSoftware/muscle/$PoPipe_muscleVersion/muscle3.8.31_i86linux64 -$PoPipe_muscleAlgorithm $PoPipe_srcBaseDir/orthoGroupsFasta/$f -output $PoPipe_srcBaseDir/muscleFasta/$f 1>>$std.out 2>>$std.err;
else
$PoPipe_srcBaseDir/externalSoftware/muscle/$PoPipe_muscleVersion/muscle5.1.linux_intel64 -$PoPipe_muscleAlgorithm $PoPipe_srcBaseDir/orthoGroupsFasta/$f -output $PoPipe_srcBaseDir/muscleFasta/$f 1>>$std.out 2>>$std.err;
fi
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
done
polog "$sname end";

sname=gapFilter
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/GapFilter.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

polog "$myname end";

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then return 0; else exit 0; fi;