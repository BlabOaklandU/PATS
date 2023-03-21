#!/bin/sh
fullname=$(readlink -f $0)
myname=$(basename $fullname)

if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;

if [ -e "$PoPipe_srcBaseDir/DataStandardizerLogs/goldNum.txt" ]
then
#Run fastaInitCheck
sname=fastaInitCheck
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/fastaInitCheck.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";

polog "Skipping FastaStandardizer, pipeline is ready to run";

else

#Step One: Remove ALL '*' and '-' from .fasta files.

#Change '*', ' ', and '-' to '_' within the .fasta file names
for file in fasta/*.fasta; do mv "$file" "${file//-/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file// /_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\*/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\'/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\)/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\(/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\]/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\[/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\}/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\{/_}"; done

#Change special characters to '_' from within .fasta files (accession line and sequences)
perl -i.bak -p -e 's/\*/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\-/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\+/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\)/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\(/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\[/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\]/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\{/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\}/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\=/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\$/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\%/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\^/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\&/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\#/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\@/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\!/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\?/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\,/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\|/_/g;' fasta/*.fasta
perl -i.bak -p -e "s/\'/_/g;" fasta/*.fasta
perl -i.bak -p -e 's/ //g;' fasta/*.fasta

#Change trailing special characters
perl -i.bak -p -e 's/\_+$//g;' fasta/*.fasta
perl -i.bak -p -e 's/\.+$//g;' fasta/*.fasta


#Remove duplicate special characters
perl -i.bak -p -e 's/([^a-zA-Z0-9])\1+/$1/g;' fasta/*.fasta 


#Step Two: Standardize all accession line lenghts within fasta files
sname=fastaStandardizer
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/fastaStandardizer.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";


#Step Three: Recheck File Sizes and cleanup.
sname=fastaStandardizerChecker
polog "$sname begin";
std=$PoPipe_stdDir/$sname
perl $PoPipe_srcBaseDir/scriptsPerl/fastaStandardizerCheck.pl 1>$std.out 2>$std.err;
res=$?; if [ $res -ne 0 ]; then polog "$sname error $res"; exit $res; fi; 
polog "$sname end";
fi

polog "$myname end";
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then return 0; else exit 0; fi;