#!/bin/sh
fullname=$(readlink -f $0)
myname=$(basename $fullname)

if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;



#Step One: Remove ALL '*' and '-' from .fasta files.

#Change '*', ' ', and '-' to '_' within the .fasta file names
for file in fasta/*.fasta; do mv "$file" "${file//-/_}"; done
for file in fasta/*.fasta; do mv "$file" "${file// /_}"; done
for file in fasta/*.fasta; do mv "$file" "${file//\*/_}"; done

#Change '*' and '-' to '_' from within .fasta files (acession line and sequences)
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
perl -i.bak -p -e 's/\,/_/g;' fasta/*.fasta
perl -i.bak -p -e 's/\=/_/g;' fasta/*.fasta



#Step Two: Standardize all accession line lenghts within fasta files

perl $PoPipe_srcBaseDir/scriptsPerl/fastaStandardizer.pl;



#Step Three: Recheck File Sizes and cleanup.

perl $PoPipe_srcBaseDir/scriptsPerl/fastaStandardizerCheck.pl;

