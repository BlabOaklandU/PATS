#!/usr/bin/perl
  
use strict;
use warnings;
use File::Copy;
use File::Path;
use File::Find;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);


my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $fastaDir = "$baseDir/fasta";


#Create directory
my $directoryfBU = "$fastaDir/fastaOriginal";
unless(-d $directoryfBU)
{
mkdir $directoryfBU;
}

my $directoryFileStand = "$baseDir/DataStandardizerLogs";
unless(-d $directoryFileStand)
{
mkdir $directoryFileStand;
}

my $directoryFastaEdit = "$fastaDir/fastaEdit";
unless(-d $directoryFastaEdit)
{
mkdir $directoryFastaEdit;
}

my $directoryFastaEdit2 = "$fastaDir/fastaEdit2";
unless(-d $directoryFastaEdit2)
{
mkdir $directoryFastaEdit2;
}

my $directoryFastaEdit3 = "$fastaDir/fastaChangedBackUp";
unless(-d $directoryFastaEdit3)
{
mkdir $directoryFastaEdit3;
}


my $fastaBUDir = "$fastaDir/fastaOriginal";
my $fastaEditDir = "$fastaDir/fastaEdit";
my $fastaEditDir2 = "$fastaDir/fastaEdit2";
my $fastaEditDir3 = "$fastaDir/fastaChangedBackUp";


#generate logs
open my $logfastaFH, '>', "$baseDir/DataStandardizerLogs/fastaStandardLog.txt" or die "Cannot open fastaStandardLog.txt: $!";
open my $logfastaCountFH, '>', "$baseDir/DataStandardizerLogs/fastaCountLog.txt" or die "Cannot open fastaCountLog.txt: $!";
open my $logfastaDataFH, '>', "$baseDir/DataStandardizerLogs/fastaDataLog.txt" or die "Cannot open fastaDataLog.txt: $!";
open my $logfastaDebugFH, '>', "$baseDir/DataStandardizerLogs/fastaDebugLog.txt" or die "Cannot open fastaDebugLog.txt: $!";
open my $logcharCountFH, '>', "$baseDir/DataStandardizerLogs/dataTypeLog.txt" or die "Cannot open dataTypeLog.txt: $!";

opendir(DIR,"$fastaDir") or die "Couldn't open $fastaDir\n";

my @fastaFiles = grep(/\.fasta$/,readdir(DIR));
my $numFasta = scalar(@fastaFiles);




foreach my $filename (@fastaFiles) {
print $logfastaFH "$filename\n";
}

closedir(DIR);


my $itter = 0;
my @fastaArrayTemp;
my $fileStart = '>';
my @fastaCountArray;
my @fastaFinalCountArray;
my @dataTypeArray;
my $dataInputType;

#Nucleotides
my $adenineChar = "a";
my $cytosineChar = "c";
my $thymineChar = "t";
my $guanineChar = "g";
my $uracilChar = "u";
my $lowComplexMaskChar = "x";
my $ambiguousChar = "n";
my $adenineCharDataCounter = 0;
my $cytosineCharDataCounter = 0;
my $thymineCharDataCounter = 0;
my $guanineCharDataCounter = 0;
my $uracilCharDataCounter = 0;
my $lowComplexMaskCharDataCounter = 0;
my $ambiguousCharDataCounter = 0;
my $nucleotideSum;

#Amino Acids
my $alaChar = "A";
my $argChar = "R";
my $asnChar = "N";
my $aspChar = "D";
my $cysChar = "C";
my $glnChar = "Q";
my $gluChar = "E";
my $glyChar = "G";
my $hisChar = "H";
my $ileChar = "I";
my $leuChar = "L";
my $lysChar = "K";
my $metChar = "M";
my $pheChar = "F";
my $proChar = "P";
my $pylChar = "O";
my $serChar = "S";
my $secChar = "U";
my $thrChar = "T";
my $trpChar = "W";
my $tyrChar = "Y";
my $valChar = "V";
my $asxChar = "B";
my $glxChar = "Z";
my $xaaChar = "X";
my $xleChar = "J";
my $alaCharDataCounter = 0;
my $argCharDataCounter = 0;
my $asnCharDataCounter = 0;
my $aspCharDataCounter = 0;
my $cysCharDataCounter = 0;
my $glnCharDataCounter = 0;
my $gluCharDataCounter = 0;
my $glyCharDataCounter = 0;
my $hisCharDataCounter = 0;
my $ileCharDataCounter = 0;
my $leuCharDataCounter = 0;
my $lysCharDataCounter = 0;
my $metCharDataCounter = 0;
my $pheCharDataCounter = 0;
my $proCharDataCounter = 0;
my $pylCharDataCounter = 0;
my $serCharDataCounter = 0;
my $secCharDataCounter = 0;
my $thrCharDataCounter = 0;
my $trpCharDataCounter = 0;
my $tyrCharDataCounter = 0;
my $valCharDataCounter = 0;
my $asxCharDataCounter = 0;
my $glxCharDataCounter = 0;
my $xaaCharDataCounter = 0;
my $xleCharDataCounter = 0;
my $aminoAcidSum;

foreach (1 .. $numFasta) {

open my $fastaFH, '<', "$fastaDir/$fastaFiles[$itter]" or die "Cannot open fasta file: $fastaFiles[$itter]: $!";
while(<$fastaFH>) {
push @fastaArrayTemp, $_;
}

close $fastaFH;

my @characterCountArray;

@characterCountArray = @fastaArrayTemp;

foreach my $charSearch (@characterCountArray) {
if ($charSearch =~ /^$fileStart/) {
#do nothing
} else {
#Nucleotides
$adenineCharDataCounter++ while($charSearch =~ m/$adenineChar/ig);
$cytosineCharDataCounter++ while($charSearch =~ m/$cytosineChar/ig);
$thymineCharDataCounter++ while($charSearch =~ m/$thymineChar/ig);
$guanineCharDataCounter++ while($charSearch =~ m/$guanineChar/ig);
$uracilCharDataCounter++ while($charSearch =~ m/$uracilChar/ig);
$lowComplexMaskCharDataCounter++ while($charSearch =~ m/$lowComplexMaskChar/ig);
$ambiguousCharDataCounter++ while($charSearch =~ m/$ambiguousChar/ig);

#Amino Acids
$alaCharDataCounter++ while($charSearch =~ m/$alaChar/ig);
$argCharDataCounter++ while($charSearch =~ m/$argChar/ig);
$asnCharDataCounter++ while($charSearch =~ m/$asnChar/ig);
$aspCharDataCounter++ while($charSearch =~ m/$aspChar/ig);
$cysCharDataCounter++ while($charSearch =~ m/$cysChar/ig);
$glnCharDataCounter++ while($charSearch =~ m/$glnChar/ig);
$gluCharDataCounter++ while($charSearch =~ m/$gluChar/ig);
$glyCharDataCounter++ while($charSearch =~ m/$glyChar/ig);
$hisCharDataCounter++ while($charSearch =~ m/$hisChar/ig);
$ileCharDataCounter++ while($charSearch =~ m/$ileChar/ig);
$leuCharDataCounter++ while($charSearch =~ m/$leuChar/ig);
$lysCharDataCounter++ while($charSearch =~ m/$lysChar/ig);
$metCharDataCounter++ while($charSearch =~ m/$metChar/ig);
$pheCharDataCounter++ while($charSearch =~ m/$pheChar/ig);
$proCharDataCounter++ while($charSearch =~ m/$proChar/ig);
$pylCharDataCounter++ while($charSearch =~ m/$pylChar/ig);
$serCharDataCounter++ while($charSearch =~ m/$serChar/ig);
$secCharDataCounter++ while($charSearch =~ m/$secChar/ig);
$thrCharDataCounter++ while($charSearch =~ m/$thrChar/ig);
$trpCharDataCounter++ while($charSearch =~ m/$trpChar/ig);
$tyrCharDataCounter++ while($charSearch =~ m/$tyrChar/ig);
$valCharDataCounter++ while($charSearch =~ m/$valChar/ig);
$asxCharDataCounter++ while($charSearch =~ m/$asxChar/ig);
$glxCharDataCounter++ while($charSearch =~ m/$glxChar/ig);
$xaaCharDataCounter++ while($charSearch =~ m/$xaaChar/ig);
$xleCharDataCounter++ while($charSearch =~ m/$xleChar/ig);

}
}

$nucleotideSumDNA = $adenineCharDataCounter + $cytosineCharDataCounter + $thymineCharDataCounter + $guanineCharDataCounter + $ambiguousCharDataCounter + $lowComplexMaskCharDataCounter;

$nucleotideSumRNA = $adenineCharDataCounter + $cytosineCharDataCounter + $guanineCharDataCounter + $uracilCharDataCounter + $ambiguousCharDataCounter + $lowComplexMaskCharDataCounter;

$aminoAcidSum = $alaCharDataCounter + $argCharDataCounter + $asnCharDataCounter + $aspCharDataCounter + $cysCharDataCounter + $glnCharDataCounter + $gluCharDataCounter + $glyCharDataCounter + $hisCharDataCounter + $ileCharDataCounter + $leuCharDataCounter + $lysCharDataCounter + $metCharDataCounter + $pheCharDataCounter + $proCharDataCounter + $pylCharDataCounter + $serCharDataCounter + $secCharDataCounter + $thrCharDataCounter + $trpCharDataCounter + $tyrCharDataCounter + $valCharDataCounter + $asxCharDataCounter + $glxCharDataCounter + $xaaCharDataCounter + $xleCharDataCounter;

undef @characterCountArray;

if (-f "$fastaBUDir/$fastaFiles[$itter].bak")
{
print $logfastaFH "The file: $fastaFiles[$itter].bak already exists in this directory\n";
}else{
fmove("$fastaDir/$fastaFiles[$itter].bak","$fastaBUDir/$fastaFiles[$itter].bak");
}


open my $fastaChangeFH, '>', "$fastaEditDir/$fastaFiles[$itter]" or die "Cannot open $fastaFiles[$itter]: $!";

foreach my $lineSearch (@fastaArrayTemp) {
if ($lineSearch =~ /^$fileStart/) {
for ($lineSearch){
s/$fileStart/$fileStart$fastaFiles[$itter]|*|/g;
}
} else {
}
}

foreach my $lineSearchDupe (@fastaArrayTemp) {
if ($lineSearchDupe =~ /\n\r|\r\n/) {
#Do Nothing
} else {
for ($lineSearchDupe){
s/[\r]/\r\n/;
}
}
}

foreach my $lineSearchDupe2 (@fastaArrayTemp) {
if ($lineSearchDupe2 =~ /\n\r|\r\n/) {
#Do Nothing
} else {
for ($lineSearchDupe2){
s/[\n]/\r\n/;
}
}
}

print $fastaChangeFH @fastaArrayTemp;
close $fastaChangeFH;


print $logfastaDataFH "Fasta File: $fileStart$fastaFiles[$itter]\n";
print $logfastaDataFH "Accession Sizes (number of characters): \n";

my @countArrayTemp;

foreach my $lineSearch2 (@fastaArrayTemp) {
if ($lineSearch2 =~ /^$fileStart/) {
my $fastaLength = length($lineSearch2);
my $fastaLengthFix = length($lineSearch2) - 1;
push @countArrayTemp, "$fastaLength\n";
print $logfastaDataFH "$fastaLengthFix\n";
} else {

}
}

my @sortedCountArrayTemp = sort { $a <=> $b } @countArrayTemp;
my $maxCharCount = $sortedCountArrayTemp[-1];

my $maxCharCountFix = $maxCharCount - 1;

print $logfastaCountFH "Fasta File: $fastaFiles[$itter]\n";
print $logfastaCountFH "$maxCharCountFix\n";

push @fastaFinalCountArray, $maxCharCount;


undef @fastaArrayTemp;


#####################################
#Print fasta character counts to log#
#####################################

print $logcharCountFH "##################################################\n";
print $logcharCountFH "$fastaFiles[$itter]\n";
print $logcharCountFH "##################################################\n\n";

#Determine Data Type:
if ($nucleotideSumDNA < $aminoAcidSum && $nucleotideSumRNA){
$dataInputType = "Amino Acids / Protien";
}elsif ($aminoAcidSum = $nucleotideSumDNA && $uracilCharDataCounter = 0) { 
$dataInputType = "Nucleotide / DNA";
}elsif ($aminoAcidSum = $nucleotideSumRNA && $uracilCharDataCounter > 0){
$dataInputType = "Nucleotide / RNA";
}
}else{
die ("Data contains both amino acids and nucleotides");	
}

print $logfastaDebugFH "Data Type: $dataInputType\n\n";

if ($dataInputType eq "Nucleotide / DNA"){
#Nucleotides DNA
print $logcharCountFH "Nucleotide Counts\n\n";
print $logcharCountFH "Adenine: $adenineCharDataCounter\n";
print $logcharCountFH "Cytosine: $cytosineCharDataCounter\n";
print $logcharCountFH "Thymine: $thymineCharDataCounter\n";
print $logcharCountFH "Guanine: $guanineCharDataCounter\n\n";
print $logcharCountFH "Total Nucleotide Length: $nucleotideSumDNA\n\n\n";

}elsif ($dataInputType eq "Nucleotide / RNA"){
#Nucleotides RNA
print $logcharCountFH "Nucleotide Counts\n\n";
print $logcharCountFH "Adenine: $adenineCharDataCounter\n";
print $logcharCountFH "Cytosine: $cytosineCharDataCounter\n";
print $logcharCountFH "Guanine: $guanineCharDataCounter\n\n";
print $logcharCountFH "Uracil: $uracilCharDataCounter\n\n";
print $logcharCountFH "Total Nucleotide Length: $nucleotideSumRNA\n\n\n";
}

} elsif ($dataInputType eq "Amino Acids / Protien"){
#Amino Acids
print $logcharCountFH "Amino Acid Counts\n\n";
print $logcharCountFH "Alanine: $alaCharDataCounter\n";
print $logcharCountFH "Arginine: $argCharDataCounter\n";
print $logcharCountFH "Asparagine: $asnCharDataCounter\n";
print $logcharCountFH "Aspartic Acid: $aspCharDataCounter\n";
print $logcharCountFH "Cysteine: $cysCharDataCounter\n";
print $logcharCountFH "Glutamine: $glnCharDataCounter\n";
print $logcharCountFH "Glutamic Acid: $gluCharDataCounter\n";
print $logcharCountFH "Glycine: $glyCharDataCounter\n";
print $logcharCountFH "Histidine: $hisCharDataCounter\n";
print $logcharCountFH "Isoleucine: $ileCharDataCounter\n";
print $logcharCountFH "Leucine: $leuCharDataCounter\n";
print $logcharCountFH "Lysine: $lysCharDataCounter\n";
print $logcharCountFH "Methionine: $metCharDataCounter\n";
print $logcharCountFH "Phenylalanine: $pheCharDataCounter\n";
print $logcharCountFH "Proline: $proCharDataCounter\n";
print $logcharCountFH "Pyrrolysine: $pylCharDataCounter\n";
print $logcharCountFH "Serine: $serCharDataCounter\n";
print $logcharCountFH "Selenocysteine: $secCharDataCounter\n";
print $logcharCountFH "Threonine: $thrCharDataCounter\n";
print $logcharCountFH "Tryptophan: $trpCharDataCounter\n";
print $logcharCountFH "Tyrosine: $tyrCharDataCounter\n";
print $logcharCountFH "Valine: $valCharDataCounter\n";
print $logcharCountFH "Aspartic Acid / Asparagine: $asxCharDataCounter\n";
print $logcharCountFH "Glutamic Acid / Glutamine: $glxCharDataCounter\n";
print $logcharCountFH "Any Amino Acid (non-specified): $xaaCharDataCounter\n";
print $logcharCountFH "Leucine or Isoleucine: $xleCharDataCounter\n\n";
print $logcharCountFH "Total Amino Acid Length: $aminoAcidSum\n\n\n\n";
} else {
print $logcharCountFH "Check .fasta input files, data type not detected.";
die ("Error, check DataTypeLog.txt");
}
$itter ++;
}


my @sortedMaxArray = sort { $a <=> $b } @fastaFinalCountArray;
my $goldenNumber = $sortedMaxArray[-1];

my $goldenNumberFix = $goldenNumber - 1;

print $logfastaCountFH "\n\n\n";
print $logfastaCountFH "The longest acession string is:\n";
print $logfastaCountFH $goldenNumberFix;


#Reopen newly created fasta edit files and adjust accession length
my @fastaEditArrayTemp;

my @charPlaceHold = ( 'X' );

opendir(DIR2,"$fastaEditDir") or die "Couldn't open $fastaEditDir\n";

my @fastaEditFiles = grep(/\.fasta/,readdir(DIR2));
my $numEditFasta = scalar(@fastaEditFiles);

closedir(DIR2);

my $itter2 = 0;

foreach (1 .. $numEditFasta) {

open my $fastaFH2, '<', "$fastaEditDir/$fastaEditFiles[$itter2]" or die "Cannot open fasta file: $fastaEditFiles[$itter2]: $!";
while(<$fastaFH2>) {
push @fastaEditArrayTemp, $_;
}
close $fastaFH2;

foreach my $lineEditSearch (@fastaEditArrayTemp) {

if ($lineEditSearch =~ /^$fileStart/) {

for ($lineEditSearch){

my $addNumOfChar = $goldenNumber - length($lineEditSearch);
my $stringSize = join '' => map $charPlaceHold[rand @charPlaceHold], 1 .. $addNumOfChar;

s/\*/$stringSize/g;

}
} else {

}
}


open my $fastaChangeFH2, '>', "$fastaEditDir2/$fastaFiles[$itter2]" or die "Cannot open $fastaFiles[$itter2]: $!";
print $fastaChangeFH2 @fastaEditArrayTemp;
close $fastaChangeFH2;

fcopy("$fastaEditDir2/$fastaFiles[$itter2]","$fastaEditDir3/$fastaFiles[$itter2]");

fmove("$fastaEditDir2/$fastaFiles[$itter2]","$fastaDir/$fastaFiles[$itter2]");

undef @fastaEditArrayTemp;

$itter2 ++;
}



#Clean Up


rmtree $fastaEditDir;
rmtree $fastaEditDir2;









