#!/bin/perl.exe
use strict;
use warnings;

my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $mcBaseDir =  $ENV{'MolecularClockBaseDir'};
my $mcScriptsDir = $ENV{'MolecularClockScriptDir'};
my $mcsNullDir = $ENV{'nullDir'};


open my $logMolClockDebug, '>>', "$mcBaseDir/molClockDebugLog.txt" or die "Cannot open molClockDebugLog.txt: $!";


print $logMolClockDebug "Running nul_checker.pl... \n";

my $directory = $mcsNullDir;
chdir($directory) or die;

#############################
# USER INPUTS
#############################
#
#

#Location of SeqGen Sequences
my $genes_loc = $ENV{'seqGenSequencesDir'};
chomp($genes_loc);

#Name of output file
my $output_name = "$mcsNullDir/null_chatacters.txt";
chomp($output_name);
#
#
##############################

#CHARACTER TO LOOK FOR
my $nul_char = "\0";

#OPEN DIRECTORY WHERE GENES ARE LOCATED 
opendir(my $gene_dh, $genes_loc) or die $!;

#open text file to write faulty genes to
open(my $out_fh, ">", $output_name) or die "Could not open file $output_name $!";

#ITERATE THROUGH ALL GENE FILES AND OPEN
while(my $gene_file_name = readdir($gene_dh)){
	if ($gene_file_name eq ".." or $gene_file_name eq "."){
		next;
	}
	my $gene_path = $genes_loc."/".$gene_file_name;
	open my $gene_fh, '<', $gene_path or die "Can't open file $!";
	my $text = do { local $/; <$gene_fh> };
	# IF NUL CHARACTER IS FOUND IN FILE, PRINTS FILE NAME TO TXT DOCUMENT
	if (index($text, $nul_char) != -1){
		print $out_fh $gene_file_name, "\n";
	}
	close $gene_fh;
}

close $out_fh;
closedir($gene_dh);


print $logMolClockDebug "nul_checker completed.\n";