#!/bin/perl.exe
use strict;
use warnings;

use Bio::TreeIO;
use Bio::Tree::Node;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::AnnotatableNode;
use List::MoreUtils qw(firstidx);
use List::Util qw(first);
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;


my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $mcBaseDir =  $ENV{'MolecularClockBaseDir'};
my $mcScriptsDir = $ENV{'MolecularClockScriptDir'};
my $mcsplitGenesDir = $ENV{'splitGenesDir'};

open my $logMolClockDebug, '>>', "$mcBaseDir/molClockDebugLog.txt" or die "Cannot open molClockDebugLog.txt: $!";

print $logMolClockDebug "\nRunning gammaCalc.pl... \n";

my $phy_string = ".phy";
my $mccout = "megacc_Out";

#Path of .mao file
my $mao_loc = $mcBaseDir/MEGA-X_PATS_MAO.mao;
chomp($mao_loc);

#Location of Directory of .meg files
my $genes_loc = $mcsplitGenesDir;
chomp($genes_loc);

#Ingroup to Tip Distance (A_tree)
my $ingroup_tip_dist_A = $ENV{'ingroupTipDistanceTreeA'};
chomp($ingroup_tip_dist_A);

#Ingroup to Tip Distance (B_tree)
my $ingroup_tip_dist_B = $ENV{'ingroupTipDistanceTreeB'};
chomp($ingroup_tip_dist_B);

#Name of root output gene folder (ex. 765phyl_Gene_)
my $output_gene_root = $ENV{'rootOutputDirName'};
chomp($output_gene_root);

#Number of Genes
my $num_genes = $ENV{'numberOfGenesMC'};
chomp($num_genes);

#Genes to be skipped, seperated by spaces
my $skip_genes = $ENV{'numberOfGenesSkippedMC'};
chomp($skip_genes);

#Outfile Name A tree
my $outfile_nameA = $ENV{'outFileTreeA'};
chomp($outfile_nameA);

#Outfile Name B tree
my $outfile_nameB = $ENV{'outFileTreeB'};
chomp($outfile_nameB);

# SPLITS MULTIPLE GENES TO BE SKIPPED INTO ARRAY 
my @genes_to_skip = split(/ /, $skip_genes);

my $gene_count = 1;

#read each GENE file in location (iterates through all gene files in folder)
my $genes_dir = $genes_loc;
opendir (my $genes_dh, $genes_dir) or die $!;
my $skip_num = 0;
for(my $i = 1; $i <= $num_genes; $i++){
	if($skip_num < scalar(@genes_to_skip)){
		if($i == $genes_to_skip[$skip_num]){
			$skip_num += 1;
			next;
		}
	}
	
	my $tree_gene_loc = $genes_loc."/".$output_gene_root.$i;
	my $A_Tree_dir = $tree_gene_loc."/"."A_tree";
	my $B_Tree_dir = $tree_gene_loc."/"."B_tree";
	my $outfile_name_pathA = $A_Tree_dir."/".$outfile_nameA;
	my $outfile_name_pathB = $B_Tree_dir."/".$outfile_nameB;
	
	opendir (my $tree_genes_dh, $tree_gene_loc) or die $!;
	#A-tree and B-tree (iterates twice, through a-tree and b-tree folder)
		while (my $tree_file_name = readdir($tree_genes_dh)){
		if ($tree_file_name eq ".." or $tree_file_name eq "."){
			next;
		}
		
		my $meg_file_loc = $tree_gene_loc."/".$tree_file_name;
		opendir(my $meg_dh, $meg_file_loc) or die "Cannot open $meg_file_loc\n";

	#Each .meg file (iterates through twice, for 1 .meg file and .phy file)
		my @tree_files = readdir $meg_dh;
		foreach my $meg_file_name (@tree_files)
		{
			if ($meg_file_name eq ".." or $meg_file_name eq "."){
				next;
			}
			elsif(index($meg_file_name, ".txt") != -1){
				next;
			}
			elsif(index($meg_file_name, ".csv") != -1){
				next;
			}
			elsif(index($meg_file_name, ".ctl") != -1){
				next;
			}
			elsif(index($meg_file_name, ".pl") != -1){
				next;
			}
			#to skip .phy files 
			if (index($meg_file_name, $phy_string) != -1 or index($meg_file_name, $mccout) != -1) {
				next;
			}
			my $whole_meg_path = $meg_file_loc."/".$meg_file_name;
			my $outfile_root =  $meg_file_loc."/".$meg_file_name."_OMeanDist";
			
			#run mega cc 
			system("$baseDir/externalSoftware/MEGA/11.0.13/megacc -a ${mao_loc} -d ${whole_meg_path} -o ${outfile_root}");
#EDIT!!!!!			
			#OPEN FILE WITH OVERALL MEAN DISTANCE (name_of_current_file.csv)
			my $csv_filename = $meg_file_name;
			my $find_meg = ".meg";
			my $replace_csv = ".csv";
			$csv_filename =~ s/$find_meg/$replace_csv/g;
			$csv_filename = $meg_file_loc."/".$csv_filename;
			open(my $csv_fh, '<:encoding(UTF-8)', $csv_filename) or die "Could not open file '$csv_filename' $!";
			
			#read all contents into array (line by line)
			chomp(my @csv_lines = <$csv_fh>);
			my $find_space = " ";
			my $replace_nothing = "";
			$csv_lines[2] =~ s/$find_space/$replace_nothing/g;
			close $csv_fh;
			
			#Perform Calculation for RATE
			my $rate = 0;
			if(index($csv_filename, "_A_") != -1){
				$rate = $csv_lines[2]/$ingroup_tip_dist_A;
			}
			if(index($csv_filename, "_B_") != -1){
				$rate = $csv_lines[2]/$ingroup_tip_dist_B;
			}
			#USE 2 FOR ALPHA
			my $rgene_alpha = 2;
			#calculation for beta
			my $rgene_beta = $rate/($rate * $rate);
			
			#perform calculation for Brown Mean
			my $brown_mean = 0;
			if(index($meg_file_name, "_A_") != -1){
				$brown_mean = 1.5/$ingroup_tip_dist_A;
			}
			if(index($meg_file_name, "_B_") != -1){
				$brown_mean = 1.5/$ingroup_tip_dist_B;
			}
			#use 1 for alpha
			my $sigma_alpha = 1;
			#calculation for beta
			my $sigma_beta = $brown_mean/($brown_mean * $brown_mean);
			
			#open txt file to output all calculations 
			my $calc_out_file = $csv_filename;
			my $replace_txt = ".txt";
			$calc_out_file =~ s/$replace_csv/$replace_txt/g;
			open(my $calc_fh, ">", $calc_out_file) or die "Could not open file '$calc_out_file'";
			
			print $calc_fh "File: ", $calc_out_file, "\n\n";
			if(index($calc_out_file, "_A_") != -1){
				print $calc_fh "Inputted Ingroup-to-Tip Distance for A-tree: ", $ingroup_tip_dist_A, "\n";
			}
			if(index($calc_out_file, "_B_") != -1){
				print $calc_fh "Inputted Ingroup-to-Tip Distance for B-tree: ", $ingroup_tip_dist_B, "\n";
			}
			
			print $calc_fh "Overall Mean Distance: ", $csv_lines[2], "\n";
			print $calc_fh "Calculated Rate: ", $rate, "\n\n";
			print $calc_fh "Rgene_gamma calculations (alpha beta): \n";
			print $calc_fh $rgene_alpha, " ", $rgene_beta, "\n";
			print $calc_fh "Brown Mean: ", $brown_mean, "\n";
			print $calc_fh "Sigma2_gamma calculations (alpha beta): \n";
			print $calc_fh $sigma_alpha, " ", $sigma_beta, "\n";
			
			close $calc_fh;	
		}
		closedir($meg_dh);
	}
	$gene_count += 1;	
}