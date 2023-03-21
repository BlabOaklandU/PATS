#!/bin/perl.exe
use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Find;
use File::Copy;
use File::Path;



my $baseDir = $ENV{'PoPipe_srcBaseDir'};
my $moleClockBaseDir =  "$baseDir/molecularClocks";
my $moleClockScriptsDir = "$moleClockBaseDir/Scripts";

use lib "$moleClockScriptsDir";
use moleClockModule;


print "Running ctlMakerScript.pl... \n";

##################################################################
# USER INPUTS                                                    #
##################################################################
#
#

#Location of Calibration Tree A
my $caliA_loc = $splitarray[30];
chomp($caliA_loc);

#Location of Calibration Tree B
my $caliB_loc = $splitarray[31];
chomp($caliB_loc);

#Location of Genes
my $genes_loc = $splitarray[20];
chomp($genes_loc);

#Outfile Name A tree
my $outfile_nameA = $splitarray[24];
chomp($outfile_nameA);
#print "outfileA: $outfile_nameA\n";

#Outfile Name B tree
my $outfile_nameB = $splitarray[25];
chomp($outfile_nameB);
#print "outfileB: $outfile_nameB\n";

#Root Age
my $root_age = $splitarray[29];
chomp($root_age);

#Name of root output gene folder (ex. 765phyl_Gene_)
my $output_gene_root = $splitarray[21];
chomp($output_gene_root);

#Number of Genes
my $num_genes = $splitarray[22];
chomp($num_genes);

#Genes to be skipped, seperated by spaces
my $skip_genes = $splitarray[23];
chomp($skip_genes);

# SPLITS MULTIPLE GENES TO BE SKIPPED INTO ARRAY 
my @genes_to_skip = split(/ /, $skip_genes);
#
#
#DONE OBTAINING USER INPUTS
##################################################################


##################################################################
# ACCESSING FILES WITHIN EACH FOLDER                             #
##################################################################
#
#
#
#open location of gene folders
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
		my $phy_file_loc = $tree_gene_loc."/".$tree_file_name;
		opendir(my $phy_dh, $phy_file_loc) or die $!;
		
		#Each .phy file (only iterates through once, for 1 .phy file)
		while (my $phy_file_name = readdir($phy_dh)){
			#skipping every file that is not the .phy file
			if ($phy_file_name eq ".." or $phy_file_name eq "."){
				next;
			}
			elsif(index($phy_file_name, ".txt") != -1){
				next;
			}
			elsif(index($phy_file_name, ".csv") != -1){
				next;
			}
			elsif(index($phy_file_name, ".meg") != -1){
				next;
			}
			elsif(index($phy_file_name, ".ctl") != -1){
				next;
			}
			elsif(index($phy_file_name, ".pl") != -1){
				next;
			}
			
			my $seq_file_name = $phy_file_loc."/".$phy_file_name;
			
			###################################################
			# READ AND STORE GAMMA VALUES FROM TXT 		      #
			###################################################
			#
			#
			my $gamma_txt_name = $phy_file_loc."/".$phy_file_name;
			#replace .phy with .txt
			my $find_phy = "[.]phy";
			my $rep_txt = ".txt";
			$gamma_txt_name =~ s/$find_phy/$rep_txt/g;
				
			open(my $txt_fh, "<", $gamma_txt_name) or die "Could not open file $gamma_txt_name";
			chomp(my @gamma_lines = <$txt_fh>);
			my $rgene_values = $gamma_lines[7];
			my $sigma_values = $gamma_lines[10];
			
			close $txt_fh;
			#
			#
			#DONE READING AND STORING GAMMA VALUES FROM TXT
			####################################################
			
			
			###########################################################
			# CREATE CTL FILE 1 AND 2 AND WRITE TO THEM SIMULTANEOUSLY#
			###########################################################
			#
			#
			my $ctl_1_name = $phy_file_loc."/"."mcmctree1.ctl";
			my $ctl_2_name = $phy_file_loc."/"."mcmctree2.ctl";
			
			open(my $ctl1_fh, ">", $ctl_1_name) or die "Could not open file '$ctl_1_name'";
			open(my $ctl2_fh, ">", $ctl_2_name) or die "Could not open file '$ctl_2_name'";
			
			#seqfile 
			print $ctl1_fh "        seqfile = ",$seq_file_name,"\n";
			print $ctl2_fh "        seqfile = ",$seq_file_name,"\n";
			#treefile, must specify between a and b tree
			if(index($seq_file_name,"A_tree") != -1){
				print $ctl1_fh "      treefile = ",$caliA_loc,"\n";
				print $ctl2_fh "      treefile = ",$caliA_loc,"\n";
				
				print "outfilefilenameA: $outfile_name_pathA\n";
				print $ctl1_fh "      outfile = ",$outfile_name_pathA,"\n";
				print $ctl2_fh "      outfile = ",$outfile_name_pathA,"\n";
			}
			else{
				print $ctl1_fh "      treefile = ",$caliB_loc,"\n";
				print $ctl2_fh "      treefile = ",$caliB_loc,"\n";
				
				print "outfilefilenameB: $outfile_name_pathB\n";
				print $ctl1_fh "      outfile = ",$outfile_name_pathB,"\n";
				print $ctl2_fh "      outfile = ",$outfile_name_pathB,"\n";
			}
			
			#seed 
			print $ctl1_fh "          seed = -1","\n\n";
			print $ctl2_fh "          seed = -1","\n\n";
			
			#ndata 
			print $ctl1_fh "         ndata = 1","\n";
			print $ctl2_fh "         ndata = 1","\n";
			
			#seqtype 
			print $ctl1_fh "       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs","\n";
			print $ctl2_fh "       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs","\n";
			
			#usedata
			print $ctl1_fh "       usedata = 3    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)","\n";
			print $ctl2_fh "       usedata = 2    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)","\n";
			
			#clock 
			print $ctl1_fh "         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates","\n";
			print $ctl2_fh "         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates","\n";
			
			#rootage 
			print $ctl1_fh "       RootAge = '<$root_age'  * safe constraint on root age, used if no fossil for root.","\n\n";
			print $ctl2_fh "       RootAge = '<$root_age'  * safe constraint on root age, used if no fossil for root.","\n\n";
			
			#model 
			print $ctl1_fh "         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85","\n";
			print $ctl2_fh "         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85","\n";
			
			#alpha 
			print $ctl1_fh "         alpha = 0    * alpha for gamma rates at sites","\n";
			print $ctl2_fh "         alpha = 0    * alpha for gamma rates at sites","\n";
			
			#ncatG 
			print $ctl1_fh "         ncatG = 5    * No. categories in discrete gamma","\n\n";
			print $ctl2_fh "         ncatG = 5    * No. categories in discrete gamma","\n\n";
			
			#cleandata 
			print $ctl1_fh "     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?","\n\n";
			print $ctl2_fh "     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?","\n\n";
			
			#BDparas 
			print $ctl1_fh "       BDparas = 2 2 0.1   * birth, death, sampling","\n";
			print $ctl2_fh "       BDparas = 2 2 0.1   * birth, death, sampling","\n";
			
			#rgene_gamma 
			print $ctl1_fh "   rgene_gamma = ",$rgene_values,"\n";
			print $ctl2_fh "   rgene_gamma = ",$rgene_values,"\n";
			
			#kappa_gamma 
			print $ctl1_fh "   kappa_gamma = 6 2      * gamma prior for kappa","\n";
			print $ctl2_fh "   kappa_gamma = 6 2      * gamma prior for kappa","\n";
			
			#alpha_gamma
			print $ctl1_fh "   alpha_gamma = 1 1      * gamma prior for alpha","\n\n";
			print $ctl2_fh "   alpha_gamma = 1 1      * gamma prior for alpha","\n\n";
			
			#note
			print $ctl1_fh "     * gamma prior for rate for genes\n";
			print $ctl2_fh "     * gamma prior for rate for genes\n";
			
			#sigma2_gamma
			print $ctl1_fh "  sigma2_gamma = $sigma_values    * gamma prior for sigma^2     (for clock=2 or 3)\n\n\n";
			print $ctl2_fh "  sigma2_gamma = $sigma_values    * gamma prior for sigma^2     (for clock=2 or 3)\n\n\n";
			
			#finetune
			print $ctl1_fh "      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, rates, mixing, paras, RateParas, FossilErr\n\n";
			print $ctl2_fh "      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, rates, mixing, paras, RateParas, FossilErr\n\n";
			
			#print
			print $ctl1_fh "         print = 1\n";
			print $ctl2_fh "         print = 1\n";
			#burnin
			print $ctl1_fh "        burnin = 4000\n";
			print $ctl2_fh "        burnin = 4000\n";
			#sampfreq
			print $ctl1_fh "      sampfreq = 2\n";
			print $ctl2_fh "      sampfreq = 2\n";
			#nsampple
			print $ctl1_fh "       nsample = 40000\n\n";
			print $ctl2_fh "       nsample = 40000\n\n";
			
			#note
			print $ctl1_fh "*** Note: Make your window wider (100 columns) before running the program.\n";
			print $ctl2_fh "*** Note: Make your window wider (100 columns) before running the program.\n";
			
			close $ctl1_fh;
			close $ctl2_fh;
			#
			#
			#DONE WRITING TO BOTH CTL FILES
			##########################################################################################			
		};
		closedir($phy_dh);
	}
	closedir($tree_genes_dh);	
}
#
#
#DONE ACCESSING FILES IN EACH FOLDER
###################################################################





