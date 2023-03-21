# PATS Version 1.2.0

# Instructions for Use
This software requires sudo privileges
Please make sure to follow all the install instructions

Within the Example directory, open up the example.txt file as a guide and follow the README

This software comes with Muscle 5.1.0, ProteinoOrtho 6.0.14 (STABLE) & 6.1.7 (EXPERIMENTAL), FastTree 2.1.11, RAxML 8.2.12, BLAST+ 2.13.0+, and Diamond 2.0.15 

Software listed above do not need to be installed
If you are to install a different version of the software listed above, this software will not work
These are the most up-to-date versions of the software listed above

*This software utilizes the threaded version of PERL*



# STEP ONE 
Install 7Zip for archive functionality (install newest version)

	Debian 
		sudo apt install p7zip-full

	Fedora
		sudo dnf install p7zip

	RedHat
		sudo yum install p7zip



# STEP TWO
Load Perl Modules

If you are familiar with CPAN or CPANM, use it to install the other needed packages  
If you prefer to use something other than CPAN, that is fine, just make sure the following packages are installed  



*Update CPAN
	CPAN
	install CPAN
	reload CPAN

*Install CPANM
	cpan App:cpanminus

*Install (update) modules via CPANM
	cpan File::Find
	cpan List::MoreUtils
	cpan File::Copy
	cpan File::Copy::Recursive

*Get PATS ready
Unzip the PATS_v1.0.zip file

	unzip PATS_v1.0.zip

*Assign read/write/execute permissions
The two main directories that need read/write/execute permissions is the External Software and Scripts directories
It is easier and quicker to just chmod the whole pipeline directory

	chmod -R 777 PATS_v1.0

*Navigate to the newly created directory

	cd PATS_v1.0/PATS-master/FullPipeLinePackage
	
	
	
# FASTA STANDARDIZER
This software package comes with a .fasta stadardizer bash script that will get all of your .fasta file data structures ready for the pipeline
Please note that all acession lines within and across all .fasta files have to be the same length (i.e. letters + numbers + spaces + newlines)
Run the fastaStandardizer.sh before utilizng the pipeline

To use the fasta standardizer:
	1) Copy all of the .fasta files that will be apart of the analysis into the fasta folder located within PATS_v1.0/PATS-master/FullPipeLinePackage
	2) Open up a terminal and navigate to the base directory; cd PATS_V1.0/Pats-master/FullPipeLinePackage
	3) Source the poPipeSettings.sh file; source ./poPipeSettings.sh
	4) Run the fasta standardizer bash scrip; bash ./fastaStandardizer.sh
	
This generates two folders within the fasta directory:
	1) fastaOriginal - Houses .bak of the origial .fasta files (However, the .fasta file names have been changed to not include '*', '-', and white space. These three characters have been changed to '_'
	2) fastaChangedBackUp - Houses a copy of the completely edited .fasta files that are ready for the pipeline

The fasta standardizer also pushes a copy of the fully edited fasta files that are ready for analysis into the base fasta folder which allows for eased transition into utilizing the pipeline



# FASTA STANDARDIZER ERRORS
Running of the FASTA STANDARDIZER generates a log folder in the base directory; DataStandardizedLogs

	1) This folder will house the outputs needed to check the validity of the fasta standardizer
	2) This also houses a log file that is needed for the pipeline to run. Do NOT edit or adjust these files.
	
On completion, the fasta standardizer will also either fail (warning message) or pass (success message)



# IMPORTANT INFORMATION
PLEASE NOTE: Anything below that has 'export' is a variable and is important for setup of the pipeline  
In STEP 3 (which is below this section), these will be explained in more detail

This pipeline has three main analysis routes: 
	
	1) Standard run (no permutations)
	2) Standard run (no permutations) + Permutations with ONLY tree creation
	3) Standard run with permutations
	
Either tree building software (RAxML or FastTree) can be selected, dependent on needs 
RAxML takes more time to run where as FastTree is much quicker  
NOTE: Only select one, both cannot be selected

		'export PoPipe_fastTree'
		'export PoPipe_RAxMLTree'


# 1) Standard run (no permutations)
This option will run the pipeline fully and output a single tree (either RAxML or FastTree, user defined)  

 run this option please make sure the settings are as follows:
	
		'export PoPipe_createArchive=1'
		'export PoPipe_folderCleanUp=0'
		'export PoPipe_rerunTreesOnly=0'
		'export PoPipe_rerunFullAnalysis=0'
		'export PoPipe_keepOne=0'
		'export PoPipe_removeOne=0'
		'export poPipe_removeGroup=0'
	
NOTE: If you set folderCleanUp=1 you will not be able to implement option 2.	
			


# 2) Standard run (no permutations) + Permutations with ONLY tree creation
First you must run a standard run with folderCleanUp=0 (same as option 1)  
After the standard run, you can then implement permutations at the tree creation level (this will not re-align the sequences)  
Set either PoPipe_keepOne, PoPipe_removeOne, or poPipe_removeGroup equal to 1 to implement it into the analysis  
NOTE: only one permutation (keepone, removeone, or removegroup) can be active during a single run  

To run this option please make sure the settings are as follows:  
	
		'export PoPipe_createArchive=1'
		'export PoPipe_folderCleanUp=0'
		'export PoPipe_rerunTreesOnly=1'
		'export PoPipe_rerunFullAnalysis=0'
		'export PoPipe_keepOne=1'
		'export PoPipe_removeOne=0'
		'export poPipe_removeGroup=0'
	
NOTE: Make sure to set folderCleanUp=0 or this will not work.  	
Once the run is complete, you can then run the next permutation OR switch to a different tree building software  



# 3) Standard run with permutations
This option will run the pipeline fully and output trees based on the permutations selected  
NOTE: PoPipe_folderCleanUp has to be set to 1, always archive this or you will not see any outputs  
Set either PoPipe_keepOne, PoPipe_removeOne, or poPipe_removeGroup equal to 1 to implement it into the analysis  
NOTE: only one permutation (keepone, removeone, or removegroup) can be active during a run  

To run this option please make sure the settings are as follows:  
	
		'export PoPipe_createArchive=1'
		'export PoPipe_folderCleanUp=1'
		'export PoPipe_rerunTreesOnly=0'
		'export PoPipe_rerunFullAnalysis=1'
		'export PoPipe_keepOne=0'
		'export PoPipe_removeOne=0'
		'export poPipe_removeGroup=0'
			
	
	
# INPUT FILE FORMAT
Input files should be in standard FASTA format. 

The fastaStandardizer script will get all fasta files ready for the pipeline.

# STEP THREE
Editing the poPipeSettings.sh File (config file)

*FIRST and MOST IMPORTANTLY*  
Set your base directory path within poPipeSettings.sh  
You can use any text editor you like (e.g., nano, vim)  

Edit poPipeSettings.sh with nano and change PoPipe_srcBaseDir

	nano poPipeSettings.sh
	'export PoPipe_srcBaseDir=/home/USER/Desktop/PATS_v1.0/PATS-master/FullPipeLinePackage'

*Next set the OrthoSpecies filter threshold*  
This is a value between 0 - 100 that depicts the minimum percentage of species required in an orthologous group  
Example: 60 -- This creates orthologous groups based on 60 - 100 percent of species present in the analysis   
Setting at 100 -- Each orthologous group will include ALL species  
Setting at 0 -- Allows any number of species in an orthologous group (2 - max number of species present)  

	'export PoPipe_speciesFilterCutoff=90'

*Next set the GAP filter threshold*  
This is a value between 0 - 100 that depicts the percentage of gaps allowed in the analysis  
Example: 25 -- This allows up to 25 percent of gaps per column in the analysis   
Setting at 100 -- Removes none of the gaps (all allowed)  
Setting at 0 -- Removes any column where a gap is present (NO gaps allowed)  

	'export PoPipe_gapFilterCutoff=75'
	
*Next set the BLAST clustering connectivity*  
The default value is 0.1
Increasing this value leads to more and smaller clusters
  
	'export PoPipe_proteinOrthoConn=0.1'

*Create an archive of completed pipe (compressed .7z)*  
This will compress all populated files via the pipeline  
When set to 1 -- creates archive; when set to 0 -- does not create archive  
Recommended to always set this to 1  
See section below for more information (Archiving)  

	'export PoPipe_createArchive=1'

*Folder cleanup will delete all folders populated during the analysis* 
This should always be paired with the archiver  
This functionality is needed when wanting to run an analysis via step 3 (see above)  
See section below for more information (File Cleanup)  

	'export poPipe_folderCleanUp=0'

*Run only the tree creation (can be combined with permutations)*  
This will ONLY rerun the trees based off the orthologs previously created from a full run  
A standard run needs to be finished first then this can be set to 1  
Never run this with folderCleanUp=1  
See section below for more information (Rerunning ONLY tree creation)  

	'export poPipe_rerunTreesOnly=0'

*Run full analysis from start to finish (can be combined with permutations)*  
This will permutate then generate the orthologs  
Need to run with folderCleanUp=1  

	'export poPipe_rerunFullAnalysis=0'

*Tree Selection Protocol (fastTree or RAxML)*  
Only select one tree creation protocol (1 = TRUE; 0 = FALSE)  
By default, FastTree will only use the LG model with gamma  

	'export PoPipe_fastTree=1'
	'export PoPipe_RAxMLTree=0'

*Permutation Settings (1 = TRUE; 0 = FALSE)*
**Only one permutation allowed per run**
**No permutations = standard run**  
See section below for more information (Permutations)  

KeepOne permutation will read in the specified list of species and will loop the analysis, iterating through each of the species list, keeping one of the species in the analysis, generating multiple outputs  

	'export poPipe_keepOne=1'

RemoveOne permutation will read in the specified list of species and will loop the analysis, iterating through each of the species list, removing one per analysis, generating multiple outputs  

	'export poPipe_removeOne=1'

RemoveGroup permutation will read in the specified list of species and will run the analysis with the group of species listed, removed  

	'export poPipe_removeGroup=1'



# STEP FOUR
*Source settings file (make sure you are in the directory of poPipeSettings)  
'std' directory CANNOT be present before sourcing poPipeSettings  

To run multiple analyses, ALWAYS make sure to re-source poPipeSettings.sh  

	cd ./FullPipeLinePackage

	source ./poPipeSettings.sh

Once sourced, a folder should be created in the working directory ('std')  

*Run poPipeSettings.sh  

	bash ./poPipeSettings.sh

*Run poPipeFull.sh  
This will start the nalysis  

	bash ./poPipeFull.sh



# Permutations
When running the analysis with any permutations, keepOne.txt/removeOne.txt/removeGroup.txt need to be populated to run that specific analysis  
If you want to run a removeOne permutation, open the removeOne.txt and paste the species file name (including extensions) one per new line  
The species proper names can be found (if you run a standard analysis first) in the logFiles directory in the file SpeciesList.txt  
Then you can simply copy and paste each species on a new line in the removeOne.txt file  

**Keep in mind that one permutation can only be run at a single time**



# Rerunning ONLY tree creation
Rerunning trees only will allow for permutations post-alignment  
Make sure you run the standard run (all permutations set to 0, PoPipe_rerunFullAnalysis=0, and export PoPipe_rerunTreesOnly=0)  
This will allow you to run either RAxML or fastTree on the same species without a problem as well as each permutation  
No overwrite will occur if you run each permutation once  



#Orthologs, Co-Orthologs, and Paralogs
All orthologs, co-orthologs, and paralogs are calculated based upon the users OrthoSpecies filter threshold input

This software will generate a log for each of the following, located in the logFiles directory:
	1) Orthologs    -- Number of species = number of genes AND number of species >= OrthoSpecies filter threshold input
	2) Co-Orthologs -- Number of species < number of genes AND number of species >= OrthoSpecies filter threshold input
	3) Paralogs     -- Number of species < number of genes AND number of species < OrthoSpecies filter threshold input
	
	*Example:
	100 species in analysis, OrthoSpecies filter threshold input set at 60% (cut off at 60 species).
	
	                  Num Species	Num Genes
	1) Ortholog    -- 	 60           60
	2) Co-Ortholog --    60           61
	3) Paralog     --    59           60



# Archiving
Archiving is a very useful tool that will zip (7zip) your completed run and all populated files  
The archive will be named based on the type of analysis run and will be date and time stamped  
To view which species were in the analysis, explore the archive and open up archiveLog.txt  



# File Cleanup
File cleanup is a way to delete all populated folders to minimize amount of space used  
This is needed for a pre-alignment permutation and should ALWAYS be paired with the archive function  
ALL populated folders will be deleted; thus, this cannot be paired with PoPipe_rerunTreesOnly  



# Data Structure 
**Folders:**  

*Pre-created folders:*  

	'fasta'                -- Houses all of the .fasta files that will be utilized in the analysis  
	'permutations'         -- Houses the three .txt files used for permutating (keepone, removeone, and removegroup) the analysis  
	'externalSoftware'     -- Houses proteinOrtho, Muscle, RAxML, and FastTree software  
	'scriptsPerl'          -- Houses the Perl scripts used in the pipeline 
    'DataStandardizerLogs' -- Houses the log outputs after the fastaStandardizer.sh is run.	

*Folders generated during analysis:* 

	'std'                  -- Houses error and out file logs, initial folder created when sourcing the poPipeSettings.sh file   
	'logFiles'             -- Houses all of the logs generated throughout the analysis 
	'poWork'               -- Houses the work files from proteinOrtho  
	'orthologs'            -- Houses the output from protienOrtho (myproject.proteinortho)  
	'og'                   -- Houses multiple outputs from the speciesfilter which are individual ortholog based on the speciesfilter threshold per OG file (does not contain sequences)  
	'orthoGroupsFasta'     -- Houses multiple outputs from the OGFilter, these outputs are the same as the species filter outputs but populated with the genomic data (contains sequences)  
	'muscleFasta'          -- Houses the Muscle aligned outputs, same as the OGFilter outputs, but now they are aligned and contain gaps (contains sequences)  
	'OGGapFilter'          -- Houses the outputs from the GapFilter which are a single orthologous group that has gaps removed based on the gapfilter threshold  
	'cleanedGap'           -- Houses single orthologous group files that have been cleand and properly set with gaps
	'concatFilterIndv'     -- Houses the output from the ConcatFilter, this creates a file per species in the analysis and populates that file with the edited genomic sequences in order  
	'concatFilterGroup'    -- Houses the concatFilterIndv files in a single file based on species  
	'treeRAxML'            -- Houses the RAxML tree output files  
	'fastTree'             -- Houses the FastTree tree output files  

*Files:*

The only file that should be edited is the poPipeSettings.sh file  
If you edit any of the shell or Perl scripts, it could break the program  



# External Software Versions
FastTree Version:      2.1.11
RAxML Version:         8.2.12
ProtienOrtho Version:  6.0.14 (Galaxy/Stable) & 6.1.3 (Experimental/Newest)
MUSCLE Version:        3.8.31 & 5.1.0
BLAST+:		           2.13.0+
Diamond:			   2.0.15