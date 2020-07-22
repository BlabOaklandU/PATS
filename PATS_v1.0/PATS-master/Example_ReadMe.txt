######################
#Example for PipeLine#
######################



#########################################
# Standard Run Example - NO Permutation #
#########################################

# Example Folder: 
# Move the fasta files to the pipeline fasta folder

# If you have followed the readMe, the pipeline should be setup and ready to run. 

# To start the example, open up the poPipeSettines.sh:

# 1) Set PoPipe_createArchive=1
# 2) Set PoPipe_folderCleanUp=0
# 3) Set PoPipe_rerunTreesOnly=0
# 4) Set PoPipe_rerunFullAnalysis=0
# 5) Set PoPipe_keepOne=0
# 6) Set poPipe_removeGroup=0
# 7) Set PoPipe_removeOne=0
# 8) Set PoPipe_fastTree=1 (This is a much quicker tree creator than is RAxML)


# Start the PipeLine:

# 1) Open your terminal and navigate to the FullPipeLinePackage directory 
# 2) Source the poPipeSettings.sh file
# 3) Bash the poPipeSettings.sh file
# 4) Bash the poPipeFull.sh file

# The pipe should now be running. This will run a standard run without any permutations




################################
# Keep One Permutation Example #
################################

# Example Folder: 
# Move the fasta files to the pipeline fasta folder
# Move the permutation file (keepOne.txt) to the pipeline permutations folder (overwrite if asked)
# View the permutation file to see the formatting for how to add species to each list
# Each specie addition should be on it own newline


# If you have followed the readMe, the pipeline should be setup and ready to run

# To start the example, open up the poPipeSettines.sh:

# 1) Set PoPipe_createArchive=1
# 2) Set PoPipe_folderCleanUp=1
# 3) Set PoPipe_rerunTreesOnly=0
# 4) Set PoPipe_rerunFullAnalysis=1
# 5) Set PoPipe_keepOne=1
# 6) Set poPipe_removeGroup=0
# 7) Set PoPipe_removeOne=0
# 8) Set PoPipe_fastTree=1 (This is a much quicker tree creator than is RAxML)


# Start the PipeLine:

# 1) Open your terminal and navigate to the FullPipeLinePackage directory. 
# 2) Source the poPipeSettings.sh file
# 3) Bash the poPipeSettings.sh file
# 4) Bash the poPipeFull.sh file

# The pipe should now be running. This will run a Keep One permutation with a rerun of the whole analysis per permutation iteration.





##################################
# Remove One Permutation Example #
##################################

# Example Folder: 
# Move the fasta files to the pipeline fasta folder
# Move the permutation file (removeOne.txt) to the pipeline permutations folder (overwrite if asked)
# View the permutation file to see the formatting for how to add spcies to each list
# Each specie addition should be on it own newline.


# If you have followed the readMe, the pipeline should be setup and ready to run. 

# To start the example, open up the poPipeSettines.sh:

# 1) Set PoPipe_createArchive=1
# 2) Set PoPipe_folderCleanUp=1
# 3) Set PoPipe_rerunTreesOnly=0
# 4) Set PoPipe_rerunFullAnalysis=1
# 5) Set PoPipe_keepOne=0
# 6) Set poPipe_removeGroup=0
# 7) Set PoPipe_removeOne=1
# 8) Set PoPipe_fastTree=1 (This is a much quicker tree creator than is RAxML)


# Start the PipeLine:

# 1) Open your terminal and navigate to the FullPipeLinePackage directory. 
# 2) Source the poPipeSettings.sh file
# 3) Bash the poPipeSettings.sh file
# 4) Bash the poPipeFull.sh file

# The pipe should now be running. This will run a Remove One permutation with a rerun of the whole analysis per permutation iteration.





####################################
# Remove Group Permutation Example #
####################################

# Example Folder: 
# Move the fasta files to the pipeline fasta folder
# Move the permutation file (removeGroup.txt) to the pipeline permutations folder (overwrite if asked)
# View the permutation file to see the formatting for how to add spcies to each list
# Each specie addition should be on it own newline.


# If you have followed the readMe, the pipeline should be setup and ready to run. 

# To start the example, open up the poPipeSettines.sh:

# 1) Set PoPipe_createArchive=1
# 2) Set PoPipe_folderCleanUp=1
# 3) Set PoPipe_rerunTreesOnly=0
# 4) Set PoPipe_rerunFullAnalysis=1
# 5) Set PoPipe_keepOne=0
# 6) Set poPipe_removeGroup=1
# 7) Set PoPipe_removeOne=0
# 8) Set PoPipe_fastTree=1 (This is a much quicker tree creator than is RAxML)


# Start the PipeLine:

# 1) Open your terminal and navigate to the FullPipeLinePackage directory. 
# 2) Source the poPipeSettings.sh file
# 3) Bash the poPipeSettings.sh file
# 4) Bash the poPipeFull.sh file

# The pipe should now be running. This will run a Remove Group permutation with a rerun of the whole analysis per permutation iteration.





