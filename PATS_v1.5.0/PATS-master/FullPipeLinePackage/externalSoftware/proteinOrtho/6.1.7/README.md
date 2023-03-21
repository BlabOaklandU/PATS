# Proteinortho

Proteinortho is a tool to detect orthologous genes within different species.

**Input**: Multiple fasta files (orange boxes) with many proteins/genes (circles). 
**Output**: Groups (\*.proteinortho) and pairs (\*.proteinortho-graph) of orthologs proteins/genes.

For doing so, it compares similarities of given gene sequences and clusters them to find significant groups. 
The algorithm was designed to handle large-scale data and can be applied to hundreds of species at one. 
Details can be found in ([doi:10.1186/1471-2105-12-124](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-124)).
To enhance the prediction accuracy, the relative order of genes (synteny) can be used as additional feature for the discrimination of orthologs. The corresponding extension, namely PoFF ([doi:10.1371/journal.pone.0105015](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105015)), is already build in Proteinortho. The general workflow of proteinortho: 

<img src="https://www.uni-marburg.de/de/fb16/ipc/ag-lechner/graph.png/@@images/image/unimr_lead_image_sd" alt="proteinortho.workflow.png" height="250">

First an initial all vs. all comparison between all proteins of all species is performed to determine protein similarities (upper right image). <br>
The second stage is the clustering of similar genes to meaningful co-orthologous groups (lower right image). <br>
Connected components within this graph can be considered as putative co-orthologous groups in theory and are returned in the output (lower left image).

# New Features of Proteinortho Version 6

  - Implementation of various Blast alternatives for step (for -step=2 the -p= options): Diamond, MMseqs2, Last, Topaz, Rapsearch2, Blat, Ublast and Usearch
  - Multithreading support for the clustering step (-step=3)
  - Integration of the LAPACK Fortran Library for a faster clustering step (-step=3)
  - Integration of the bitscore weights in the connectivity calculation for more data dependant splits (-step=3)
  - Continuous Integration & Continuous Development [![pipeline status](https://gitlab.com/paulklemm_PHD/proteinortho/badges/master/pipeline.svg)](https://gitlab.com/paulklemm_PHD/proteinortho/pipelines) 
<details>
  <summary>New minor features: (Click to expand)</summary>

  - Output now supports OrthoXML (-xml) and HTML.
  - [proteinortho_history.pl](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Tools%20and%20additional%20programs) a new tool for tracking proteins (or pairs of proteins) in the workflow of proteinortho.
  - [proteinortho_summary.pl](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Tools%20and%20additional%20programs)
  - Various test routines (make test).
  - New heuristics for connectivity calculation (-step=3).
</details><details>
  <summary>6.0.12: (Click to expand)</summary>

  - improved [proteinortho_history.pl](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Tools%20and%20additional%20programs) : now the program is "smarter" in detecting files automatically
  - added [proteinortho_summary.pl](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Tools%20and%20additional%20programs) : a tool for summarizing the proteinortho-graph on species level. With the output it is easy to identify weak connected species.   
  - removed the diamond spam
</details>
 <details><summary>6.0.13: (Click to expand)</summary>

  - added -p=autoblast : this option alows the comparison of aminoacid and nucleotide sequences. E.g. Proteom-vs-Genome: find the protein that corresponds to a given gene (/cluster). 
  - added -isoform={ncbi,uniprot,trinity} option : The reciprocal best hit graph is build using isoform information (isoforms are treated equivalent). [more information about --isoform](https://gitlab.com/paulklemm_PHD/proteinortho/-/wikis/FAQ#how-does-the-isoform-work)
</details>

6.0.14 : public release to https://usegalaxy.eu/

<details><summary>6.1.1: (Click to expand)</summary>

  - main overhaul of the multithreading system of proteinortho_clustering.cpp now using c++11 worker threads instead of a plain openMP. The new system now uses the provided cores more efficiently
  - added -core parameter to proteinortho_clustering.cpp (and the main perl script): stop clustering if a group would split into groups that dont span all species of the inital connected component
</details>

A more detailed list of all changes: [CHANGELOG](https://gitlab.com/paulklemm_PHD/proteinortho/blob/master/CHANGELOG)

# Table of Contents
1. [Installation](#installation)
2. [Synopsis and Description](#synopsis)
3. [Options/Parameters](#options)
4. [PoFF synteny extension](#poff)
5. [Output description](#output)
6. [Examples](#examples)

# [Proteinortho-Wiki](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/) Table of Contents

1. [Tools and additional programs](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Tools%20and%20additional%20programs)
2. [Error Codes and Troubleshooting](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error-Codes) <- look here if you cannot compile/run proteinortho
3. [Large compute jobs example](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Large-compute-jobs-(the--jobs-option))
3. [Synteny + Core Proteome Example](https://gitlab.com/paulklemm_PHD/proteinortho/-/wikis/synteny-example)
4. [FAQ](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/FAQ) <br>
[(...)](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/)

Bug reports: Please have a look at chapter [2.](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error-Codes) first or send a mail to incoming+paulklemm-phd-proteinortho-7278443-issue-@incoming.gitlab.com. (please include the 'parameter-vector' that is printed for all errors)
You can also send mails to lechner@staff.uni-marburg.de. Any suggestions, feedback and comments are welcome!


# Installation

 **Proteinortho comes with precompiled binaries of all executables (Linux/x86) so you should be able to run perl proteinortho6.pl in the downloaded directory for Linux/x86.**
You could also move all executables to a local bin directory (e.g. with make install PREFIX=\~/bin).
If you cannot execute the src/BUILD/Linux_x86_64/proteinortho_clustering, then you have to recompile with make, see the section 2. Building and installing proteinortho from source.

<br>

#### Installation with (bio)conda (for Linux + OSX) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/proteinortho/README.html) [![alt](https://img.shields.io/conda/dn/bioconda/proteinortho.svg?style=flat)](https://bioconda.github.io/recipes/proteinortho/README.html)

    conda install -c bioconda proteinortho

If you need conda (see [here](https://docs.anaconda.com/anaconda/install/)) and the bioconda channel: `conda config --add channels bioconda`.

<br>

#### Installation with brew (for OSX) [![install with brew](https://img.shields.io/badge/install%20with-brew-brightgreen.svg?style=flat)](https://formulae.brew.sh/formula/proteinortho) [![dl](https://img.shields.io/badge/dynamic/json.svg?label=downloads&query=$[%27analytics%27][%27install%27][%27365d%27][%27proteinortho%27]&url=https%3A%2F%2Fformulae.brew.sh%2Fapi%2Fformula%2Fproteinortho.json&color=green)](https://formulae.brew.sh/formula/proteinortho)

    brew install proteinortho

If you need brew (see [here](https://brew.sh/index_de))

<br>

#### Deploy with docker [![install with docker](https://img.shields.io/badge/install%20with-docker-brightgreen.svg?style=flat)](https://quay.io/repository/biocontainers/proteinortho)

    docker pull quay.io/biocontainers/proteinortho:TAG

you can find the TAG [here](https://quay.io/repository/biocontainers/proteinortho?tab=tags) (e.g. 6.0.23--hfd40d39_0).




<details>
  <summary>how to docker (Click to expand)</summary>

  <br>

  First define a [TAG](https://quay.io/repository/biocontainers/proteinortho?tab=tags) with:

  ```export TAG='put the version tag here'```

  To start a simple bash shell with the proteinortho container use:

  ```docker run --rm -it quay.io/biocontainers/proteinortho:$TAG bash ```

  Here you can start/use proteinortho.
  You can change "6.0.22--hfd40d39_0" with any tag/version that is available [here](https://quay.io/repository/biocontainers/proteinortho?tab=tags). Sadly there is no ":latest" tag available ...

  ### Now lets try to mount your home in the proteinortho container

  This is neccessary if you want to access your local files:

  ```docker run --rm --mount "type=bind,src=/home/$(id -un),dst=/home/$(id -un)" -u $(id -u):$(id -g) -it quay.io/biocontainers/proteinortho:$TAG bash```

  now you have your home directory mounted to /home/YOURNAME. (load your bashrc within the container : ```source /home/YOURNAME/.bashrc```)


</details>

<br>

#### Available at Galaxy Europe

Simply go to the european galaxy server and search for proteinortho:

    https://usegalaxy.eu

Or you can integrate proteinortho into your own galaxy instance using: [proteinortho (iuc repository)](https://toolshed.g2.bx.psu.edu/view/iuc/proteinortho/4850f0d15f01)

<br>

#### Installation with dpkg (root privileges are required)

Disclamer: Be aware that this method usually lacks 6-12 months behind the latest version

The deb package can be downloaded here: [unstable](https://packages.debian.org/unstable/proteinortho) or [stable](https://packages.debian.org/stable/proteinortho).
Afterwards the deb package can be installed with `sudo dpkg -i proteinortho*deb`.

<br>

#### *Installation with apt-get*

Disclamer: Be aware that this method usually lacks 6-12 months behind the latest version ([current version](https://packages.debian.org/stable/proteinortho))

*proteinortho is released to stable Debian 11 (2021), so you can install it with `(sudo) apt install proteinortho`*

<br>

#### Prerequisites for compiling proteinortho from source

Proteinortho uses standard software which is often installed already or is part of then package repositories and can thus easily be installed. The sources come with a precompiled version of Proteinortho for 64bit Linux x86.

<details>
  <summary>To <b>run</b> Proteinortho, you need: (Click to expand)</summary>


   - At least one of the following the following programs (default is diamond):

     - NCBI BLAST+ or NCBI BLAST legacy (to test this, type tblastn. apt-get install ncbi-blast+)
     - Diamond (apt-get install diamond, brew install diamond, conda install diamond, https://github.com/bbuchfink/diamond)
     - Last (http://last.cbrc.jp/)
     - Rapsearch (https://github.com/zhaoyanswill/RAPSearch2)
     - Topaz (https://github.com/ajm/topaz)
     - usearch (https://www.drive5.com/usearch/download.html)
     - ublast (is part of usearch)
     - blat (http://hgdownload.soe.ucsc.edu/admin/)
     - mmseqs2 (conda install mmseqs2, https://github.com/soedinglab/MMseqs2)
   - Perl v5.08 or higher (to test this, type perl -v in the command line)
   - (optional) Python v3.0 or higher to include synteny analysis (to test this, type 'python -V' in the command line)
   - Perl standard modules (these should come with Perl): Thread::Queue, File::Basename, Pod::Usage, threads (if you miss one just install with `cpan install ...` )
</details>

<br>
<details>
  <summary>To <b>compile</b> Proteinortho (linux/osx), you need: (Click to expand)</summary>

   - GNU make (to test this, type 'make' in the command line)
   - GNU g++ v4.1 or higher (to test this, type 'g++ --version' in the command line)
   - openmp (to test this, type 'g++ -fopenmp' in the command line)
   - (optional) gfortran for compiling LAPACK (to test this, type 'whereis gfortran' in the command line)
   - (optional) CMake for compiling LAPACK (to test this, type 'cmake' in the command line), OR you can use your own compiled version of lapack (you can get this with 'apt-get install liblapack3') and run 'make USEPRECOMPILEDLAPACK=TRUE'

</details>

<br>

#### Building and installing proteinortho from source (linux and osx)

  You need to have a working lapack library, check this e.g. with 'dpkg --get-selections | grep lapack'. Install lapack e.g. with 'apt-get install libatlas3-base' or liblapack3.

  If you dont have Lapack, then 'make' will automatically compiles an old lapack (v3.8.0) for you automatically as fallback !

  Fetch the latest source code archive downloaded from <a href="https://gitlab.com/paulklemm_PHD/proteinortho/-/archive/master/proteinortho-master.zip">here</a>
<details> <summary>or from here (Click to expand)</summary>

  > git clone https://gitlab.com/paulklemm_PHD/proteinortho

  > wget https://gitlab.com/paulklemm_PHD/proteinortho/-/archive/master/proteinortho-master.zip
</details>
<br>

  - `tar -xzvf proteinortho*.tar.gz` or `unzip proteinortho*.zip` : Extract the files
  - `cd proteinortho*` : Change directory into the extracted folder
  - You can now run proteinortho6.pl directly (linux only).
  - `make clean && make` : If you want to recompile Proteinortho. (For osx you need a newer g++ compiler to support multithreading, see below)
  - `make test` : To make sure Proteinortho works as expected. The output should look like below ('Make test output').
  - `make install` or `make install PREFIX=~/bin` if you dont have root privileges.

<details>
  <summary><b>OSX additional informations (the -fopenmp error)</b></summary>
<pre>
Install a newer g++ compiler for -fopenmp support (multithreading) with brew (get brew here https://brew.sh/index_de)

```
brew install gcc --without-multilib
```

Then you should have a g++-7 or whatever newer version that there is (g++-8,9,...).
Next you have to tell make to use this new compiler with one of the following:
```
ln -s /usr/local/bin/gcc-7 /usr/local/bin/gcc
ln -s /usr/local/bin/g++-7 /usr/local/bin/g++
```

OR(!) specify the new g++ in 'make CXX=/usr/local/bin/g++-7 all'
</pre>
</details>

<details>
  <summary>'make' successful output (Click to expand)</summary>
<pre>
[  0%] Prepare proteinortho_clustering ...
[ 20%] Building **proteinortho_clustering** with LAPACK (static/dynamic linking)
[ 25%] Building **graphMinusRemovegraph**
[ 50%] Building **cleanupblastgraph**
[ 75%] Building **po_tree**
[100%] Everything is compiled with no errors.
</pre>

The compilation of proteinortho\_clustering has multiple fall-back routines. If everything fails please look here [Troubleshooting (proteinortho wiki)](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error%20Codes).

</details>

#### Make test output

<details>
  <summary>'make test' successful output (Click to expand)</summary>
<pre>
Everything is compiled with no errors.
[TEST] 1. basic proteinortho6.pl -step=2 tests
 [1/11] -p=blastp+ test: passed
 [2/11] -p=blastp+ synteny (PoFF) test: passed
 [3/11] -p=diamond test: passed
 [4/11] -p=diamond (--moresensitive) test (subparaBlast): passed
 [5/11] -p=lastp (lastal) test: passed
 [6/11] -p=topaz test: passed
 [7/11] -p=usearch test: passed
 [8/11] -p=ublast test: passed
 [9/11] -p=rapsearch test: passed
 [10/11] -p=blatp (blat) test: passed
 [11/11] -p=mmseqsp (mmseqs) test: passed
[TEST] 2. -step=3 tests (proteinortho_clustering)
 [1/2] various test functions of proteinortho_clustering (-test): passed
 [2/2] Compare results of 'with lapack' and 'without lapack': passed
[TEST] Clean up all test files...
[TEST] All tests passed
</pre>
</details>

If you have problems compiling/running the program go to [Troubleshooting (proteinortho wiki)](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error%20Codes).

<br>

# SYNOPSIS
  > **proteinortho [options] \<fasta file(s)\>**
 
   one fasta file for each input species; at least two species are required

# DESCRIPTION
  **proteinortho** is a tool to detect orthologous genes within different
  species. 

  Proteinortho assumes, that you have all your gene sequences in FASTA
  format either represented as amino acids or as nucleotides. The source
  code archive contains some examples, namely C.faa, E.faa, L.faa, M.faa
  located in the test/ directory. **By default Proteinortho assumes amino**
  **acids sequences and thus uses diamond** (-p=diamond) to compare sequences. If you have
  nucleotide sequences, you need to change this by adding the parameter
  -p=blastn+ (or some other algorithm). (In case you have only have NCBI
  BLAST legacy installed, you need to tell this too - either by adding
  -p=blastp or -p=blastn respectively.) The full command for the example
  files would thus be
  > proteinortho6.pl -project=test test/C.faa test/E.faa

  test/L.faa test/M.faa. Instead of naming the FASTA files one by one, you
  could also use test/*.faa. Please note that the parameter
  -project=test is optional, for naming the output. With this, you can set the prefix of the output
  files generated by Proteinortho. If you skip the project parameter, the
  default project name will be myproject.

# OPTIONS graphical user interface

Open `proteinorthoHelper.html` in your favorite browser or visit [lechnerlab.de/proteinortho](http://lechnerlab.de/proteinortho/) online for an interactiv exploration of the different options of proteinortho.

# OPTIONS

 **Main parameters** (can be used with -- or -)

   - **--project**=name (default: myproject)
    prefix for all resulting file names

   - **--inproject**=name (default: same as --project)
    load data from this namespace instead (works with intermediate files for step=2 and blast-graph for step=3)

   - **--cpus**=number (default: all available)
    the number of processors to use (multicore/processor support)

  - **--ram**=number (default: 90% of free memory)
    maximal used ram threshold for LAPACK and the input graph in MB

  - **--verbose**={0,1,2} (default: 1)
    verbose level. 1:keeps you informed about the progress

  - **--silent**
    sets verbose level to 0.

  - **--temp**=directory(.)
    path to the temporary files

  - **--force**
    forces the recalculation of the blast results in any case in step=2. Also forces the recreation of the database generation in step=1

  - **--clean**
    removes all database-index-files generated by the -p algorithm afterwards

  - **--step**={0,1,2,3} (default: 0)
    0 -> all. 1 -> prepare blast (build db). 2 -> run all-versus-all
    blast. 3 -> run the clustering.

    <details>
      <summary>(Show more information)</summary>
        
        proteinortho test/*faa 
      
        # the following 3 commands are producing the same results as the command above
        proteinortho -step=1 test/*faa 
        proteinortho -step=2 test/*faa 
        proteinortho -step=3
      
    </details>   
    
  - **--keep**
    stores temporary blast results for reuse (proteinortho_cache_project directory). 
    In a second run the intermediate blast results are loaded instead of calculated.
    You can adjust the parameters e.g. a more strict -e evalue cut off and write the output to a different namespace using --inproject.

    <details>
      <summary>(Show more information)</summary>

        # 1. generate db files
        
        proteinortho -step=1 -project=test -keep infile/*fasta
        
        # 2. run the all-versus-all blast of some input files (infile/) 
        
        proteinortho -step=2 -project=test -keep infile/*fasta
        
        # now you can insert more fasta files to infile/ and reuse everything computed 
        
        proteinortho -step=2 -project=test -keep infile/*fasta
        
        # finally run clustering
        
        proteinortho -step=3 -project=test -keep

    </details>     
        
  - **--isoform**={ncbi,uniprot,trinity} [more information about --isoform](https://gitlab.com/paulklemm_PHD/proteinortho/-/wikis/FAQ#how-does-the-isoform-work)
   
    Merge isoforms to a single entity. 

    <details><summary>ncbi</summary> 
        
        isoforms are specified in ncbi style 
        
        ---
        >**ENSMUSP00000021091.8** pep chromosome:GRCm38:11:74673949:74724670:-1 **gene:ENSMUSG00000020745.15** transcript:ENSMUST00000021091.14 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Pafah1b1 description:platelet-activating factor acetylhydrolase, **isoform** 1b, subunit 1 [Source:MGI Symbol;Acc:MGI:109520]
        >**ENSMUSP00000099578.2** pep chromosome:GRCm38:11:74673950:74723858:-1 **gene:ENSMUSG00000020745.15** transcript:ENSMUST00000102520.8 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Pafah1b1 description:platelet-activating factor acetylhydrolase, **isoform** 1b, subunit 1 [Source:MGI Symbol;Acc:MGI:109520]  
        ---
        
        Different protein identifier (ENSMUSP00000021091.8, ENSMUSP00000099578.2) but the same gene id (ENSMUSG00000020745.15). The word 'isoform' is also mandatory!
        
    </details>
    <details><summary>uniprot</summary> 
        
        isoforms are specified in uniprot style using the '*_additional.fa' files 
        
        E.g. C.fa: 
        
        ---
        >tr|ADHA2|R4GDP1_DANRE Gamma-aminobutyric
        (...)
        ---
        
        C_additional.fa: 
        
        ---
        >tr|QDHQ4|R4GDP1_DANRE isoform of ADHA2
        (...)
        ---
        
        QDHQ4 is the isoform of ADHA2. Please simply add the *_additional.fa files to the proteinortho call!
        
    </details>
    <details><summary>trinity</summary> 
        
        isoforms are specified in trinity style:
        
        ---
        >TRINITY_DN1000_c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]
        (...)
        ---
        
        The protein id is TRINITY_DN1000_c115_g5a and the isoform id is specified with i1
        
    </details>

 **Search options (step 1-2)**
  (output: <myproject>.blast-graph)

  - **--p**=algorithm (default: diamond) 

    <details>
      <summary>show all options (Click to expand)</summary>

        - autoblast : automatically detects the blast+ program (blastp,blastn,tblastn,blastx) depending on the input (can also be mixed together!)

        - blastn_legacy,blastp_legacy,tblastx_legacy : legacy blast family (shell commands: blastall -) family. The suffix 'n' or 'p' indicates nucleotide or protein input files.

        - blastn+,blastp+,tblastx+ : standard blast family (shell commands: blastn,blastp,tblastx)
        family. The suffix 'n' or 'p' indicates nucleotide or protein input files.

        - diamond : Only for protein files! standard diamond procedure and for
        genes/proteins of length >40 with the additional --sensitive flag
        Warning: Please use version 0.9.29 or later to avoid this known bug: #24

        - lastn,lastp : lastal. -n : dna files, -p protein files (BLOSUM62 scoring matrix)!

        - rapsearch : Only for protein files!

        - mmseqsp,mmseqsn : mmseqs2. -n : dna files, -p protein files

        - topaz : Only for protein files!

        - usearch : usearch_local procedure with -id 0 (minimum identity
        percentage).

        - ublast : usearch_ublast procedure.

        - blatp,blatn : blat. -n : dna files, -p protein files
    </details>
    <br>

  - **--sim**=float (default: 0.95)
    min. reciprocal similarity for additional hits. 1 : only the best reciprocal hits are reported, 0 : all possible reciprocal blast matches (within the -e) are reported.

<details>
  <summary>More (Click to expand)</summary>

  - **--e**=evalue (default: 1e-05)
    E-value for blast
    (column 11 of blast outfmt 6 output)

  - **--selfblast**
    apply selfblast, detects paralogs without orthologs

  - **--identity**=number (default: 25)
    min. percent identity of best blast hits
    (column 3 (pident) of blast outfmt 6 output)

  - **--cov**=number (default: 50)
    min. coverage of best blast alignments in %
    coverage between protein A and B = min ( alignment_length_A_B/length_A, alignment_length_A_B/length_B )
    (alignment_length_A_B = column 4 of blast outfmt 6 output)

  - **--subparaBlast**='options'
    additional parameters for the search tool (-p=blastp+,diamond,...) example -subpara='-seg no'
    or -subpara='--more-sensitive' for diamond

  - **--identical**
    only return entries that are 100% identical

  - **--range**=number (default:disabled)
    maximal length difference for any blast hit. e.g. 0 = filter for hits between proteins of same length

</details>
<br>

 **Synteny options (optional, step 2)**
  (output: <myproject>.ffadj-graph, <myproject>.poff-graph, <myproject>.poff.tsv (tab separated file with groups))

<details>
  <summary>More (Click to expand)</summary>

  - **--synteny**
    activate PoFF extension to separate similar by contextual adjacencies
    (requires .gff for each .fasta)

  - **--dups**=number (default: 0)
    PoFF: number of reiterations for adjacencies heuristic, to determine
    duplicated regions

  - **--cs**=number (default: 3)
    PoFF: Size of a maximum common substring (MCS) for adjacency matches

  - **--alpha**=number (default: .5)
    PoFF: weight of adjacencies vs. sequence similarity
</details>
<br>

 **Clustering options (step 3)**
  (output: <myproject>.proteinortho.tsv, <myproject>.proteinortho.html, <myproject>.proteinortho-graph)

  - **--conn**=float (default: 0.1)
    min. algebraic connectivity. <b>This is the main parameter for the clustering step.</b> Choose larger values then more splits are done, resulting in more and smaller clusters. (There are still cluster with an alg. conn. below this given threshold allowed if the protein to species ratio is good enough, see -minspecies option below)

<details>

  <summary>More (Click to expand)</summary>

  - **--singles**
    report singleton genes without any hit

  - **--purity**=float (default: 1e-7)
    avoid spurious graph assignments

  - **--minspecies**=float (default: 1, must be >=0)
    min. number of genes per species. If a group is found with up to (minspecies) genes/species, it wont be split again (regardless of the connectivity).

  - **--nograph**
    do not generate \*-graph file (pairwise orthology relations)

  - **--subparaCluster**='options'
    additional parameters for the clustering algorithm (proteinortho_clustering) example -subparaCluster='-maxnodes 10000'.
    Note: -rmgraph cannot be set. All other parameters of subparaCluster are replacing the default values (like -cpus or -minSpecies)

  - **--xml**
    do generate an orthologyXML file (see http://www.orthoxml.org for more information). You can also use proteinortho2xml.pl <myproject.proteinortho>.

  - **--core**
  stop clustering if a split would result in groups that do not span across all species of the inital connected component. Overrules the -conn threshold.

  - **--coreMinSpecies**
  sets the minimal number of species for the -core option (default:0)

  - **--coreMaxProts**
  sets the maximal number of proteins per species for the -core option (default:100)

</details>
<br>

 **Misc options**

  - **--checkfasta**
    checks input fasta files if the given algorithm can process the given fasta file.

<details>
  <summary>(Click to expand)</summary>

  - **--cleanblast**
    cleans blast-graph with proteinortho_cleanupblastgraph

  - **--desc**
    write description files (for NCBI FASTA input only)

  - **--binpath**=directory (default: $PATH)
    path to your local executables (blast, diamond, mcl, ...)

  - **--debug**
    gives detailed information for bug tracking

</details>
<br>

 **Large compute jobs**
  - **--jobs**=M/N
    If you want to involve multiple machines or separate a Proteinortho
    run into smaller chunks, use the -jobs=**M**/**N** option. First, run
    'proteinortho6.pl -steps=1 ...' to generate the indices. Then you can
    run 'proteinortho6.pl -steps=2 -jobs=**M**/**N** ...' to run small chunks
    separately. Instead of **M** and **N** numbers must be set representing the
    number of jobs you want to divide the run into (**M**) and the job
    division to be performed by the process. E.g. to divide a Proteinortho
    run into 4 jobs to run on several machines, use 'proteinortho6.pl -steps=2 -jobs=1/4', 'proteinortho6.pl -steps=2 -jobs=1/4', 'proteinortho6.pl -steps=2 -jobs=2/4', 'proteinortho6.pl -steps=2 -jobs=3/4', 'proteinortho6.pl -steps=2 -jobs=4/4'.

    See [Large compute jobs, the --jobs option (proteinortho wiki)](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Large-compute-jobs-(the--jobs-option)) for more details.

<br>

# PoFF

  The PoFF extension allows you to use the relative order of genes (synteny)
  as an additional criterion to disentangle complex co-orthology relations.
  To do so, add the parameter -synteny. You can use it to either come closer
  to one-to-one orthology relations by preferring synthetically conserved
  copies in the presence of two very similar paralogs (default), or just to
  reduce noise in the predictions by detecting multiple copies of genomic
  areas (add the parameter -dups=3). Please note that you need additional
  data to include synteny, namely the gene positions in GFF3 format.
  AsProteinortho is primarily made for proteins, it will only accept GFF
  entries of type CDS (column #3 in the GFF-file). The attributes column
  (#9) must contain Name=GENE IDENTIFIER where GENE IDENTIFIER corresponds
  to the respective identifier in the FASTA format. It may not contain a
  semicolon (;)! Alternatively, you can also set ID=GENE IDENTIFIER. Example
  files are provided in the source code archive. Hence, we can run
  proteinortho6.pl -project=test -synteny test/A1.faa test/B1.faa test/E1.faa
  test/F1.faa to add synteny information to the calculations. Of course,
  this only makes sense if species are sufficiently similar. You won't gain
  much when comparing e.g. bacteria with fungi. When the analysis is done
  you will find an additional file in your current working directory, namely
  test.poff.tsv (tab separated file). This file is equivalent to the test.proteinortho.tsv file (above) but
  can be considered more accurate as synteny was involved for its
  construction.

  A full example is described here: [Synteny example (proteinortho wiki)](https://gitlab.com/paulklemm_PHD/proteinortho/-/wikis/synteny-example).

# Output
 **BLAST Search (step 1-2)**

<details>
  <summary>myproject.blast-graph (Click to expand)</summary>

    filtered raw blast data based on adaptive reciprocal best blast
    matches (= reciprocal best match matches within a range of 95% by default) 
     
    A line starting with # indicates the two species that are analysed below. E.g. '# M.faa L.faa' tells that the next lines are for species M versus species L.
    
    All matches are reciprocal matches. If
    e.g. a match for M_15 L_15 is shown, L_15 M_15 exists implicitly.
    
    E-Values and bit scores for both directions are given behind each
    match. 
    
    The 4 numbers below the species (e.g. '# 3.8e-124        434.9...') are representing the median values for this comparison.

      # file_a    file_b
      # a   b     evalue_ab     bitscore_ab   evalue_ba     bitscore_ba
      # E.faa     C.faa   
      # 3.8e-124        434.9   2.8e-126        442.2
      E_11  C_11  5.9e-51 190.7   5.6e-50 187.61
      E_10  C_10  3.8e-124    434.9   2.8e-126    442.2
      ...
 </details>
 <br>

 **Clustering (step 3)**

<details>
  <summary>myproject.proteinortho-graph (Click to expand)</summary>
    
    clustered version of the myproject.blast-graph.
    
    Its connected components are represented in myproject.proteinortho.tsv / myproject.proteinortho.html.
    
    The format of myproject.blast-graph is the equivalent to the myproject.blast-graph (see above).

      # file_a    file_b
      # a   b     evalue_ab     bitscore_ab   evalue_ba     bitscore_ba
      # E.faa     C.faa
      E_10  C_10  3.8e-124    434.9   2.8e-126    442.2
      E_11  C_11  5.9e-51 190.7   5.6e-50 187.6
      ...
 </details>
 <br>

 <details>
  <summary> myproject.proteinortho.tsv (Click to expand)</summary>

    The connected components of myproject.proteinortho-graph. 
    
    The very first column indicates the number of species covered by this group. 
    The second column indicates the number of genes included in this group. 
    
    If the number of genes is bigger than the number of species, there are co-orthologs present. 
    
    The third column gives the algebraic connectivity of the respective group. This indicates how densely the genes are connected
    in the orthology graph that was used for clustering. 
    A connectivity of 1 indicates a perfect dense cluster with each gene beeing connected/orthologous to each
    other gene. 
    
    By default, Proteinortho splits each group into two more dense subgroups when the connectivity is below 0.1 (can be user defined).
    
    Hint: you can open this file in Excel / Numbers / Open Office.

      # Species   Genes   Alg.-Conn.    C.faa   C2.faa  E.faa   L.faa   M.faa
      2   5     0.16  *     *     *     L_643,L_641   M_649,M_640,M_642
      3   6     0.138   C_164,C_166,C_167,C_2   *     *     L_2   M_2
      2   4     0.489   *     *     *     L_645,L_647   M_644,M_646

 </details>
 <br>

[myproject.proteinortho-graph.summary](https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Tools-and-additional-programs#proteinortho-graphblast-graph-species-summary-table)

 <br>
 <details>
  <summary> myproject.proteinortho.html (Click to expand)</summary>
    The html version of the myproject.proteinortho.tsv file
 </details>
 <br>

 **POFF (-synteny)**

  The synteny based graph files (myproject.ffadj-graph and
  myproject.poff-graph) have two additional columns: same_strand and
  simscore. The first one indicates if two genes from a match are located at
  the same strands (1) or not (-1). The second one is an internal score
  which can be interpreted as a normalized weight ranging from 0 to 1 based
  on the respective e-values. Moreover, a second comment line is followed
  after the species lines, e.g.

    # M.faa L.faa
    # Scores: 4   39    34.000000     39.000000

  <details>
  <summary>myproject.ffadj-graph (Click to expand)</summary>

    filtered blast data based on adaptive reciprocal best blast matches
    and synteny (only if -synteny is set)

 </details>
 <br>

  <details>
  <summary>myproject.poff-graph (Click to expand)</summary>

    clustered ffadj graph. Its connected components are represented in
    myproject.poff.tsv (tab separated file) (only if -synteny is set)

 </details>
 <br>


# EXAMPLES
 **Calling proteinortho**
  Sequences are typically given in plain fasta format like the files in
  test/

  test/C.faa:

    >C_10
    VVLCRYEIGGLAQVLDTQFDMYTNCHKMCSADSQVTYKEAANLTARVTTDRQKEPLTGGY
    HGAKLGFLGCSLLRSRDYGYPEQNFHAKTDLFALPMGDHYCGDEGSGNAYLCDFDNQYGR
    ...

   test/E.faa:

    >E_10
    CVLDNYQIALLRNVLPKLFMTKNFIEGMCGGGGEENYKAMTRATAKSTTDNQNAPLSGGF
    NDGKMGTGCLPSAAKNYKYPENAVSGASNLYALIVGESYCGDENDDKAYLCDVNQYAPNV
    ...

  To run proteinortho for these sequences, simply call

    perl proteinortho6.pl test/C.faa test/E.faa test/L.faa test/M.faa

  To give the outputs the name 'test', call

    perl proteinortho6.pl -project=test test/*faa

  To use blast instead of the default diamond, call

    perl proteinortho6.pl -project=test -p=blastp+ test/*faa

  If installed with make install, you can also call

    proteinortho -project=test -p=blastp+ test/*faa


# Hints
  Using .faa to indicate that your file contains amino acids and .fna to
  show it contains nucleotides makes life much easier.

  Sequence IDs must be unique within a single FASTA file. Consider renaming
  otherwise. Note: Till version 5.15 sequences IDs had to be unique among
  the whole dataset. Proteinortho now keeps track of name and species to
  avoid the necessissity of renaming.

  You need write permissions in the directory of your FASTA files as
  Proteinortho will create blast databases. If this is not the case,
  consider using symbolic links to the FASTA files.

  The directory src contains useful tools, e.g. proteinortho_grab_proteins.pl which
  fetches protein sequences of orthologous groups from Proteinortho output
  table. (These files are installed during 'make install')

# Kmere Heuristic

## Example 1

In the following example a huge blast graph is used for step 3 (clustering).
The first connected component contains 7410694 nodes, hence the kmere heuristic is activated.
Since the fiedler vector would result in a good split, the kmere heuristic is then deactivated immediatly.

<details>
  <summary>as fallback (Click to expand)</summary>

    ...
    [CRITICAL WARNING]   Failed to partition subgraph with 6929 nodes into (6929,0,0) sized groups, now using kmere heuristic as fall-back.
    ...

</details>

<details>
<summary>working example for large graphs (Click to expand)</summary>

    ...
    17:32:15 [DEBUG] (kmere-heuristic) The current connected component is so large that the k-mere heuristic can be used. First: Testing if a normal split would result in a good partition (|.|>20%) of the CC.
     [WARNING] (kmere-heuristic) A normal split would NOT result in a good partition (|.|>20%) of the CC, therefore  the k-mere heuristic is now used. The current connected component will be split in 3.85373 (= number of proteins <6929> / ( n
    odes per species <1> * number of species <1798>)) groups greedily accordingly to the fiedler vector.
    ...

</details>

<details>
<summary>example for large graphs, where kmere is tested but not needed (Click to expand)</summary>

    ...
    20:27:07 [DEBUG] (kmere-heuristic) The current connected component is so large that the k-mere heuristic can be used. First: Testing if a normal split would result in a good partition (|.|>20%) of the CC.
    20:27:09 [DEBUG] (kmere-heuristic) A normal split would result in a good partition (|.|>20%) of the CC, therefore returning now to the normal algorithm (no k-mere heuristic).
    ...

</details>

# Credit where credit is due

 - The all-versus-all BLAST-analysis (-step=2) is only possible with (one of) the following underlying algorithms:
   - NCBI BLAST+ or NCBI BLAST legacy (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
   - Diamond (doi:10.1038/nmeth.3176, https://github.com/bbuchfink/diamond)
   - Last (doi:10.1101/gr.113985.110, http://last.cbrc.jp/)
   - Rapsearch2 (doi:10.1093/bioinformatics/btr595, https://github.com/zhaoyanswill/RAPSearch2)
   - Topaz (doi:10.1186/s12859-018-2290-3, https://github.com/ajm/topaz)
   - usearch,ublast (doi:10.1093/bioinformatics/btq461, https://www.drive5.com/usearch/download.html)
   - blat (http://hgdownload.soe.ucsc.edu/admin/)
   - mmseqs2 (doi:10.1038/nbt.3988 (2017). https://github.com/soedinglab/MMseqs2)
 - The clustering step (-step=3) got a huge speedup with the integration of LAPACK (Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd., http://www.netlib.org/lapack/)
 - The html output of the *proteinortho.tsv (orthology groups) is enhanced by clusterize (https://github.com/NeXTs/Clusterize.js), reducing the scroll lag.

# ONLINE INFORMATION
  For download and online information, see
  <https://www.bioinf.uni-leipzig.de/Software/proteinortho/>
  or
  <https://gitlab.com/paulklemm_PHD/proteinortho>

# REFERENCES
  Lechner, M., Findeisz, S., Steiner, L., Marz, M., Stadler, P. F., &
  Prohaska, S. J. (2011). Proteinortho: detection of (co-) orthologs in
  large-scale analysis. BMC bioinformatics, 12(1), 124.
