#!/bin/bash --login

module load slurm

fullname=$(readlink -f $0)
myname=$(basename $fullname)
if [ -z "$PoPipe_baseDir" ]; then echo "$myname PoPipe_baseDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;
if [ -z "$PoPipe_fastaDir" ]; then echo "$myname PoPipe_fastaDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;
if [ -z "$PoPipe_protOrthoDir" ]; then echo "$myname PoPipe_protOrthoDir not set (did you source poPipeSettings.sh?)"; exit 10; fi;

polog "$myname begin";
mkdir $PoPipe_protOrthoDir;
if [[ ! -d $PoPipe_protOrthoDir ]]; then echo "could not create $PWD/$PoPipe_protOrthoDir"; exit 10; fi;
pushd $PoPipe_protOrthoDir >/dev/null

sname=indices
polog "$sname begin"
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=1 -verbose=1 -checkfasta -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi;
polog "$sname end"

if [[ "$POLargeJobSplit" -eq 1 && "$POLargeJobSplitVarification" == "YES" ]] 
then

if [ "$PONumNodes" -eq 1 ]
then
sname=blaGraph
polog "$sname begin"
polog "Running on a single node (default)"
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=2 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi;
polog "$sname end"

sname=cluster
polog "$sname begin"
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=3 -cleanblast -nograph -debug -clean -verbose=1 -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi;
polog "$sname end"

elif [ "$PONumNodes" -eq 2 ]
then
polog "Multiple Node Run"
sname=blaGraph
polog "$sname begin"
polog "Two nodes"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=1/2 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_1of2.out 2>${std}_1of2.err & 
jobPID1=$!
polog "Job1 PID: $jobPID1"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=2/2 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_2of2.out 2>${std}_2of2.err &
jobPID2=$!
polog "JobB PID: $jobPID2"
wait
polog "$sname end"

sname=cluster
polog "$sname begin"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=3 -debug -cleanblast -nograph -clean -verbose=1 -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err &
clusterJobPID=$!
polog "Cluster Node PID: $clusterJobPID"
wait
polog "$sname end"

elif [ "$PONumNodes" -eq 5 ]
then
polog "Multiple Node Run"
sname=blaGraph
polog "$sname begin"
polog "Five nodes"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=1/5 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_1of5.out 2>${std}_1of5.err & 
jobPID1=$!
polog "Job1 PID: $jobPID1"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=2/5 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_2of5.out 2>${std}_2of5.err & 
jobPID2=$!
polog "Job2 PID: $jobPID2"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=3/5 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_3of5.out 2>${std}_3of5.err & 
jobPID3=$!
polog "Job3 PID: $jobPID3"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=4/5 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_4of5.out 2>${std}_4of5.err & 
jobPID4=$!
polog "Job4 PID: $jobPID4"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=5/5 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_5of5.out 2>${std}_5of5.err & 
jobPID5=$!
polog "Job5 PID: $jobPID5"
wait
polog "$sname end"

sname=cluster
polog "$sname begin"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=3 -debug -cleanblast -nograph -clean -verbose=1 -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err &
clusterJobPID=$!
polog "Cluster Node PID: $clusterJobPID"
wait
polog "$sname end"

elif [ "$PONumNodes" -eq 10 ]
then
polog "Multiple Node Run"
sname=blaGraph
polog "$sname begin"
polog "Ten nodes"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=1/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_1of10.out 2>${std}_1of10.err & 
jobPID1=$!
polog "Job1 PID: $jobPID1"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=2/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_2of10.out 2>${std}_2of10.err & 
jobPID2=$!
polog "Job2 PID: $jobPID2"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=3/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_3of10.out 2>${std}_3of10.err & 
jobPID3=$!
polog "Job3 PID: $jobPID3"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=4/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_4of10.out 2>${std}_4of10.err & 
jobPID4=$!
polog "Job4 PID: $jobPID4"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=5/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_5of10.out 2>${std}_5of10.err & 
jobPID5=$!
polog "Job5 PID: $jobPID5"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=6/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_6of10.out 2>${std}_6of10.err & 
jobPID6=$!
polog "Job6 PID: $jobPID6"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=7/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_7of10.out 2>${std}_7of10.err & 
jobPID7=$!
polog "Job7 PID: $jobPID7"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=8/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_8of10.out 2>${std}_8of10.err & 
jobPID8=$!
polog "Job8 PID: $jobPID8"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=9/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_9of10.out 2>${std}_9of10.err & 
jobPID9=$!
polog "Job9 PID: $jobPID9"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=10/10 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_10of10.out 2>${std}_10of10.err & 
jobPID10=$!
polog "Job10 PID: $jobPID10"

wait

polog "$sname end"

sname=cluster
polog "$sname begin"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=3 -cleanblast -nograph -debug -clean -verbose=1 -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err &
clusterJobPID=$!
polog "Cluster Node PID: $clusterJobPID"

wait

polog "$sname end"

elif [ "$PONumNodes" -eq 20 ]
then
polog "Multiple Node Run"
sname=blaGraph
polog "$sname begin"
polog "Twenty nodes"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=1/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_1of20.out 2>${std}_1of20.err & 
jobPID1=$!
polog "Job1 PID: $jobPID1"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=2/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_2of20.out 2>${std}_2of20.err & 
jobPID2=$!
polog "Job2 PID: $jobPID2"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=3/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_3of20.out 2>${std}_3of20.err & 
jobPID3=$!
polog "Job3 PID: $jobPID3"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=4/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_4of20.out 2>${std}_4of20.err & 
jobPID4=$!
polog "Job4 PID: $jobPID4"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=5/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_5of20.out 2>${std}_5of20.err & 
jobPID5=$!
polog "Job5 PID: $jobPID5"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=6/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_6of20.out 2>${std}_6of20.err & 
jobPID6=$!
polog "Job6 PID: $jobPID6"


srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=7/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_7of20.out 2>${std}_7of20.err & 
jobPID7=$!
polog "Job7 PID: $jobPID7"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=8/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_8of20.out 2>${std}_8of20.err & 
jobPID8=$!
polog "Job8 PID: $jobPID8"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=9/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_9of20.out 2>${std}_9of20.err & 
jobPID9=$!
polog "Job9 PID: $jobPID9"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=10/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_10of20.out 2>${std}_10of20.err & 
jobPID10=$!
polog "Job10 PID: $jobPID10"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=11/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_11of20.out 2>${std}_11of20.err & 
jobPID11=$!
polog "Job11 PID: $jobPID11"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=12/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_12of20.out 2>${std}_12of20.err & 
jobPID12=$!
polog "Job12 PID: $jobPID12"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=13/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_13of20.out 2>${std}_13of20.err & 
jobPID13=$!
polog "Job13 PID: $jobPID13"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=14/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_14of20.out 2>${std}_14of20.err & 
jobPID14=$!
polog "Job14 PID: $jobPID14"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=15/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_15of20.out 2>${std}_15of20.err & 
jobPID15=$!
polog "Job15 PID: $jobPID15"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=16/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_16of20.out 2>${std}_16of20.err & 
jobPID16=$!
polog "Job16 PID: $jobPID16"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=17/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_17of20.out 2>${std}_17of20.err & 
jobPID17=$!
polog "Job17 PID: $jobPID17"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=18/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_18of20.out 2>${std}_18of20.err & 
jobPID18=$!
polog "Job18 PID: $jobPID18"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=19/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_19of20.out 2>${std}_19of20.err & 
jobPID19=$!
polog "Job19 PID: $jobPID19"

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -jobs=20/20 -step=2 $PoPipe_fastaDir/*.fasta 1>${std}_20of20.out 2>${std}_20of20.err & 
jobPID20=$!
polog "Job20 PID: $jobPID20"
wait

polog "$sname end"

sname=cluster
polog "$sname begin"
std=$PoPipe_stdDir/$sname

srun -N 1 -n 1 --mem=$memUsage100 perl /$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -cpus=$cpusPerTask -ram=$memUsage100 -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=3 -cleanblast -nograph -debug -clean -verbose=1 -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err &
clusterJobPID=$!
polog "Cluster Node PID: $clusterJobPID"
wait
polog "$sname end"

else
sname=blaGraph
polog "$sname begin"
polog "Selection of node number error reverting to a single node"
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=2 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi;
polog "$sname end";

sname=cluster
polog "$sname begin"
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=3 -cleanblast -nograph -debug -clean -verbose=1 -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err &
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi;
polog "$sname end"
fi

#Not a 'large job': default
else
sname=blaGraph
polog "$sname begin"
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=2 -verbose=1 -cleanblast -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi;
polog "$sname end"

sname=cluster
polog "$sname begin"
std=$PoPipe_stdDir/$sname
$PoPipe_srcBaseDir/externalSoftware/proteinOrtho/$PoPipe_proteinOrthoVersion/proteinortho6.pl -binpath=$PoPipe_srcBaseDir/externalSoftware/blast/bin -binpath=$PoPipe_srcBaseDir/externalSoftware/diamond -step=3 -cleanblast -nograph -debug -clean -verbose=1 -conn=$PoPipe_proteinOrthoConn $PoPipe_fastaDir/*.fasta 1>${std}.out 2>${std}.err
if [ $? -ne 0 ]; then polog "$sname error $?"; exit $?; fi;
polog "$sname end"
fi


popd >/dev/null

polog "$myname end"
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then return 0; else exit 0; fi;