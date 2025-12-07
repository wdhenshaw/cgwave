eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to create a SLURM batch script for running on the DCS cluster at the CCI of RPI
#
# Usage:
#  makeBatch.p -batchFile=<script-name> -solverName=<s> -np=<i> -nodes=<i> -gpus=<i> -timeLimit=hh:mm:ss <commands to run>
# 
# Examples:
# cgWave: 
#  makeBatch.p -batchFile=batch.sh -solverName=cgWave $cgWave gaussian.cmd -noplot -g=square256.order2 -tp=0.1 -tf=1.0 -beta=15 -x0=.45 -y0=.55 -z0=0 -omega=15 -amp=-225 -ts=explicit -cfl=.95 -upwind=0 -orderInTime=-1 -bc=d -debug=0 -show=gaussianOmega15.show -go=go
# 
# cgwh: (WaveHoltz)
#  makeBatch.p -batchFile=batch.sh -solverName=cgwh $cgwh -noplot gaussian.cmd -g=square64.order2 -x0=0.6 0.4 .55 .5 -y0=0.45 0.55 .5 .5 -nf=3 -freq=5.1 10.1 15.1 -amp=25 100 225 -beta=15 15 15 -adjustOmega=1 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=1 -matlab=square128OptFilter -useOptFilter=1 -solveri=best -rtoli=1e-12 -atoli=1e-12 -solverh=best -restart=2000 -iluh=20 -maxith=10000 -rtolh=1.0e-12 -atolh=1.0e-12 -show=square128OptFilter.show -go=dfks
#
# EigenWave:
#  makeBatch.p -batchFile=batch.sh -solverName=cgwh $cgwh noplot eig.cmd -g=square128.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -show=squareEigs.show -orderOfExtrapolation=-1 -freq=10 -numEigs=10  -matlab=square  -solveri=best -rtoli=1e-12 -atoli=1e-12  -go=ks 
# 

printf("Create a SLURM batch script. Usage:\n");
printf("makeBatch.p -batchFile=<script-name> -solverName=<s> -np=<i> -nodes=<i> -gpus=<i> -timeLimit=hh:mm:ss <commands to run>\n");

# $cgWave = "/gpfs/u/home/PCM2/PCM2hnsh/barn/cgwave/bin/cgWave";

# defaults:

$solverName = "cgWave";
$batchFile = "batch.sh";


$timeLimit = "00:05:00";  # time-limit in hh:mm:ss
$gpus = 1;  # request this many gpus per node

$np=4;    # total number of MPI tasks
$nodes=1; # number of nodes 

$cmd = "";   # string of commands
$grid = "";

foreach $arg ( @ARGV )
{
  if(    $arg =~ /-batchFile=(.*)/ ){ $batchFile=$1; }
  elsif( $arg =~ /-solverName=(.*)/ ){ $solverName=$1; }  
  elsif( $arg =~ /-np=(.*)/ ){  $np=$1; }
  elsif( $arg =~ /-nodes=(.*)/ ){ $nodes=$1; }  
  elsif( $arg =~ /-gpus=(.*)/ ){ $gpus=$1; }   
  elsif( $arg =~ /-timeLimit=(.*)/ ){ $timeLimit=$1;}   
  else
  {
    # add the arg to the main command line 
    if( $arg =~ /-g=(.*)/ ){ $grid=$1; }  

    $cmd .= " $arg"; # add this arg to the command
  }
}

$jobName = $solverName;

printf("Setting: solverName=$solverName, np=$np, nodes=$nodes, timeLimit=$timeLimit, gpus=$gpus (gpus per node)\n");
printf("cmd=[$cmd]\n");

$logFile = "output/$jobName" . "Grid$grid" . "Np$np.log";
$errFile = "output/$jobName" . ".err";

open($outfile, '>', $batchFile) or die "Cannot open $batchFile for writing: $!";

print $outfile "#!/bin/bash -xe\n";
print $outfile "#\n";
print $outfile "# Usage: (over-ride options on the command line)\n";
print $outfile "#    sbatch --ntasks=8 --nodes=1 batch.sh\n";
print $outfile "#\n";
print $outfile "#SBATCH --job-name=$jobName              # Job name\n";
print $outfile "#SBATCH --gres=gpu:$gpus                   # Always ask for a GPU on DCS or the job wll be refused, set to 6 to avoid node sharing\n";
print $outfile "#SBATCH --ntasks=$np                     # Number of MPI tasks (i.e. processes)\n";
print $outfile "#SBATCH --nodes=$nodes                      # Maximum number of nodes to be allocated\n";
print $outfile "#SBATCH --time=$timeLimit                # Wall time limit (days-hrs:min:sec)\n";
print $outfile "#SBATCH --output=$logFile  # log file (stdout)\n";
print $outfile "#SBATCH --error=$errFile   # error file (stderr)\n";
print $outfile "\n";
print $outfile "echo \"Date              = \$(date)\"\n";
print $outfile "echo \"Hostname          = \$(hostname -s)\"\n";
print $outfile "echo \"Working Directory = \$(pwd)\"\n";
print $outfile "echo \"\"\n";
print $outfile "echo \"Number of Nodes Allocated      = \$SLURM_JOB_NUM_NODES\"\n";
print $outfile "echo \"Number of Tasks Allocated      = \$SLURM_NTASKS\"\n";
print $outfile "echo \"Number of Cores/Task Allocated = \$SLURM_CPUS_PER_TASK\"\n";
print $outfile "\n";
# Run mod p in tcsh to set the environment
print $outfile "tcsh\n";
# print $outfile "mod p\n";
# 
print $outfile "mpirun $cmd";
# print $outfile "mpirun $cgWave gaussian.cmd -noplot -numParGhost=$numParGhost -g=$grid  -tp=$tp -tf=$tf -beta=$beta -x0=$x0 -y0=$y0 -z0=$z0 -omega=$omega -amp=$amp -ts=explicit -cfl=.95 -upwind=1 -orderInTime=-1 -bc=d -debug=0 -implicitUpwind=0 -show=$show -go=go\n";

print $outfile "\n";

close($outfile);

printf("Wrote batch script = [$batchFile]\n");
printf("To run and over-ride SBATCH options on the command line:\n");
printf("sbatch [--ntasks=8] [--nodes=1] $batchFile\n");
printf("To run:\n");
printf("sbatch $batchFile\n");


exit(0);


