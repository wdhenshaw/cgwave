eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to submit a batch job
# 
# Examples:
#  submit.p -grid=square64.order2 -tf=1 -tp=0.1 -show=gaussian.show -np=4 -timeLimit=00:05:00
#  submit.p -grid=square256.order4.ng3 -tf=1 -tp=0.1 -x0=.55 -y0=0.55 -omega=20 -beta=20 -amp=-400 -show=gaussian.show -np=4 -timeLimit=00:05:00
#  submit.p -grid=square128.order4.ng3 -numParGhost=3 -tf=1 -tp=0.1 -x0=.55 -y0=0.55 -omega=20 -beta=20 -amp=-400 -show=gaussian.show -np=4 -timeLimit=00:05:00

printf("Submit a batch job. Usage:\n");
printf("submit.p -jobName=<s> -grid=<s> -tf=<f> -debug=%d -np=<i> -nodes=<i> -timeLimit=hh:mm:ss -gpus=<i>\n");

$cgWave = "/gpfs/u/home/PCM2/PCM2hnsh/barn/cgwave/bin/cgWave";

# defaults:
$jobName = "cgWave";
$grid = "square64.order2";

$tf=1.0;
$tp=0.1;
$debug=0; 
$show=" ";
$x0=0.5; $y0=0.5; $z0=0.;
$omega=10; $beta=10; $amp=-200; 
$numParGhost=-1; # -1 = use default (automatically guessed)

$timeLimit = "00:05:00";  # time-limit in hh:mm:ss
$gpus = 1;  # request this many gpus per node

$np=4;    # total number of MPI tasks
$nodes=1; # number of nodes 

foreach $arg ( @ARGV )
{
  if( $arg =~ /-grid=(.*)/ ){ $grid=$1; }
  elsif( $arg =~ /-tf=(.*)/ ){ $tf=$1; }
  elsif( $arg =~ /-tp=(.*)/ ){ $tp=$1; }
  elsif( $arg =~ /-x0=(.*)/ ){ $x0=$1; }
  elsif( $arg =~ /-y0=(.*)/ ){ $y0=$1; }
  elsif( $arg =~ /-z0=(.*)/ ){ $z0=$1; }
  elsif( $arg =~ /-omega=(.*)/ ){ $omega=$1; }
  elsif( $arg =~ /-beta=(.*)/ ){ $beta=$1; }
  elsif( $arg =~ /-amp=(.*)/ ){ $amp=$1; }
  elsif( $arg =~ /-tp=(.*)/ ){ $tp=$1; }
  elsif( $arg =~ /-numParGhost=(.*)/ ){ $numParGhost=$1; }  
  elsif( $arg =~ /-np=(.*)/ ){  $np=$1; }
  elsif( $arg =~ /-nodes=(.*)/ ){ $nodes=$1; }  
  elsif( $arg =~ /-debug=(.*)/ ){ $debug=$1; }  
  elsif( $arg =~ /-show=(.*)/ ){ $show=$1; }  
  elsif( $arg =~ /-gpus=(.*)/ ){ $gpus=$1; }   
  elsif( $arg =~ /-timeLimit=(.*)/ ){ $timeLimit=$1;}   
  else
  {
    printf("Unknown arg=$rag\n");
  }
}

printf("Setting: tf=$tf, tp=$tp, grid=$grid, debug=$debug\n");
printf("         jobName=$jobName, np=$np, nodes=$nodes, timeLimit=$timeLimit, gpus=$gpus (gpus per node)\n");

$batchFile = "batch.sh";

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
print $outfile "mpirun $cgWave gaussian.cmd -noplot -numParGhost=$numParGhost -g=$grid  -tp=$tp -tf=$tf -beta=$beta -x0=$x0 -y0=$y0 -z0=$z0 -omega=$omega -amp=$amp -ts=explicit -cfl=.95 -upwind=1 -orderInTime=-1 -bc=d -debug=0 -implicitUpwind=0 -show=$show -go=go\n";

print $outfile "\n";

close($outfile);

printf("Wrote batch script = [$batchFile]\n");
printf("To run and over-ride SBATCH options on the command line:\n");
printf("sbatch [--ntasks=8] [--nodes=1] batch.sh\n");
printf("To run:\n");
printf("sbatch batch.sh\n");


exit(0);

