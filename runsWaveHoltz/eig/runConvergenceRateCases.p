eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to run WaveHoltz cases
#
#    EIGENWAVE: RUN CASES TO COMPARE CONVERGENCE RATES TO THEORY 
#
#  usage: 
#         runConvergentRatesCases.p caseNames
# 
# runConvergentRateCases.p -grid=square 
#
# @caseNames = @ARGV;
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
#

$grid = "square";
$order =2; 
foreach $arg ( @ARGV )
{
  if( $arg =~ /-grid=(.*)/ )
  {
    $grid=$1;
  }
  elsif( $arg =~ /-order=(.*)/ )
  {
    $order=$1;
  } 
}

# printf("grid=[$grid]\n");

# @krylovTypes = ("gmres","bicgstab");
# @krylovTypes = ("gmres", "bicgstab" );
$cgwh = "/home/henshw/Dropbox/research/cgwave/bin/cgwh";



# square128 freq=40
$cmd = "setsid $cgwh -noplot eig.cmd -g=square128.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=/home/henshw/runs/eig/square128O2EigsEv512.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=NP -show=junk.show  -orderOfExtrapolation=-1 -freq=40 -numEigs=NR  -matlab=square128O2Freq40NrNRNpNP -go=k > square128O2Freq40NrNRNpNP.out";

# square256 freq=60 
$cmd = "setsid $cgwh -noplot eig.cmd -g=square256.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=100 -eigenVectorFile=/home/henshw/runs/eig/square256O2Ev512.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=NP -show=junk.show  -orderOfExtrapolation=-1 -freq=60 -numEigs=NR  -matlab=square256O2Freq60NrNRNpNP -go=k > square256O2Freq60NrNRNpNP.out";

# $grid = "square128.order2";


@npv = (2,4,6,8,10,12,14,16);  # Np : number of periods
# @nrv = (5,10,15,20);  # number of requested eigenvalues 
@nrv = (4,6,8,10,12,14,16,18,20);  # number of requested eigenvalues 
for $np (@npv) 
{
  for $nr (@nrv)  
  {
    printf("\n =================== np = $np, nr=$nr ==========================\n");

     $myCmd = $cmd;
     $myCmd =~ s/NR/$nr/g;
     $myCmd =~ s/NP/$np/g;

     printf("Run [$myCmd]\n");
     # Execute the command: 
     my $output = `$myCmd`; 

  }
}