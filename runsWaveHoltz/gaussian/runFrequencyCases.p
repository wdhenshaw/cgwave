eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to run WaveHoltz cases
#
#   RUN CASES TO SHOW ITERAION COUNTS AS A FUNCTION OF FREQUENCY
#
#  usage: 
#         runFrequencyCases.p -startCase=<i> -endCase=<i> -testRun=[0|1] -parallel=[0|1] -bc=[rrrr|drdr|rrdd] -go=[f|k|fk]
# 
# Wave guide: 
#   runFrequencyCases.p -bc=rrdd -startCase=0 -endCase=10 -maxIterations=20000 -testRun=1
#

$testRun=0; # set to 1 for a test run with no big computations
$startCase=0;
$endCase=10000; 
$parallel=0;    # set to 1 to run the job in parallel 
$show=0;        # set to 1 to save a show file 
$bc = "rrrr"; 
$maxIterations=1000; 
$go = "fk";

foreach $arg ( @ARGV )
{
  if( $arg =~ /-startCase=(.*)/        ){ $startCase=$1;     }
  elsif( $arg =~ /-endCase=(.*)/       ){ $endCase=$1;       }
  elsif( $arg =~ /-testRun=(.*)/       ){ $testRun=$1;       } 
  elsif( $arg =~ /-parallel=(.*)/      ){ $parallel=$1;      } 
  elsif( $arg =~ /-show=(.*)/          ){ $show=$1;          }  
  elsif( $arg =~ /-bc=(.*)/            ){ $bc=$1;            }        
  elsif( $arg =~ /-go=(.*)/            ){ $go=$1;            }        
  elsif( $arg =~ /-maxIterations=(.*)/ ){ $maxIterations=$1; }        
  else
  {
    printf("Unknown arg=[$arg]\n");
  }
}



$cwh  = "/home/henshw/Dropbox/complexWaveHoltz";
if( $parallel ==0 )
{
  $cgwh = "/home/henshw/Dropbox/research/cgwave/bin/cgwh";
}
else
{
  $cgwh = "/home/henshw/Dropbox/research/cgwavep/bin/cgwh";
}

sub nChooseK {
    my ($n, $k) = @_;
    return 0 if $k > $n || $k < 0;
    return 1 if $k == 0 || $k == $n;
    $k = $n - $k if $k > $n / 2; # Symmetry property

    my $result = 1;
    for (1 .. $k) {
        $result = $result * ($n - $_ + 1) / $_;
    }
    return $result;
}

sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }


$c=1;
$domainSize=1;
$pi = atan2(1.,1.)*4;
$epsr=1e-2; 



@gridName = ();
@freqArray=();
@ds = ();

$i=0;
$gridName[$i] = "square128";   $freqArray[$i]=30;   $ds[$i]=1./128;  $np[$i]= 1; $i=$i+1;
$gridName[$i] = "square128";   $freqArray[$i]=50;   $ds[$i]=1./128;  $np[$i]= 1; $i=$i+1;
$gridName[$i] = "square256";   $freqArray[$i]=80;   $ds[$i]=1./256;  $np[$i]= 1; $i=$i+1;
$gridName[$i] = "square512";   $freqArray[$i]=120;  $ds[$i]=1./512;  $np[$i]= 1; $i=$i+1;
$gridName[$i] = "square512";   $freqArray[$i]=160;  $ds[$i]=1./512;  $np[$i]= 1; $i=$i+1;
$gridName[$i] = "square1024";  $freqArray[$i]=250;  $ds[$i]=1./1024; $np[$i]= 2; $i=$i+1;
$gridName[$i] = "square1024";  $freqArray[$i]=300;  $ds[$i]=1./1024; $np[$i]= 2; $i=$i+1;
$gridName[$i] = "square2048";  $freqArray[$i]=400;  $ds[$i]=1./2048; $np[$i]= 4; $i=$i+1;
$gridName[$i] = "square2048";  $freqArray[$i]=500;  $ds[$i]=1./2048; $np[$i]= 4; $i=$i+1;

$gridName[$i] = "square4096";  $freqArray[$i]=700;  $ds[$i]=1./4096; $np[$i]= 8; $i=$i+1;
$gridName[$i] = "square4096";  $freqArray[$i]=900;  $ds[$i]=1./4096; $np[$i]= 8; $i=$i+1;

$gridName[$i] = "square8192";  $freqArray[$i]=1200; $ds[$i]=1./8192; $np[$i]=16; $i=$i+1;
$gridName[$i] = "square8192";  $freqArray[$i]=1500; $ds[$i]=1./8192; $np[$i]=16; $i=$i+1;

$orderOfAccuracy=4;

if( $show==1 )
{  # save a show file 
  $go = $go . "s";
}

$numCases=$i;

$endCase = min($endCase,$numCases-1);

for( $i=$startCase; $i<=$endCase; $i=$i+1 )  
{
  $grid = $gridName[$i];
  $freq = $freqArray[$i];
  $amp = $freq*$freq;

  $matlab = "$grid" ."O4" . "B$bc" . "Freq$freq";

  $outFile = $matlab . ".out"; # pipe output from run here 
  $dx = $ds[$i];

  #  ----- Compute points-per-wavelength --- 
  $kWaveNumber = $freq/$c;                 # k = omega/c
  $lambdaWaveLength = 2*$pi/$kWaveNumber;   # lambda = 2*pi/k
  $ppw[$i] = $lambdaWaveLength/$dx;             # PPW = lambda/dx 
  $NLambda[$i] = $domainSize/$lambdaWaveLength; # size of domain in wavelengths
  $mu = $orderOfAccuracy/2; 
  $bp = 2./( ($mu+1)*($mu+1) * nChooseK(2*$mu+2,$mu+1) );
  $ppwRuleOfThumb[$i] = 2.*$pi*( $pi*$bp*$NLambda[$i]/$epsr)**(1./$orderOfAccuracy);
  # printf(" ppwRuleOfThumb=$ppwRuleOfThumb\n"); 

  # ---------------------------------------- 

  # printf("nChooseK(%d,%d)=%d, pi=$pi\n",2*$mu+2,$mu+1,nChooseK(2*$mu+2,$mu+1));

  printf("\n Case $i:\n");
  printf(" omega=%10.3e k=%8.2e lambda=%8.2e dx=%8.2e domainSize=%8.2e NLambda=%g ppw=lambda/dx=%5.1f ppw(rule-of-thumb)=%5.1f (epsr=%g)\n",
           $freq,$kWaveNumber,$lambdaWaveLength,$dx,$domainSize,$NLambda[$i],$ppw[$i],$ppwRuleOfThumb[$i],$epsr); 

  $prefix = "";
  if( ($parallel ==1) && ($np[$i] > 1) )
  { # run this job in parallel
    $prefix = "mpirun -np $np[$i] ";
  }
  

  $showFile = "";
  if( $show==1 )
  { # save a show file 
    $showFile="-show=$matlab.show";
  }

  if( $bc eq "rrrr" )
  { # open domain 
    $bcOption = "-bc=a";
  }
  elsif( $bc eq "rrdd" )
  { # wave guide
    $bcOption = "-bc=a -bc3=d -bc4=d";
  }
  elsif( $bc eq "drdr" )
  { # quarter plane 
    $bcOption = "-bc=a -bc1=d -bc3=d";
  }
  else
  {
    printf("ERROR: Unknown bc=[$bc]\n");
    exit(1);
  }

  $cmd = "$prefix$cgwh -noplot gaussian.cmd -g=$grid.order4 -ts=explicit -cfl=.95 -x0=0.55 -y0=0.6 -nf=1 -freq=$freq -amp=$amp -beta=$freq -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=$maxIterations -numPeriods=1 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 $bcOption -implicitUpwind=0 -takeImplicitFirstStep=1 -filterD0t=1 -matlab=$matlab $showFile -go=$go > $outFile";

  printf("Run [$cmd]\n");

  if( $testRun == 0 )
  {
      $rt = system($cmd);
      printf("Done: rt=$rt\n");
      # Use this next option to capture the output from the run
      # my $output = `$cmd`;

      $cpCmd = "cp $matlab*.m $cwh/matlab/results"; 
      printf("Copy files: [$cpCmd]\n");
      $rt = system($cpCmd);
  }

}

printf("\n ===== SUMMARY OF PPW ======\n");
for( $i=$startCase; $i<=$endCase; $i=$i+1 )  
{
  printf("case %2d : omega=%10.3e domainSize=%10.2e Nlambda=%10.2e, ppw=lambda/dx=%5.1f ppw(rule-of-thumb)=%5.1f (epsr=%g) grid=$gridName[$i]\n",
           $i,$freqArray[$i],$domainSize,$NLambda[$i],$ppw[$i],$ppwRuleOfThumb[$i],$epsr);  
}
