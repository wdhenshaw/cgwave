eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to run WaveHoltz cases
#
#    POLLUTION CONVERGENCE 
#
#  usage:
#   runsCases.p -grid=square -freq=10 -beta=3 -numRes=1
#
# @caseNames = @ARGV;
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
#

sub factorial {
    my $num = shift;
    my $result = 1;
    for my $i (1 .. $num) {
        $result *= $i;
    }
    return $result;
}

sub nchoosek {
    my ($n, $k) = @_;
    if ($k < 0 || $k > $n) {
        return 0; # Invalid input
    }
    if ($k == 0 || $k == $n) {
        return 1; # Base cases
    }
    if ($k > $n / 2) {
        $k = $n - $k; # Optimization: C(n,k) = C(n, n-k)
    }
    return factorial($n) / (factorial($k) * factorial($n - $k));
}

$runCodes = 1; # set to zero to not re-run the codes

$grid = "square";
$known = "modulatedGaussian"; 
# $known = "boxHelmholtz"; 
$freq=40; $beta=1; $x0=.55; $y0=.55; 

$order =2; 
$numRes=3;
$resStart=1;
foreach $arg ( @ARGV )
{
  if( $arg =~ /-grid=(.*)/ )
  {
    $grid=$1;
  }
  elsif( $arg =~ /-order=(.*)/ ) {  $order=$1; } 
  elsif( $arg =~ /-numRes=(.*)/ ) { $numRes=$1; } 
  elsif( $arg =~ /-resStart=(.*)/ ){ $resStart=$1; }   
  elsif( $arg =~ /-freq=(.*)/ ) {  $freq=$1; }   
  elsif( $arg =~ /-beta=(.*)/ ) {  $beta=$1; }   
  elsif( $arg =~ /-x0=(.*)/ ) {  $x0=$1; }   
  elsif( $arg =~ /-y0=(.*)/ ) {  $y0=$1; }   
  elsif( $arg =~ /-runCodes=(.*)/ ){ $runCodes=$1; }             
  elsif( $arg =~ /-known=(.*)/ ){ $known=$1; }             
}

printf("grid=[$grid]\n");

# @krylovTypes = ("gmres","bicgstab");
# @krylovTypes = ("gmres", "bicgstab" );
$cgwh = "/home/henshw/Dropbox/research/cgwave/bin/cgwh";
$comp = "/home/henshw/Overture.g/bin/comp";


# Modulated Gaussian
if( $known eq "modulatedGaussian")
{
  $CASE = "modGaussFreq$freq";
}
else
{
  $CASE = "trigFreq$freq";
}

$c=1;
$amp=1;  # Solution amplitude
$L = 1.; # domain size
if( $grid =~ "rectangle4By1" ){ $L = 4; }
if( $grid =~ "rectangle8By1" ){ $L = 8; }
if( $grid =~ "rectangle16By1" ){ $L = 16; }

if( $grid =~ "sic" ){ $L = 2; }  # disk


#choose kx^2 + ky^2 = k^2 
$kx = - $freq/sqrt(2.);
$ky = $kx;
$cmd = "setsid $cgwh -noplot trigHelmholtz.cmd -known=$known -g=GRID -nf=1 -freq=$freq -k0=-$freq -beta=$beta -x0=$x0 -y0=$y0 -kx=$kx -ky=$ky -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -matlab=MATLAB -show=SHOW -go=ds > OUTPUTFILE";

$g=1; 
if( $grid eq "square" || $grid eq "nonSquare" )
{
  $ds0 = 1./32.; # initial grid space
  $grids[$g]= $grid ."32.order"; $g++;
  $grids[$g]= $grid ."64.order"; $g++;
  $grids[$g]= $grid ."128.order"; $g++;
  $grids[$g]= $grid ."256.order"; $g++;
  $grids[$g]= $grid ."512.order"; $g++;
  $grids[$g]= $grid ."1024.order"; $g++;
  $grids[$g]= $grid ."2048.order"; $g++;
}
elsif( $grid eq "box" || $grid eq "nonBox"  )
{
  $ds0= 1./10.; 
  $grids[$g]= $grid ."1.order"; $g++;
  $grids[$g]= $grid ."2.order"; $g++;
  $grids[$g]= $grid ."4.order"; $g++; 
  $grids[$g]= $grid ."8.order"; $g++;  
  $grids[$g]= $grid ."16.order"; $g++;      
  $grids[$g]= $grid ."32.order"; $g++;      
  $grids[$g]= $grid ."64.order"; $g++;      
  $grids[$g]= $grid ."128.order"; $g++;      
  $grids[$g]= $grid ."256.order"; $g++;      
}
else
{
  # DEFAULT GRID NAMES 

  # $ds0= 1./10.; 
  # $grids[$g]= $grid ."1.order"; $g++;
  # $grids[$g]= $grid ."2.order"; $g++;
  $ds0= 1./40.; 
  $grids[$g]= $grid ."4.order"; $g++; 
  $grids[$g]= $grid ."8.order"; $g++;  
  $grids[$g]= $grid ."16.order"; $g++;      
  $grids[$g]= $grid ."32.order"; $g++;      
  $grids[$g]= $grid ."64.order"; $g++;      
  $grids[$g]= $grid ."128.order"; $g++;      
  $grids[$g]= $grid ."256.order"; $g++;      

}


$orderStart=$order; $orderEnd=$order; 
for( $order=$orderStart; $order<=$orderEnd; $order=$order+2 )
{
  printf("\n =================== ORDER = $order ==========================\n");

  $i=0;

  for( $g=$resStart; $g<=$numRes; $g++ )
  {
    $myGrid = $grids[$g] . $order; 
    $myGridName = $myGrid;
    $myGridName =~ s/\.//g; # remove "."
    printf("  ----------- GRID = $myGrid (g=$g, resStart=$resStart, numRes=$numRes) --------------\n");

    $myCmd = $cmd;
    $myCmd =~ s/MATLAB/$CASE$myGridName/g;
    $myCmd =~ s/GRID/$myGrid/g;
    $myCmd =~ s/CASE/$CASE/g;

    $showFileName = "$CASE$myGridName.show";
    $myCmd =~ s/SHOW/$showFileName/g;

    $showFile[$i]=$showFileName; # save for later

    $outputFile = "$CASE$myGrid.out";
    $myCmd =~ s/OUTPUTFILE/$outputFile/g;

    if( $runCodes )
    {
      printf("Run [$myCmd]\n");
      # Execute the command: 
      my $output = `$myCmd`; 
    }

    # EXTRACT ERRORS:
    $fileName = $outputFile;
    open(FILE,"$fileName") || die "cannot open file $fileName!" ;
    while( <FILE> )
    {
      $line = $_;
      chop($line);

      $subString = "CgWave::getErrors:";
      if( $line =~ /$subString/ )
      {
        # printf("STRING FOUND\n");
        # printf("line=[%s]\n",$line);

        $l2RelErr = $line;
        $l2RelErr =~ s/.*l2RelErr=(.*)/\1/;
        # printf("grid=$myGrid, l2RelErr=$l2RelErr\n");

        $maxRelErr = $line;
        $maxRelErr =~ s/.*maxRelErr=(\S+)/\1/;
        printf("grid=$myGrid, l2RelErr=$l2RelErr, maxRelErr=$maxRelErr\n");

        $err[$i]=$l2RelErr;

        $maxErr[$i] = $maxRelErr;

        break;

      }
      else
      {
        # printf("NOT FOUND line=[%s]\n",$line);
      }


    }
    close(FILE);


    $i = $i+1;


  } # end for res

  $num = $i; 
  $pi = atan2(1.,1.)*4.;
  # printf("pi=$pi\n");




  # --------- RUN COMP TO PERFORM RICHARDSON EXTRAP ----------
  $iw = $order+1; # interp width 
  $outFile = "comp$CASE.cmd";
  open(COMPFILE,">$outFile") || die "cannot open output file $outFile!" ;

  $compLogFile = "comp$CASE.log";
  print COMPFILE "#\n";
  print COMPFILE "# comp command file written by runCases.p\n";
  print COMPFILE "output file name: $compLogFile\n";
  print COMPFILE "interpolation width: $iw\n";
  print COMPFILE "specify files (coarse to fine)\n";
  for( $i=0; $i<$num; $i++ )
  {
    print COMPFILE "$showFile[$i]\n";
  }
  print COMPFILE "exit\n";
  print COMPFILE "choose a solution\n";
  print COMPFILE "1\n";
  print COMPFILE "compute errors\n";
  print COMPFILE "exit\n";

  close(COMPFILE);
  printf("Wrote comp file [$outFile]\n");

  if( $runCodes )
  {
    sleep(1); # wait for comp file to be closed

    $myCmd = "$comp -noplot $outFile";
    printf("Run [$myCmd]\n");
    # Execute the command: 
    my $output = `$myCmd`; 

    sleep(1); # wait for compLogFile file to be closed
  }

  printf("READ $compLogFile...\n");
  open(FILE,"$compLogFile") || die "cannot open file $compLogFile!" ;
  $found=0; 
  while( <FILE> )
  {
    $line = $_; chop($line);
    # printf("line=[%s]\n",$line);

    $subString = "l2 norm results";
    if( $line =~ /$subString/ )
    {
      $found=1; 
      printf("FOUND line=[%s]\n",$line);
      while( <FILE> )
      {
        $line = $_; chop($line);
        printf("line=[%s]\n",$line);
        $ss="uDHS0";
        if( $line =~ /$ss/ )
        {
          # printf("ERROR LINE=[%s]\n",$line);
          for( $i=0; $i<$num; $i++ )
          {
             $compErr = $line;
             $compErr =~ s/.* ee\($i\) = (\S+)\, .*/\1/;
             $errRE[$i] = $compErr; 
             printf("errRE[$i]=$errRE[$i]\n");
          }
          last;
        }
        if( $line =~ /$subString/ )
        {
          printf("\n DONE READING comp.log\n");
          last;
        }
      }
    }
    if( $found eq "1" ){ last; }
  }
  close(FILE);
  # ------------- END COMP ------------------


  printf("\n ====================================\n");
  printf("     SUMMARY: known=$known, grid=$myGrid, L=$L\n");
  for( $i=0; $i<$num; $i++ )
  {
    printf(" i=$i: l2RelErr=$err[$i]");

    $k = $freq/$c; 
    $lambda=2*$pi/$k;
    $Nlambda = $L/$lambda;
    $ds=$ds0/2**($i+$resStart-1);
    $ppw[$i] = $lambda/$ds;

    $p = $order;
    $mu=$p/2; 
    $bp = 2./( ($mu+1)**2 * nchoosek(2*$mu+2,$mu+1) );
    $errPPW[$i] = .5*$bp*$k*$L*( $k*$ds )**$p;

    $theoryRatio = $errPPW[$i]/$err[$i];

    if( $i>0 )
    { 
      $ratio = $err[$i-1]/$err[$i]; 
      $rate = log($ratio)/log(2);
      printf(" (theory=%9.2e, theory/comp=%5.2f) ratio=%5.2f, rate=%5.2f, PPW=%6.1f (ds=%9.2e, lambda=%5.2f, Nlambda=%5.1f)",
             $errPPW[$i],$theoryRatio,$ratio,$rate,$ppw[$i],$ds,$lambda,$Nlambda); 
    }

    printf("\n");
  }  

  printf("Format for plotErrorVersusPPW.m\n");

   #  order{m}=4; 
   #  L{m}=4; 
   #  omega{m}=25; 
   #  Nv{m} = [4,8,16,32];
   #  amp{m} = 2.6;  % amplitiude of solution
   #  errv{m} = [8.6e-1, 5.2e-2, 3.2e-3, 2.0e-4 ]./amp{m}; % L2 norm errors from Richardson extrapolation  

  printf("    m=m+1; \n");
  printf("    %% known = $known;\n");
  printf("    order{m}=$order;\n");
  printf("    omega{m}=$freq;\n");
  printf("    amp{m}=$amp;\n");
  printf("    L{m}=$L;\n");
  printf("    x0{m}=$x0; y0{m}=$y0;\n");
  printf("    dsv{m}    = [");
  for( $i=0; $i<$num; $i++ ){ $ds=$ds0/2**($i+$resStart-1); printf("%12.4e ",$ds); } printf("];\n");

  printf("    errv{m}   = [");
  for( $i=0; $i<$num; $i++ ){ printf("%12.4e ",$err[$i]); } printf("];\n");

  printf("    errPPW{m} = [");
  for( $i=0; $i<$num; $i++ ){ printf("%12.4e ",$errPPW[$i]); } printf("];\n");


  printf("    errRE{m}  = [");
  for( $i=0; $i<$num; $i++ ){ printf("%12.4e ",$errRE[$i]); } printf("];\n");

  printf("    maxErrv{m}= [");
  for( $i=0; $i<$num; $i++ ){ printf("%12.4e ",$maxErr[$i]); } printf("];\n");  

}
