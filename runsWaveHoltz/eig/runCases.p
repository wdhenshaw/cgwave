eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to run WaveHoltz cases
#
#    EIGENWAVE: COMPARE EXPLICIT TO IMPLICIT TIME STEPPING
#
#  usage: 
#         runCases.p caseNames
# 
# runCases.p -grid=square -order=2 -numRes=4
# runCases.p -grid=square -order=4 -numRes=4
# runCases.p -grid=nonSquare -order=2 -numRes=3
#
# @caseNames = @ARGV;
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
#

$grid = "square";
$order =2; 
$numRes=3;
$resStart=1;
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
  elsif( $arg =~ /-numRes=(.*)/ )
  {
    $numRes=$1;
  } 
  elsif( $arg =~ /-resStart=(.*)/ )
  {
    $resStart=$1;
  }       
}

printf("grid=[$grid]\n");

# @krylovTypes = ("gmres","bicgstab");
# @krylovTypes = ("gmres", "bicgstab" );
$cgwh = "/home/henshw/Dropbox/research/cgwave/bin/cgwh";


# Explicit time-stepping:
$cmd[0] = "setsid $cgwh -noplot eig.cmd -g=GRID -nf=1 -adjustOmega=0 -ts=explicit -cfl=0.95 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16 -matlab=MATLAB -go=k > GRID_TSExp.out";

# Implicit time-stepping + Yale
$cmd[1] = "setsid $cgwh -noplot eig.cmd -g=GRID -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=yale -orderOfExtrapolation=-1 -freq=12 -numEigs=16 -matlab=MATLAB -go=k > GRID_TSImpYale.out";

# Implicit time-stepping + MG solver
$cmd[2] = "setsid $cgwh -noplot eig.cmd -g=GRID -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16 -matlab=MATLAB -go=k > GRID_TSImpMG.out";

# Implicit time-stepping + bi-CG-Stab solver
$cmd[3] = "setsid $cgwh -noplot eig.cmd -g=GRID -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=bcgs -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16 -matlab=MATLAB -go=k  > GRID_TSImpBiCGStab.out";

$numSolvers=4; 

@impSolver = ("ExpTS","ImpTSDirect","ImpTSMG","ImpTSBiCGStab");
$g=1; 
if( $grid eq "square" || $grid eq "nonSquare" )
{
  $grids[$g]= $grid ."128.order"; $g++;
  $grids[$g]= $grid ."256.order"; $g++;
  $grids[$g]= $grid ."512.order"; $g++;
  $grids[$g]= $grid ."1024.order"; $g++;
  $grids[$g]= $grid ."2048.order"; $g++;
}
elsif( $grid eq "box" || $grid eq "nonBox" )
{
  $grids[$g]= $grid ."1.order"; $g++;
  $grids[$g]= $grid ."2.order"; $g++;
  $grids[$g]= $grid ."4.order"; $g++; 
  $grids[$g]= $grid ."8.order"; $g++;  
  $grids[$g]= $grid ."16.order"; $g++;      
}
else
{
  printf("UNKNOWN grid=[$grid]\n");
  exit(1);

}


$orderStart=$order; $orderEnd=$order; 
for( $order=$orderStart; $order<=$orderEnd; $order=$order+2 )
{
  printf("\n =================== ORDER = $order ==========================\n");
  for( $g=$resStart; $g<=$numRes; $g++ )
  {
    $myGrid = $grids[$g] . $order; 
    $myGridName = $myGrid;
    $myGridName =~ s/\.//g; # remove "."
    printf("  ----------- GRID = $myGrid (g=$g, resStart=$resStart, numRes=$numRes) --------------\n");
    for( $i=0; $i<$numSolvers; $i++ )  
    {
       $myCmd = $cmd[$i];
       $myCmd =~ s/MATLAB/$myGridName$impSolver[$i]/g;
       $myCmd =~ s/GRID/$myGrid/g;

       printf("Run [$myCmd]\n");
       # Execute the command: 
       # my $output = `$myCmd`; 

    }
  }
}
# foreach $caseName ( @caseNames )  # process all cases
# {

#   printf("caseName=$caseName\n");

#   if( $caseName eq "square64" )
#   {
#     @deflateValues = (1,2,4,8);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=square64.order2 -x0=.55 -y0=.6 -nf=1 -freq=11 -amp=100 -beta=11 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=1000 -numPeriods=1 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig/square64O2Ev128.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=square64DeflateDEFLATE -krylovType=KRYLOVTYPE -gmresRestartLength=200 -go=ka";


#   }
#   elsif( $caseName eq "square128" )
#   {
#     @deflateValues = (0,1,16,32,64,128);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=square128.order4 -x0=.65 -y0=.7 -nf=1 -freq=40 -amp=1600 -beta=40 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig2/square128O4Ev512 -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=square128O4freq40DeflateDEFLATE  -krylovType=KRYLOVTYPE -gmresRestartLength=200 -go=ka";

#   }
#   elsif( $caseName eq "square128c" )
#   {
#     # Use coarse grid EVs
#     # @deflateValues = (0,1,16,32,64,128);
#     @deflateValues = (0,1,4,8,16,32,64);
#     # @deflateValues = (64);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=square128.order4 -x0=.65 -y0=.7 -nf=1 -freq=40 -amp=1600 -beta=40 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig/square64O4Ev512 -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=square128O4freq40EV64DeflateDEFLATE  -krylovType=KRYLOVTYPE -gmresRestartLength=200 -go=a";

#   }  
#   elsif( $caseName eq 'square256')
#   {
#     @deflateValues = (1,2,4,8,16,32,64,128,200);
#     # @deflateValues = (2,200);
#     # @deflateValues = (1);
#     # @deflateValues = (0,4,8,16,32,64,128,200);
#     # @deflateValues = (200);
#     $cmd = "setsid $cgwh -noplot gaussian.cmd -g=square256.order4 -x0=.62 -y0=.73 -nf=1 -freq=60 -amp=3600 -beta=60 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig2/square256O4Ev512 -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=square256O4freq60DeflateDEFLATE -krylovType=KRYLOVTYPE -gmresRestartLength=400 -go=ka";
#   }
#   elsif( $caseName eq 'square256c')
#   {
#     @deflateValues = (1,2,4,8,16,32,64);
#     @deflateValues = (100);
#     # @deflateValues = (128);
#     $cmd = "setsid $cgwh -noplot gaussian.cmd -g=square256.order4 -x0=.62 -y0=.73 -nf=1 -freq=60 -amp=3600 -beta=60 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig2/square128O4Ev512 -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=square256O4freq60EV128DeflateDEFLATE -krylovType=KRYLOVTYPE -gmresRestartLength=400 -go=a";
#   }  
#   elsif( $caseName eq "square512" )
#   {
#     @deflateValues = (1,8,16,32,64,128,200);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=square512.order4 -x0=.62 -y0=.73 -nf=1 -freq=100 -amp=10000 -beta=100 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig2/square512O4Ev1024 -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=square512O4freq100DeflateDEFLATE -krylovType=KRYLOVTYPE -gmresRestartLength=400 -deflateForcing=0 -go=ka";
#   }
#   elsif( $caseName eq "disk8" )
#   {
#     @deflateValues = (0,2,8,16,32,64,128,200);
#     # @deflateValues = (128,200);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=sice8.order4 -x0=-0.25 -y0=-0.2 -nf=1 -freq=30 -amp=900 -beta=30 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=1000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig/sice8O4EigsEv512.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=diskG8O4Freq30DeflateDEFLATE -krylovType=KRYLOVTYPE -gmresRestartLength=400 -go=ka";  
#   }
#   elsif( $caseName eq "disk8c" )
#   {
#     @deflateValues = (0,2,8,16,32,64);
#     # @deflateValues = (0);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=sice8.order4 -x0=-0.25 -y0=-0.2 -nf=1 -freq=30 -amp=900 -beta=30 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=1000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig/sice4O4EigsEv256.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=diskG8O4Freq30EVG4DeflateDEFLATE -krylovType=KRYLOVTYPE -gmresRestartLength=400 -go=a";  
#   }  
#   elsif( $caseName eq "disk16" )
#   {
#     @deflateValues = (0,2,8,16,32,64,128,170,200,250,300);
#     @deflateValues = (350);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=sice16.order4 -x0=-0.25 -y0=-0.2 -nf=1 -freq=50 -amp=2500 -beta=50 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=1000 -numPeriods=2 -cfl=10000 -eigenVectorFile=/home/henshw/runs/eig/sice16O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -matlab=diskG16O4Freq50DeflateDEFLATE -krylovType=KRYLOVTYPE -gmresRestartLength=400 -go=ka";  
#   }  
#   elsif( $caseName eq "knifeEdgeG32Np4" )
#   {
#     @deflateValues = (0,1,8,16,32,64,128,150);
#     # @deflateValues = (150);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=tipGridSmalle32.order4.ng3.hdf -ts=implicit -cfl=10000 -x0=-0.2 -y0=0.2 -nf=1 -freq=80 -amp=6400 -beta=80 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=4 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/tipGridSmallG32Order4EigsEv512.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -rtoli=1e-10 -atoli=1e-10 -matlab=knifeEdgeG32O4Freq80Np4DeflateDEFLATE -krylovType=KRYLOVTYPE -go=ka";   
#   }
#   elsif( $caseName eq "knifeEdgeG32Np2" )
#   {
#     @deflateValues = (0,1,8,16,32,64,128,150);
#     # @deflateValues = (150);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=tipGridSmalle32.order4.ng3.hdf -ts=implicit -cfl=10000 -x0=-0.2 -y0=0.2 -nf=1 -freq=80 -amp=6400 -beta=80 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=2 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/tipGridSmallG32Order4EigsEv512.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -rtoli=1e-10 -atoli=1e-10 -matlab=knifeEdgeG32O4Freq80Np2DeflateDEFLATE -krylovType=KRYLOVTYPE -go=ka";   
#   }  
#   elsif( $caseName eq "knifeEdgeG32Np1" )
#   {
#     @deflateValues = (0,1,8,16,32,64,128,150);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=tipGridSmalle32.order4.ng3.hdf -ts=implicit -cfl=10000 -x0=-0.2 -y0=0.2 -nf=1 -freq=80 -amp=6400 -beta=80 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=1 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/tipGridSmallG32Order4EigsEv512.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -rtoli=1e-10 -atoli=1e-10 -matlab=knifeEdgeG32O4Freq80Np1DeflateDEFLATE -krylovType=KRYLOVTYPE -go=ka";   
#   } 
#   elsif( $caseName eq "knifeEdgeG32ExplicitNp1" )
#   {
#     # explicit time-stepping
#     @deflateValues = (0,1,8,16,32,64,128,150);
#     @deflateValues = (100);

#     $cmd = "$cgwh -noplot gaussian.cmd -g=tipGridSmalle32.order4.ng3.hdf -ts=explicit -cfl=.95 -x0=-0.2 -y0=0.2 -nf=1 -freq=80 -amp=6400 -beta=80 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=1 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/tipGridSmallG32Order4EigsEv512.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -rtoli=1e-10 -atoli=1e-10 -matlab=knifeEdgeG32O4Freq80ExplicitNp1DeflateDEFLATE -krylovType=KRYLOVTYPE -go=ka";   
#   } 
#   elsif( $caseName eq "knifeEdgeG32Np32" )
#   {
#     @deflateValues = (0,1,8,16,32,64,128,150);
#     $cmd = "$cgwh -noplot gaussian.cmd -g=tipGridSmalle32.order4.ng3.hdf -ts=implicit -cfl=10000 -x0=-0.2 -y0=0.2 -nf=1 -freq=80 -amp=6400 -beta=80 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=10000 -numPeriods=32 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/tipGridSmallG32Order4EigsEv512.show -deflateWaveHoltz=1 -numToDeflate=DEFLATE -rtoli=1e-10 -atoli=1e-10 -matlab=knifeEdgeG32O4Freq80Np32DeflateDEFLATE -krylovType=KRYLOVTYPE -go=ka";   
#   }  
#   elsif( $caseName eq "doubleEllipseG4")    
#   {
#     # setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride4.order4.hdf -ts=implicit -cfl=1000 -x0=-1.4 -y0=-5.6 -nf=1 -amp=-400 -beta=10 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG4O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G4Order4Deflate64 -show=darkCornerRoomFreq10Np8G4Order4Deflate64.show -go=dfks > darkCornerRoomFreq10Np8G4Order4Deflate64.out
#   }
#   else
#   {
#     printf("ERROR: Unknown caseName=[$caseName]\n");
#   }


#   foreach $krylovType ( @krylovTypes )
#   {
#     foreach $deflate ( @deflateValues )
#     {
#       printf("--- krylovType=$krylovType, deflate=$deflate ---\n");
#       $myCmd = $cmd;
#       $myCmd =~ s/DEFLATE/$deflate/g;
#       printf("krylovType=$krylovType\n");
#       $myCmd =~ s/KRYLOVTYPE/$krylovType/;
#       printf("Run [$myCmd]\n");
#       # $rt = system($myCmd);
#       # printf("Done: rt=$rt\n");
#       my $output = `$myCmd`;

#       # printf("+++++++ Done: output: +++++++++++\n$output\n");       
#     }
#   }
#   # $rt = system($cmd);
#   # printf("Done: rt=$rt\n");
#   # my $output = `$cmd`;
#   # printf("+++++++ Done: output: +++++++++++\n$output\n");  

# }
  # # if( $fileName =~ /-cycles/ )
  # # {
  # #   $perfData=1; 
  # #   next;
  # # }
  # $idebug=0; 

  # open(FILE,"$fileName") || die "cannot open file $fileName!" ;
  # # open(OUTFILE,">junk.X") || die "cannot open output file junk.X!" ;
  
  # $it=1; # base 1 for matlab
  # while( <FILE> )
  # {
  #   $line = $_;
  #   printf("%s",$line);

  #   @tokens = split(' ',$line);
  #   # foreach my $token (@tokens) 
  #   $numTokens= $#tokens; 
  #   if( $idebug>2 )
  #   {
  #     for( $i=0; $i<=$numTokens; $i++ )
  #     {
  #       printf("token[$i]=$tokens[$i]\n");
  #     }
  #   }

  #   $nconv = $tokens[2]; # nconv=?
  #   $nconv =~ s/nconv=//; 
  #   $numConverged[$it]= $nconv;
  #   printf("numConverged[$it] = $numConverged[$it]\n");

  #   $num=1; 
  #   for( $i=5; $i<=$numTokens; $i=$i+2 )  
  #   {
  #     $eigenValueByIteration[$num][$it] = $tokens[$i];
  #     $err = $tokens[$i+1];
  #     $err =~ s/\((.*)\)/\1/; # remove surrounding ( )
  #     $eigenValueErrorByIteration[$num][$it] = $err;
  #     # printf("eigenValueByIteration[$num][$it] = $eigenValueByIteration[$num][$it], eigenValueErrorByIteration[$num][$it]=$eigenValueErrorByIteration[$num][$it]\n");
  #     $num=$num+1; 
  #   }  
  #   $numValues[$it]=$num; 

  #   $it = $it+1;

  # }
  # $numIterations=$it-1; 
  # close(FILE);

  # # Save results to a matlab file
  # $output=$fileName;
  # $output =~ s/.txt//; 
  # $output = $output . ".m";
  # open(OUTFILE,">$output") || die "cannot open file $output!" ;


  # printf(OUTFILE "numIterations=$numIterations;\n");
  # printf(OUTFILE "numConverged=[");
  # for( $it=1; $it<=$numIterations; $it++ )
  # {
  #   printf(OUTFILE " $numConverged[$it] ");
  #   if( $it<$#numConverged ){ printf(OUTFILE ","); }
  # }
  # printf(OUTFILE "];\n");

  # for( $it=1; $it<=$numIterations; $it++ )
  # {
  #   printf(OUTFILE "\n%% Iteration $it :\n");
  #   for( $num=1; $num<$numValues[$it]; $num++ )
  #   {
  #     printf(OUTFILE "eigenValueByIteration($num,$it)=$eigenValueByIteration[$num][$it];\n");
  #     printf(OUTFILE "eigenValueErrorByIteration($num,$it)=$eigenValueErrorByIteration[$num][$it];\n");
  #   }
  # }

  # close(OUTFILE);
  # printf("Wrote file =[$output]\n");


