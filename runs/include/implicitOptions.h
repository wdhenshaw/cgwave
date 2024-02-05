#
#  Include file to set options for implicit time-stepping
#
# ---- set defaults:  *finish me*
#
if( $ts eq "" ){ $ts="implicit"; }
if( $implicitUpwind eq "" ){ $implicitUpwind=0; }
if( $solveri eq "" ){ $solveri="yale"; }
# weights in implicit time-stepping:
if( $beta2 eq "" ){ $beta2=0.5; }
if( $beta4 eq "" ){ $beta4=0; }
if( $beta6 eq "" ){ $beta6=0; }
if( $beta8 eq "" ){ $beta8=0; }
if( $rtoli eq "" ){ $rtoli=1.e-8; }
if( $atoli eq "" ){ $atoli=1.e-8; }
if( $maxiti eq "" ){ $maxiti=500; }
if( $mgMaxIts eq "" ){ $mgMaxIts=15; } # max MG iterations
if( $debugmg eq "" ){ $debugmg=0; }
if( $debugOges eq "" ){ $debugOges=0; }
if( $rectangular eq "" ){ $rectangular="implicit"; } # set  $rectangular="explicit" to treat rectangular grids explicitly 
#
if( $ts eq "implicit" ){ $cmd="choose grids for implicit\n  rectangular=$rectangular\n done"; }else{ $cmd="#"; }
$cmd
# implicitUpwind = 1 : include upwinding in implicit matrix
implicit upwind $implicitUpwind
#
implicit weights $beta2 $beta4 $beta6 $beta8
#
# $solveri="multigrid";
# $solveri="yale";
implicit solver parameters
  # NOTE: bcgs = bi-CG stab
  if( $solveri ne "yale" ){ $cmd="choose best iterative solver\n $solveri"; }else{ $cmd="choose best direct solver"; }
  $cmd
  number of incomplete LU levels
    3
  number of GMRES vectors
    20
  maximum number of iterations
    #
    $maxiti
  relative tolerance
    # printf("implicitOptions: rtoli=$rtoli, atoli=$atoli\n"); 
    $rtoli
  absolute tolerance
    $atoli
  debug
    $debugOges
 # 
  multigrid parameters
    choose good parameters: 1
    residual tolerance $rtoli
    absolute tolerance $atoli
    maximum number of iterations
      $mgMaxIts
    debug
      $debugmg
    if( $debugmg > 1 ){ $cmdmg="show smoothing rates"; }else{ $cmdmg="#"; }
    $cmdmg
    # printf("debugmg=$debugmg, cmdmg = $cmdmg\n");
    # pause
     # maximum number of extra levels
     #  1
    #
    # do not average coarse grid equations
    # pause 
    # Coarse level solver:  
    Oges parameters
      choose best direct solver
      # choose best iterative solver
      relative tolerance
        $rtoli
        # 1.e-10
      absolute tolerance
        $atoli
        # 1.e-10
      number of incomplete LU levels
       3
      # debug
      #   7
     minimum number of iterations
       2        
    exit    
  exit
  #pause
exit