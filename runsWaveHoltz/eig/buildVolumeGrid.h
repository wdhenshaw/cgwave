# 
# -- Make hyperbolic grid --
#
if( $gridName eq "" ){ $gridName="myGrid"; }
if( $numberOfVolumeSmooths eq "" ){ $numberOfVolumeSmooths=20; }  #
if( $nr0 eq "" ){ $nr0=8; }  #
if( $direction eq "" ){ $direction="backward"; }
if( $bc eq "" ){ $bc="1 2 3 4"; }
if( $share eq "" ){ $share="0 0 0 0"; }
# 
  $nr = intmg( $nr0 + ($order-2)/2 );
  hyperbolic
    Start curve:$startCurve
    $direction
    $nDist=($nr-3.25)*$ds;
    distance to march $nDist
    $nrm= intmg($nr)-1; 
    lines to march $nrm
    $nTheta = intmg($arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    spacing: geometric
    geometric stretch factor 1.05 
    #
    BC: left fix x, float y and z
    BC: right fix x, float y and z
    #
    generate
    #
    # Set order of interpolation for data point mapping: 
    if( $order eq 2 ){ $dmpOrder = "second order"; }else{ $dmpOrder = "fourth order" }
    printf("hypeOrder = $dmpOrder\n");
    $dmpOrder
    #     
    boundary conditions
      $bc 
      # 1 2 3 0 0 0 
    share 
      $share
      # 1 2 0 0 0 0
    name $gridName
    # open graphics
  exit
