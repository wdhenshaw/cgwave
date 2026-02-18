#
# Configuration file for multiRodGrid.cmd 
#
#     
# Box bounds: 
$xa=-1.5; $xb=1.5; $ya=-2; $yb=2;
#
# ------ list of rod parameters: --------
#
$deltar = $rodHeight/sqrt(2); 
if( $rodWidth ge 0.1 ){ $gapy= $deltar*.55; $num=8; }else{ $gapy= $deltar*.45; $num=10; }
$gapx= 1.25*$deltar; 
$i=0; 
$thetaStart=-45; $thetaEnd=$thetaStart-90; 
$yr=-($num/2-.5)*$gapy; 
$dTheta=($thetaEnd-$thetaStart)/($num-1);
for( $j=0; $j<$num; $j++ ){\
$xShiftr[$j]=0*$gapx; $yShiftr[$j]=$yr; $angler[$j]=$thetaStart+$j*$dTheta; $yr=$yr+$gapy; $i=$i+1;\
} 
# 
#
$numRods=$i;