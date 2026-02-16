#
# Configuration file for multiRodGrid.cmd 
#
#     
# Box bounds: 
$xa=-1.5; $xb=2.5; $ya=-2; $yb=2;
#
# ------ list of rod parameters: --------
#
$deltar = $rodHeight/sqrt(2); 
$gapy= $deltar*.65; 
$gapx= 1.25*$deltar; 
$i=0;
$xShiftr[$i]=0*$gapx; $yShiftr[$i]=-3*$gapy; $angler[$i]=-45; $i=$i+1; 
$xShiftr[$i]=0*$gapx; $yShiftr[$i]=-2*$gapy; $angler[$i]=-45; $i=$i+1; 
$xShiftr[$i]=0*$gapx; $yShiftr[$i]=-1*$gapy; $angler[$i]=-45; $i=$i+1; 
$xShiftr[$i]=0*$gapx; $yShiftr[$i]= 0*$gapy; $angler[$i]=-45; $i=$i+1; 
$xShiftr[$i]=0*$gapx; $yShiftr[$i]=+1*$gapy; $angler[$i]=-45; $i=$i+1;
$xShiftr[$i]=0*$gapx; $yShiftr[$i]=+2*$gapy; $angler[$i]=-45; $i=$i+1;
$xShiftr[$i]=0*$gapx; $yShiftr[$i]=+3*$gapy; $angler[$i]=-45; $i=$i+1;
# 
$xShiftr[$i]=1*$gapx; $yShiftr[$i]=-3*$gapy;  $angler[$i]=+45; $i=$i+1; 
$xShiftr[$i]=1*$gapx; $yShiftr[$i]=-2*$gapy;  $angler[$i]=+45; $i=$i+1; 
$xShiftr[$i]=1*$gapx; $yShiftr[$i]=-1*$gapy;  $angler[$i]=+45; $i=$i+1; 
$xShiftr[$i]=1*$gapx; $yShiftr[$i]= 0*$gapy;  $angler[$i]=+45; $i=$i+1; 
$xShiftr[$i]=1*$gapx; $yShiftr[$i]= 1*$gapy;  $angler[$i]=+45; $i=$i+1; 
$xShiftr[$i]=1*$gapx; $yShiftr[$i]= 2*$gapy;  $angler[$i]=+45; $i=$i+1; 
$xShiftr[$i]=1*$gapx; $yShiftr[$i]= 3*$gapy;  $angler[$i]=+45; $i=$i+1; 
#
$numRods=$i;