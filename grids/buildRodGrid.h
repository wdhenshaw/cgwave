printf("buildRodGrid: count=$count xShiftr=$xShiftr[$count] yShiftr=$yShiftr[$count] angler=$angler[$count]\n");
# pause
$count=$count+1; $gridName="rod$count"; $gridNames .= "\n" . $gridName;
convertToNurbs($mapName,$gridName,$angler[$count-1],$rotationAxis,$xShiftr[$count-1],$yShiftr[$count-1],$zShift);
$cmds
# pause