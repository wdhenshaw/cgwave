#
#  Command file for genEigs to orthogonalize results from eigenWave
# Usage: 
#     genEigs [-noplot] orthogonalize.cmd -g=<gridName> -evFile=<s> -show=<s> -orthogonalize=1 -go=go 
#
#
$show=""; $go="og"; $matlab=1; 
$multiEigTol=1e-5; # tolerance for checking if two eigenvalues are the same
GetOptions( "show=s"=>\$show,"matlab=i"=>\$matlab,"multiEigTol=f"=>\$multiEigTol,"go=s"=>\$go );
#
multiple eigenvalue tol: $multiEigTol
read eigenvectors from a file
#
#open graphics
#
orthogonalize
# compute
# open graphics
if( $show ne "" ){ $cmd="save show file"; }else{ $cmd="#"; }
$cmd
#
# if( $matlab eq "1" ){ $cmd="save matlab file"; }else{ $cmd="#"; }
# $cmd
$cmd="#"; 
if( $go eq "og" ){ $cmd = "open graphics"; }
if( $go eq "go" ){ $cmd = "exit"; }
$cmd
