eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to extract performance data and put into Matlab arrays
#  usage: 
#         getPerfArrays.p perf.dat
# 

@fileNames = @ARGV;

foreach $fileName ( @fileNames )  # process all files
{
  
  open(FILE,"$fileName") || die "cannot open file $fileName!" ;
  # open(OUTFILE,">junk.X") || die "cannot open output file junk.X!" ;
  
  @gridNames = ();
  $gridNum=-1; 
  @solverNames = ();
  $cpuTotal  = 4; # column location of the total cpu
  $cpuSolve  = 5; # time for time solve, includes arc, bc, upw, interp
  $cpuARC    = 6; # time for advance rectangular/curvilinear
  $cpuUpw    = 7; # time for upwind
  $cpuBC     = 8; # bc
  $cpuInterp = 9; # interp
  printf("cpuTotal  =1; %% location in array with cpu total\n");
  printf("cpuSolve  =2; %% location in array with cpu for solve\n");
  printf("cpuARC    =3; %% location in array with cpu adv rectangular/curvilinear\n");
  printf("cpuUpw    =4; %% location in array with cpu upwinding/dissipation\n");
  printf("cpuBC     =5; %% location in array with cpu boundary conditions\n");
  printf("cpuInterp =6; %% location in array with cpu interpolation\n");
  while( <FILE> )
  {
    $line = $_;
    chop($line);
    if( $line =~ /&/ )
    {
      @cols = split('&', $line);

      # printf("%s",$line);
      if( 1==0 )
      {
        for( $i=0; $i< @cols; $i++ )
        {
          printf("cols[$i]=$cols[$i]\n");
        }
      }
      # Get grid name: 
      $grid=$cols[0];
      $grid =~ s/^[ ]*//;  # remove leading blanks
      $grid =~ s/[ ]*$//;  # ... and trailing    
      if( $gridNum<0 || $grid ne $gridNames[$gridNum] )
      {
        $gridNum=$gridNum+1;
        $gridNames[$gridNum]=$grid;
        printf("$grid = %d; %% grid name enumerator\n",$gridNum+1);
      }
      $order = $cols[1];

      # get scheme name 
      $solver = $cols[2];
      $solver =~ s/^[ ]*//;  # remove leading blanks
      $solver =~ s/[ ]*$//;
      $solverNum=-1;
      for( $i=0; $i<@solverNames; $i++ )
      {
        if( $solver eq $solverNames[$i] )
        {
          $solverNum=$i; break;
        }
      }
      if( $solverNum==-1 )
      { $solverNames[@solverNames]=$solver; $solverNum=@solverNames-1; 
        printf("solverName{%d}=\"$solver\";\n",$solverNum+1); 
      }

      $cols[$cpuInterp] =~ s/[ ]*\\\\//; # remove ending "//" from last entry

      # $storageValue = $cols[$storage];
      # $storageValue =~ s/\\\\//; # remove ending //
      $sn=$solverNum+1;
      printf("data(cpuTotal, $order,$sn,$gridNames[$gridNum])=$cols[$cpuTotal]; %% data(:,order,solver,grid)= value\n");
      printf("data(cpuSolve, $order,$sn,$gridNames[$gridNum])=$cols[$cpuSolve]; \n");
      printf("data(cpuARC,   $order,$sn,$gridNames[$gridNum])=$cols[$cpuARC]; \n");
      printf("data(cpuUpw,   $order,$sn,$gridNames[$gridNum])=$cols[$cpuUpw]; \n");
      printf("data(cpuBC,    $order,$sn,$gridNames[$gridNum])=$cols[$cpuBC]; \n");
      printf("data(cpuInterp,$order,$sn,$gridNames[$gridNum])=$cols[$cpuInterp]; \n");


    }

  }

  # close(OUTFILE);
  close(FILE);


}