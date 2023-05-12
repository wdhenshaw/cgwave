eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to extract the performance results from SLEPc
#  usage: 
#         extractConvergence.p file
# 

@fileNames = @ARGV;

# $perfData=0; # extract CPU data 
foreach $fileName ( @fileNames )  # process all files
{
  # if( $fileName =~ /-cycles/ )
  # {
  #   $perfData=1; 
  #   next;
  # }
  $idebug=0; 

  open(FILE,"$fileName") || die "cannot open file $fileName!" ;
  # open(OUTFILE,">junk.X") || die "cannot open output file junk.X!" ;
  
  $it=1; # base 1 for matlab
  while( <FILE> )
  {
    $line = $_;
    printf("%s",$line);

    @tokens = split(' ',$line);
    # foreach my $token (@tokens) 
    $numTokens= $#tokens; 
    if( $idebug>2 )
    {
      for( $i=0; $i<=$numTokens; $i++ )
      {
        printf("token[$i]=$tokens[$i]\n");
      }
    }

    $nconv = $tokens[2]; # nconv=?
    $nconv =~ s/nconv=//; 
    $numConverged[$it]= $nconv;
    printf("numConverged[$it] = $numConverged[$it]\n");

    $num=1; 
    for( $i=5; $i<=$numTokens; $i=$i+2 )  
    {
      $eigenValueByIteration[$num][$it] = $tokens[$i];
      $err = $tokens[$i+1];
      $err =~ s/\((.*)\)/\1/; # remove surrounding ( )
      $eigenValueErrorByIteration[$num][$it] = $err;
      # printf("eigenValueByIteration[$num][$it] = $eigenValueByIteration[$num][$it], eigenValueErrorByIteration[$num][$it]=$eigenValueErrorByIteration[$num][$it]\n");
      $num=$num+1; 
    }  
    $numValues[$it]=$num; 

    $it = $it+1;

  }
  $numIterations=$it-1; 
  close(FILE);

  # Save results to a matlab file
  $output=$fileName;
  $output =~ s/.txt//; 
  $output = $output . ".m";
  open(OUTFILE,">$output") || die "cannot open file $output!" ;


  printf(OUTFILE "numIterations=$numIterations;\n");
  printf(OUTFILE "numConverged=[");
  for( $it=1; $it<=$numIterations; $it++ )
  {
    printf(OUTFILE " $numConverged[$it] ");
    if( $it<$#numConverged ){ printf(OUTFILE ","); }
  }
  printf(OUTFILE "];\n");

  for( $it=1; $it<=$numIterations; $it++ )
  {
    printf(OUTFILE "\n%% Iteration $it :\n");
    for( $num=1; $num<$numValues[$it]; $num++ )
    {
      printf(OUTFILE "eigenValueByIteration($num,$it)=$eigenValueByIteration[$num][$it];\n");
      printf(OUTFILE "eigenValueErrorByIteration($num,$it)=$eigenValueErrorByIteration[$num][$it];\n");
    }
  }

  close(OUTFILE);
  printf("Wrote file =[$output]\n");


}