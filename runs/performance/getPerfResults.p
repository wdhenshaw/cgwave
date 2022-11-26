eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
# perl program to exttract the performance results from CgWave output
#  usage: 
#         getPerfResults [-cycles] file
# 

@fileNames = @ARGV;

$perfData=0; # extract CPU data 
foreach $fileName ( @fileNames )  # process all files
{
  if( $fileName =~ /-cycles/ )
  {
    $perfData=1; 
    next;
  }

  open(FILE,"$fileName") || die "cannot open file $fileName!" ;
  # open(OUTFILE,">junk.X") || die "cannot open output file junk.X!" ;
  
  while( <FILE> )
  {
    $line = $_;
    if( $perfData==0 )
    {
      if( $line =~ /\% PerfInfo/ )
      {
        # printf("%s",$line);
        # Clean up the line : 
        $line =~ s/ \% PerfInfo//;
        $line =~ s/\.order..*\.hdf//;  # remove .order*.hdf 

        printf("%s",$line);

      }
    }
    else
    {
      if( $line =~ /\% CyclesPerfInfo/ )
      {
        # printf("%s",$line);
        # Clean up the line : 
        $line =~ s/ \% CyclesPerfInfo//;
        $line =~ s/\.order..*\.hdf//;  # remove .order*.hdf 

        printf("%s",$line);  
      }    
    }

  }

  # close(OUTFILE);
  close(FILE);


}