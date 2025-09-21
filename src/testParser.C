
#include "Ogshow.h"  
#include "PlotStuff.h"
#include "display.h"

// int 
// getLineFromFile( FILE *file, char s[], int lim);



//=============================================================================================================
//    Test Overture parser for Jeff
// 
//=============================================================================================================
int 
main(int argc, char *argv[])
{
  int debug=0;

  Overture::start(argc,argv);  // initialize Overture
  const int myid=Communication_Manager::My_Process_Number;
  const int np = Communication_Manager::numberOfProcessors();
  
  printF("Usage: `testParser [-noplot] [-cmd=file.cmd] [-g=<gridName>]'\n");
  
  aString nameOfOGFile= "square32.order2.hdf";

  aString commandFileName="";
  aString line;
  int len=0;
  bool plotOption=true;  // by default we plot interactively

  int numParGhost=2; // number of parallel ghost : 2 for order2, 3 for order 4
  

  if( argc > 1 )
  {
    int i;
    for( i=1; i<argc; i++ )
    {
      line=argv[i];
      if( line=="-noplot" || line=="noplot" )
        plotOption=false; 
      else if( len=line.matches("-g=") )
      {
        nameOfOGFile=line(len,line.length()-1);
        // printf("\n$$$$ node %i : use grid=[%s]\n",myid,(const char*)nameOfOGFile);
      }
      else if( len=line.matches("-cmd=") )
      {
        commandFileName=line(len,line.length()-1);
        // printf("\n$$$$ node %i : read command file %s\n",myid,(const char*)commandFileName);
      }
      else if( commandFileName == "" )
      {
        commandFileName=line(len,line.length()-1);
        printF("\nread command file %s\n",(const char*)commandFileName);
      }
      else if( len=line.matches("-numParGhost=") )
      {
        sScanF(line(len,line.length()-1),"%i",&numParGhost);
        printF("Setting numParGhost=%i\n",numParGhost);
      }

      // old way
      else if( len=line.matches("-grid=") )
      {
        nameOfOGFile=line(len,line.length()-1);
        // printf("\n$$$$ node %i : use grid=[%s]\n",myid,(const char*)nameOfOGFile);
      }
 
    }
  }

  
  GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("cgWave",plotOption,argc,argv));

  PlotStuffParameters psp;

  // By default start saving the command file called "cgWave.cmd"
  aString logFile="parser.cmd";
  ps.saveCommandFile(logFile);
  printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

  if( commandFileName!="" )
    ps.readCommandFile(commandFileName);

  aString answer;
  for( int i=0; i<20; i++ )
  {
    ps.inputString(answer,"Enter a response");
    printF("answer=[%s]\n",(const char*)answer);
  }



  
  Overture::finish();          
  return 0;
}
