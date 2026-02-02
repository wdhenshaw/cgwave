#include "CgWave.h"
#include "Ogshow.h"
#include "HDF_DataBase.h"

// =================================================================================================
/// \brief save the current solution (and errors) to the show file
// =================================================================================================
int CgWave::saveShow( int current, Real t, Real dt )
{
  const bool & saveShowFile = dbase.get<bool>("saveShowFile"); 

  if( !saveShowFile )
    return 0;

  Ogshow *& pShowFile = dbase.get<Ogshow*>("showFile"); 
  assert( pShowFile !=NULL );

  Ogshow & showFile = *pShowFile;


  // Build the augmented solution (may include errors)
  realCompositeGridFunction v;
  realCompositeGridFunction & u = getAugmentedSolution(current,v,t);  // u is either solution or v


  showFile.startFrame();

  // -- save t and dt --
  HDF_DataBase *dbp=NULL;

  dbp = showFile.getFrame();
  assert( dbp!=NULL );
  // save parameters that go in this frame
  HDF_DataBase & db = *dbp;
  db.put(t,"time");
  db.put(dt,"dt");

  // -- save plot titles --
  char buffer[80]; 
  aString showFileTitle[5];

  aString methodName = getMethodName();  // FD22, or FD24s etc.
  char buff[80];
  showFileTitle[0]=sPrintF(buffer,"CgWave %s",(const char *)methodName);
  showFileTitle[1]=sPrintF(buffer,"t=%4.3f, dt=%8.2e",t,dt);
  showFileTitle[2]="";  // marks end of titles

  for( int i=0; showFileTitle[i]!=""; i++ )
    showFile.saveComment(i,showFileTitle[i]);

  showFile.saveSolution( u,"u" );  // save the grid function


  // *wdh* this may now work in parallel *check me* (see Maxwell solver cgmx)
#ifndef USE_PPP
  if( pShowFile !=NULL &&
      pShowFile->isLastFrameInSubFile() )
  {
    saveSequencesToShowFile();
  }
#endif  

  showFile.endFrame();

  return 0;
}

//=========================================================================================
/// \brief  Save info into the time history arrays.
// 
//=========================================================================================
int CgWave::
saveSequenceInfo( real t0, RealArray & sequenceData )
{
  Ogshow *& show = dbase.get<Ogshow*>("showFile");
  if( show==NULL )
    return 0;

  int & sequenceCount      = dbase.get<int>("sequenceCount");
  int & numberOfSequences  = dbase.get<int>("numberOfSequences");
  RealArray & timeSequence = dbase.get<RealArray>("timeSequence");
  RealArray & sequence     = dbase.get<RealArray>("sequence");

  if( sequenceCount >= timeSequence.getLength(0) )
  {
    int num=timeSequence.getLength(0);
    Range R(0,num-1),all;
    RealArray seq;  seq=sequence;
    num=int(num*1.5+100);
    timeSequence.resize(num);
    sequence.redim(num,numberOfSequences);
    sequence(R,all)=seq;
  }

  timeSequence(sequenceCount)=t0;
  for( int n=sequenceData.getBase(0); n<=sequenceData.getBound(0); n++ )
    sequence(sequenceCount,n)=sequenceData(n); 

  sequenceCount++;
  return 0;
}



//=========================================================================================
//=========================================================================================
int CgWave::
saveSequencesToShowFile()
//=========================================================================================
{
  Ogshow *& show = dbase.get<Ogshow*>("showFile");
  const int & sequenceCount = dbase.get<int>("sequenceCount");
  if( show==NULL || sequenceCount<=0 )
    return 0;
  
  int & numberOfSequences  = dbase.get<int>("numberOfSequences");
  RealArray & timeSequence = dbase.get<RealArray>("timeSequence");
  RealArray & sequence     = dbase.get<RealArray>("sequence");
  int & saveMaxErrors      = dbase.get<int>("saveMaxErrors"); 

  Range I(0,sequenceCount-1);
  Range N=sequence.dimension(1);
  
  aString *name = new aString [numberOfSequences]; 
  aString sequenceName;
  sequenceName="solutionNorms";

  // if( checkErrors )
  //   sequenceName="errors";
  // else
  //   sequenceName="solutionNorms";

  const int numberOfComponents = 1;
  realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
  aString buff;
  for( int n=0; n<numberOfComponents; n++ )
  {
    // if( checkErrors )
    // {
    //   if( cgerrp!=NULL )
    //   {
    //     name[n]=cgerrp[0].getName(n);
    //   }
    //   else
    //   {
    //     name[n]=sPrintF("error%i",n);
    //   }
    // }
    // else
    // {
    //   // -- save solution norms --
    //   aString buff;
    //   name[n]= sPrintF(buff,"%sNorm",(const char*)u[0].u.getName(n));
    // }

    name[n]= sPrintF(buff,"%sNorm",(const char*)u[0].getName(n));
    for( int i=0; i<name[n].length(); i++ )
    { // change blanks to underscores (for matlab)
      if( name[n][i]==' ' )
      {
        name[n][i]='_';
      }
    }
  }
  
  const int & computeEnergy  = dbase.get<int>("computeEnergy"); 
  if( computeEnergy )
  {
    name[computeEnergy] = "energy";
  }


  if( saveMaxErrors )
  {
    name[saveMaxErrors] = "maxErr";
  }

  // display(sequence(I,N),"saveSequencesToShowFile: sequence(I,N)");
  
  show->saveSequence(sequenceName,timeSequence(I),sequence(I,N),name);

  delete [] name;
  
  return 0;
}


#include "ShowFileReader.h"

// =================================================================================================
/// \brief read a solution from a check point file (used for WaveHoltz and EigenWave)
/// \param solution (input) : solution to read, 1,2,3,...
/// \Return value: 
///      1=requested solution was read in, 
///      0=requested solution is not in the file, or we are not reading the check point file.
// =================================================================================================
int CgWave::readCheckPoint( int solution, RealCompositeGridFunction & v )
{
  const bool & readCheckPointFile = dbase.get<bool>("readCheckPointFile"); 

  if( !readCheckPointFile )
    return 0;


  if( !dbase.has_key("checkPointShowFileReader") )
  {
    ShowFileReader *& showFilePointer = dbase.put<ShowFileReader*>("checkPointShowFileReader");

    showFilePointer = new ShowFileReader;  // DELETE ME 
    assert( showFilePointer !=NULL );

    ShowFileReader & showFileReader = *showFilePointer;  

    aString & nameOfCheckPointFile = dbase.get<aString>("readCheckPointFileName"); // name of the check point file

    printF("Open the check point file=[%s]\n",(const char*)nameOfCheckPointFile);

    showFileReader.open(nameOfCheckPointFile);

    int & numCheckPointSolutions = dbase.put<int>("numCheckPointSolutions");

    int numFrames          = showFileReader.getNumberOfFrames();
    numCheckPointSolutions = showFileReader.getNumberOfSolutions();

    printF("Open the old check point file=[%s], There are %d solutions and %d frames\n",
         (const char*)nameOfCheckPointFile, numCheckPointSolutions,numFrames);
  }

  const int & numCheckPointSolutions = dbase.get<int>("numCheckPointSolutions");
  if( solution>numCheckPointSolutions )
    return 0;

  ShowFileReader & showFileReader = *dbase.get<ShowFileReader*>("checkPointShowFileReader");

  printF("readCheckPointFile: reading solution=%d\n",solution);

  CompositeGrid cgs;
  showFileReader.getASolution(solution,cgs,v);

  v.updateToMatchGrid(cg);


  return 1;


  // FINISH ME -- close file when done
  // showFileReader.close();


  return 0;
}


// =================================================================================================
/// \brief save a solution to the check point file (used for WaveHoltz and EigenWave)
/// \param solution (input): not currently used
// =================================================================================================
int CgWave::saveCheckPoint( int solution, RealCompositeGridFunction & v )
{
  const bool & saveCheckPointFile = dbase.get<bool>("saveCheckPointFile"); 

  if( !saveCheckPointFile )
    return 0;

  printF("saveCheckPointFile: solution=%d\n",solution);

  Ogshow *& pCheckPointFile = dbase.get<Ogshow*>("checkPointFile"); 
  assert( pCheckPointFile !=NULL );

  Ogshow & checkPointFile = *pCheckPointFile;

  checkPointFile.startFrame();

  // -- save plot titles --
  char buffer[80]; 
  aString showFileTitle[5];

  aString methodName = getMethodName();  // FD22, or FD24s etc.
  char buff[80];
  showFileTitle[0]=sPrintF(buffer,"CgWave %s",(const char *)methodName);
  showFileTitle[1]="";  // marks end of titles

  for( int i=0; showFileTitle[i]!=""; i++ )
    checkPointFile.saveComment(i,showFileTitle[i]);

  v.setName("u");  
  v.setName("u",0);  
  checkPointFile.saveSolution( v,"u" );  // save the grid function


  checkPointFile.endFrame();

  return 0;
}
