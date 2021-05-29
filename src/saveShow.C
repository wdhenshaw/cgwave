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

  showFile.endFrame();

  return 0;
}