//
// Test the computation of the (parallel) global index for PETSc and SLEPc 
// 

#include "Overture.h"
#include "GL_GraphicsInterface.h"
#include "display.h"
#include "FileOutput.h"
#include "ColourBar.h"
#include "ParallelUtility.h"
#include "xColours.h"

#include "CgWave.h"

GL_GraphicsInterface *psPointer; // create a (pointer to a) GL_GraphicsInterface object

void selectObject(const real & x=-1., const real & y=-1.);
void getCursor( real & x, real & y );

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )


#define FOR_3(i1,i2,i3,I1,I2,I3) \
  for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  \
  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  \
  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  \

//===============================================================================================
//
// Test the computation of the (parallel) global index for PETSc and SLEPc 
// 
//===============================================================================================


int 
main(int argc, char *argv[])
{
  
  PlotIt::parallelPlottingOption=1;  // turn on distributed plotting! *******************

  const bool showComputedGeometry=false;

  Overture::start(argc,argv);  // initialize Overture

  const int myid=max(0,Communication_Manager::My_Process_Number);
  const int np = Communication_Manager::numberOfProcessors();  

  Overture::turnOnMemoryChecking(true);

  printF("Usage:testGlobalIndex [-noplot] -g=gridName -cmd=file.cmd\n");
  aString nameOfOGFile="cice2.order2.hdf";
  aString commandFileName="";
  int plotOption=false;

  int numParGhost=2; // number of parallel ghost : 2 for order2, 3 for order 4

  if( argc>1 )
  {
    aString line;
    int len=0;
    for( int i=1; i<argc; i++ )
    {
      //  int len=strlen(argv[i]);
      //        printF(" argv[%i]=%s len=%i\n",i,argv[i],len);
      
      line=argv[i];
      // if( line=="plot" || line=="-plot" )
      //   plotOption=true;
      // else if( line=="noplot" || line=="-noplot" )
      //   plotOption=false;
      if( len=line.matches("-g=") )
      {
        nameOfOGFile=line(len,line.length()-1);
      }
      else if( len=line.matches("-cmd=") )
      {
        commandFileName=line(len,line.length()-1);
      }
    }
  }
  // else
  //  ps.inputFileName(nameOfOGFile, ">> Enter the name of the (old) composite grid file:", ".hdf");

  // create and read in a CompositeGrid
  #ifdef USE_PPP
    // On Parallel machines always add at least this many ghost lines on local arrays
    // const int numGhost=2;
    MappedGrid::setMinimumNumberOfDistributedGhostLines(numParGhost); 
  #endif
  CompositeGrid cg;
  bool loadBalance=true; // turn on or off the load balancer
  getFromADataBase(cg,nameOfOGFile,loadBalance);

  cg.update(MappedGrid::THEmask);
  const int numberOfDimensions = cg.numberOfDimensions();
  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  bool useNewVersion=true; // use version in CgWave
  if( useNewVersion )
  {
    GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("testGlobalIndex",plotOption,argc,argv));  	
	  CgWave cgWave(cg,ps); 


	  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	  {
	  	MappedGrid & mg = cg[grid];
	  	for( int axis=0; axis<numberOfDimensions; axis++ )
	  		for( int side=0; side<=1; side++ )
	  	    mg.setBoundaryCondition( side,axis,CgWave::dirichlet );
	  }	   

	  cgWave.initializeGlobalIndexing();

	  // ---- NOW TEST ----

	  // totalActive = total number of active points
	  // numActiveLocal = number of active points on this processor 
	  const int & totalActive    = cgWave.dbase.get<int>("totalActive");
	  const int & numActiveLocal = cgWave.dbase.get<int>("numActiveLocal");	  
    IntegerArray *&indexVector = cgWave.dbase.get<IntegerArray*>("indexVector");	  

	  RealArray v(totalActive); // just make a serial arry with all points for testing 
	  for( int i=0; i<totalActive; i++ )
	  {
	  	v(i)=i;
	  }

	  realCompositeGridFunction u(cg);
	  u=0.;

    Index I1,I2,I3;

	  // ----------------------------------
	  // ----- vector to grid function ----
	  // ----------------------------------

	  cgWave.vectorToGridFunction( v.getDataPointer(), u );


	  if( 1==0 )
		  u.display("u after vector to grid function","%5.0f ");

	  // ----------------------------------
	  // ----- grid function to vector ----
	  // ----------------------------------

	  RealArray v2(totalActive); // just make a serial arry with all points for testing 
	  v2=0;

    cgWave.gridFunctionToVector( u, v2.getDataPointer() );	  


	  // -- check result
	  RealArray v2Buffer(totalActive);
	  ParallelUtility::getSums( v2.getDataPointer(),v2Buffer.getDataPointer(),totalActive );
	  v2 = v2Buffer;

		Real maxErr = max(fabs(v-v2)); 
		printf("Max error after GF to VECTOR and reverse = %9.2e\n",maxErr);



  }
  else
  {
  	// ******************* INITIAL VERSION *************

	  printf("testGlobalIndex: START np=%d, myid=%d, ...\n",np,myid );

	  // -- store local/global index of each grid point in these arrays:
	  IntegerArray *indexVector = new IntegerArray[numberOfComponentGrids]; // delete me when done

	  Index I1,I2,I3;
	  IntegerArray numActive(numberOfComponentGrids,np);
	  numActive=0; 

	  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	  {
	  	MappedGrid & mg = cg[grid];
	  	intArray & mask = mg.mask();
		  OV_GET_SERIAL_ARRAY(int,mask,maskLocal);

		  getIndex(mg.gridIndexRange(),I1,I2,I3);
		  bool ok=ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3);  

		  int na=0; // counts active points
	    if( ok )
	    {
			  IntegerArray & index = indexVector[grid];

			  index.redim(maskLocal.dimension(0),maskLocal.dimension(1),maskLocal.dimension(2)); 	
			  index=-1; // means in-active

		    FOR_3D(i1,i2,i3,I1,I2,I3)
		    {
		    	if( maskLocal(i1,i2,i3)>0 )
		    	{
		    		index(i1,i2,i3)=na; na++;  // active points 
		    	}

		    }
		  }

	    numActive(grid,myid)=na;

	    printf("testGlobalIndex: myid=%d: grid=%d numActive=%d\n",myid,grid,numActive(grid,myid) );


	  }

	  // -- transfer info to all processors ---
	  // -- for now just sum entries since only 1 is nonzero -- thgere should be a better way

	  IntegerArray numActiveBuffer(numberOfComponentGrids,np); numActiveBuffer=0; // save sums here
	  int numberOfSums= numberOfComponentGrids*np;
	  ParallelUtility::getSums( numActive.getDataPointer(),numActiveBuffer.getDataPointer(),numberOfSums );
	  numActive=numActiveBuffer;

	  Range all;
	  int numActiveLocal = sum(numActive(all,myid)); // number of active points on this processor

	  int totalActive = sum(numActive);

	  printf("myid=%d:  numActiveLocal=%d\n",myid,numActiveLocal);

	  printF("Summary: totalActivePoints=%d\n Active points per processor:\n",totalActive);
	  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	  {
	  	printF("grid=%d: numActive=[",grid);
	  	for( int p=0; p<np; p++ )
	  		printF("%6d,",numActive(grid,p));
	  	printF("]\n");
	  }
	  
	  // ---- compute the global offsets ----

	  IntegerArray globalOffset(numberOfComponentGrids); // offset from local to global index
	  // ------------------------------------------------------------------
	  // global index is ordered first by processor then by grid: 
	  //   [ proc=0 grid=0 ] [proc=0 grid=1] ... [proc=0,numGrids-1] [ proc=1 grid=0 ] [proc=1 grid=1] ... [proc=1,numGrids-1]  ...
	  //    
	  // offset to (myid,grid) = SUM { all previous processors, all grids}  + SUM { proc=mydi and all previous grids }
	  //----------------------------------------------------------------------
	  int baseOffsetToThisProcessor=0;
	  for( int p=0; p<myid; p++ )
		  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	      baseOffsetToThisProcessor += numActive(grid,p); 

	  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	  {
	    globalOffset(grid) = baseOffsetToThisProcessor;
			for( int gridp=0; gridp<grid; gridp++ )
			  globalOffset(grid) += numActive(gridp,myid); 

	    printf("myid=%d: grid=%d globalOffset =%d\n",myid,grid,globalOffset(grid) );

	  }

	  // -- turn local index into the global index --
	  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	  {
	  	MappedGrid & mg = cg[grid];
	  	intArray & mask = mg.mask();
		  OV_GET_SERIAL_ARRAY(int,mask,maskLocal);

		  getIndex(mg.gridIndexRange(),I1,I2,I3);
		  bool ok=ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3);  

	    if( ok )
	    {
	    	const int offset = globalOffset(grid);
			  IntegerArray & index = indexVector[grid];
		    FOR_3D(i1,i2,i3,I1,I2,I3)
		    {
		    	if( index(i1,i2,i3)>=0 )
		    	{
		    		index(i1,i2,i3)+=offset;
		    	}

		    }
		  }
	  } // end for grid 


	  RealArray v(totalActive); // just make a serial arry with all points for testing 
	  for( int i=0; i<totalActive; i++ )
	  {
	  	v(i)=i;
	  }

	  realCompositeGridFunction u(cg);
	  u=0.;


	  // ----------------------------------
	  // ----- vector to grid function ----
	  // ----------------------------------
	  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	  {
	  	MappedGrid & mg = cg[grid];
	  	intArray & mask = mg.mask();
		  OV_GET_SERIAL_ARRAY(int,mask,maskLocal);

		  OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);	  

		  getIndex(mg.gridIndexRange(),I1,I2,I3);
		  bool ok=ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3);  

	    if( ok )
	    {
			  IntegerArray & index = indexVector[grid];
		    FOR_3D(i1,i2,i3,I1,I2,I3)
		    {
		    	if( index(i1,i2,i3)>=0 )
		    	{
		    		int ig = index(i1,i2,i3);   // global index
		    		assert( ig>=0 && ig<totalActive );

	          uLocal(i1,i2,i3) = v(ig);
		    	}

		    }
		  }  
		} // end for grid 


	  if( 1==0 )
		  u.display("u after vector to grid function","%5.0f ");

	  // ----------------------------------
	  // ----- grid function to vector ----
	  // ----------------------------------
	  RealArray v2(totalActive); // just make a serial arry with all points for testing 
	  v2=0;

	  for( int grid=0; grid<numberOfComponentGrids; grid++ )
	  {
	  	MappedGrid & mg = cg[grid];
	  	intArray & mask = mg.mask();
		  OV_GET_SERIAL_ARRAY(int,mask,maskLocal);

		  OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);	  

		  getIndex(mg.gridIndexRange(),I1,I2,I3);
		  bool ok=ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3);  

	    if( ok )
	    {
			  IntegerArray & index = indexVector[grid];
		    FOR_3D(i1,i2,i3,I1,I2,I3)
		    {
		    	if( index(i1,i2,i3)>=0 )
		    	{
		    		int ig = index(i1,i2,i3);   // global index
		    		assert( ig>=0 && ig<totalActive );

	          v2(ig) = uLocal(i1,i2,i3);
		    	}

		    }
		  }  
		} // end for grid  


	  // -- check result
	  RealArray v2Buffer(totalActive);
	  ParallelUtility::getSums( v2.getDataPointer(),v2Buffer.getDataPointer(),totalActive );
	  v2 = v2Buffer;

		Real maxErr = max(fabs(v-v2)); 
		printf("Max error after GF to VECTOR and reverse = %9.2e\n",maxErr);


	  delete [] indexVector;

	}

  Overture::finish();          
  return 0;
}
