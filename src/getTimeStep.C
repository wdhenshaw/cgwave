#include "CgWave.h"
#include "CompositeGridOperators.h";    
#include "ParallelUtility.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  

#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )

// ================================================================================================
/// \brief Determine the time step
// ================================================================================================
int CgWave::
getTimeStep()
{
  CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
  const real & cfl                   = dbase.get<real>("cfl");
  real & dt                          = dbase.get<real>("dt");
  const Real & dtMax                 = dbase.get<Real>("dtMax"); 
  const real & c                     = dbase.get<real>("c");
  const int & orderOfAccuracyInTime  = dbase.get<int>("orderOfAccuracyInTime");
  const int & orderOfAccuracy        = dbase.get<int>("orderOfAccuracy");
  IntegerArray & gridIsImplicit      = dbase.get<IntegerArray>("gridIsImplicit");

  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

  printF("CgWave::getTimeStep: c=%g, cfl=%g, timeSteppingMethod=%d\n",c,cfl,(int)timeSteppingMethod);



  dt              = REAL_MAX;   // actual time-step
  Real dtExplicit = REAL_MAX;   // dtExplicit : time-step when all grids are treated explicitly
  
  int numberOfDimensions = cg.numberOfDimensions();

  // Holds min/max grid spacing for eacg grid 
  RealArray & dxMinMax = dbase.get<RealArray>("dxMinMax");
  dxMinMax.redim(cg.numberOfComponentGrids(),2);

  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];

    Real dtGrid = REAL_MAX;  // max dt for this grid 
    
    real dx[3]={1.,1.,1.};
    if( mg.isRectangular() )
    {
      mg.getDeltaX(dx);
      

      if( mg.numberOfDimensions()==2 )
        dtGrid = cfl*1./( c*sqrt( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) ) );  
      else
        dtGrid = cfl*1./( c*sqrt( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) + 1./(dx[2]*dx[2]) ) );

      printF("getTimeStep: grid=%d : dx=%9.2e, dy=%9.2e, dt=%9.3e\n",grid,dx[0],dx[1],dtGrid);

      dxMinMax(grid,0)= numberOfDimensions == 2? min(dx[0],dx[1]) : min(dx[0],dx[1],dx[2]);
      dxMinMax(grid,1)= numberOfDimensions == 2? min(dx[0],dx[1]) : max(dx[0],dx[1],dx[2]);
    }
    else
    {
      // ---- curvilinear grid ----

      mg.update(MappedGrid::THEinverseVertexDerivative);
      const realArray & rx = mg.inverseVertexDerivative();
      const intArray & mask = mg.mask();
        
      
      Index I1,I2,I3;
      getIndex( mg.indexRange(),I1,I2,I3);

      // Grid spacings on unit square:
      real dr1 = mg.gridSpacing(axis1);
      real dr2 = mg.gridSpacing(axis2);
      real dr3 = mg.gridSpacing(axis3);

      // parallel version here --- also broadcast max error in forcing.bC *************************
      OV_GET_SERIAL_ARRAY_CONST(real,rx,rxLocal);
      OV_GET_SERIAL_ARRAY_CONST(int,mask,maskLocal);


      real *rxp = rxLocal.Array_Descriptor.Array_View_Pointer3;
      const int rxDim0=rxLocal.getRawDataSize(0);
      const int rxDim1=rxLocal.getRawDataSize(1);
      const int rxDim2=rxLocal.getRawDataSize(2);
      const int rxDim3=mg.numberOfDimensions();   // note
      #undef RX
      #define RX(i0,i1,i2,i3,i4) rxp[i0+rxDim0*(i1+rxDim1*(i2+rxDim2*(i3+rxDim3*(i4))))]

      const int *maskp = maskLocal.Array_Descriptor.Array_View_Pointer2;
      const int maskDim0=maskLocal.getRawDataSize(0);
      const int maskDim1=maskLocal.getRawDataSize(1);
      const int md1=maskDim0, md2=md1*maskDim1; 
      #define MASK(i0,i1,i2) maskp[(i0)+(i1)*md1+(i2)*md2]


      int includeGhost=0;
      bool ok = ParallelUtility::getLocalArrayBounds(rx,rxLocal,I1,I2,I3,includeGhost);

      int i1,i2,i3;
      real a11Min=REAL_MAX*.001;
      real a11Max=-a11Min;
      //  **** this is a guess **** check this.
      const real alpha0=1.;

      real a11,a12,a22;
      if( ok )
      {
        if( numberOfDimensions==2 )
        {

          FOR_3D(i1,i2,i3,I1,I2,I3)
          {
              
            if( MASK(i1,i2,i3)>0 )
            {
              a11 = ( RX(i1,i2,i3,0,0)*RX(i1,i2,i3,0,0) + RX(i1,i2,i3,0,1)*RX(i1,i2,i3,0,1) );
              a12 = ( RX(i1,i2,i3,0,0)*RX(i1,i2,i3,1,0) + RX(i1,i2,i3,0,1)*RX(i1,i2,i3,1,1) )*2.;
              a22 = ( RX(i1,i2,i3,1,0)*RX(i1,i2,i3,1,0) + RX(i1,i2,i3,1,1)*RX(i1,i2,i3,1,1) );

              // we could save work by delaying the sqrt to after the loop
              a11=1./sqrt( a11 *(1./(alpha0*dr1*dr1)) 
                           +abs(a12)*(.25/(alpha0*dr1*dr2))
                           +a22 *(1./(alpha0*dr2*dr2)) 
                );

              a11Min=min(a11Min,a11);
              a11Max=max(a11Max,a11);
          
            }
          }
        }
        else
        { // ***** 3D ********

          #define rxDotRx(axis,dir) (RX(i1,i2,i3,axis,0)*RX(i1,i2,i3,dir,0) \
                         + RX(i1,i2,i3,axis,1)*RX(i1,i2,i3,dir,1) \
                         + RX(i1,i2,i3,axis,2)*RX(i1,i2,i3,dir,2))
      
          // There would be a factor of 4 for the worst case plus/minus wave but we also
          // divide by a factor of 4 for the 2nd-order time stepping.
          FOR_3D(i1,i2,i3,I1,I2,I3)
          {
            if( MASK(i1,i2,i3)>0 )
            {
              // we could save work by delaying the sqrt to after the loop
              a11=1./sqrt(   rxDotRx(0,0) *(1./(dr1*dr1)) 
                             +rxDotRx(1,1) *(1./(dr2*dr2))
                             +rxDotRx(2,2) *(1./(dr3*dr3))
                             +abs(rxDotRx(1,0))*(.5/(dr2*dr1))  
                             +abs(rxDotRx(2,0))*(.5/(dr3*dr1)) 
                             +abs(rxDotRx(2,1))*(.5/(dr3*dr2)) );

              // ** a11 =  pow(a11,-.5);
                
              a11Min=min(a11Min,a11);
              a11Max=max(a11Max,a11);

            }
          }
            
          #undef rxDotRx
        }
        
        real dxMin=a11Min; 
        real dxMax=a11Max;
        dxMinMax(grid,0)= dxMin;
        dxMinMax(grid,1)= dxMax;

        dtGrid = (cfl/c) * dxMin;

        printF("getTimeStep: grid=%d, dxMin=%9.2e, dxMax=%9.2e, dt=%9.3e\n",grid,dxMin,dxMax,dtGrid);        
      }
        
    }

    if( orderOfAccuracyInTime==2 )
    {
      Real stabilityBound=1.; 

      if( orderOfAccuracy==4  )
      {
        // Scheme FD24 has a smaller stability region
        stabilityBound = sqrt(3)/2.; // = .8660... check me 
      }
      else if( orderOfAccuracy==6 )
      {
        stabilityBound = .6; // FIX ME -- JUST A GUESS
      }
      else if( orderOfAccuracy==8 )
      {
        stabilityBound = .4; // FIX ME -- JUST A GUESS
      }      
      else if( orderOfAccuracy!=2 )
      {
        printF("cgWave::getTimeStep: ERROR: unexpected orderOfAccuracy=%d\n",orderOfAccuracy);
        OV_ABORT("ERROR");
      }
      dtGrid *= stabilityBound; 

    }    

    dtExplicit = min(dtExplicit,dtGrid);

    if( !gridIsImplicit(grid) )
      dt = min(dt,dtGrid); 
    
  } // end for grid 

  // Now enforce time-step to be at most dtMax 
  if( dtMax>0 )
  {
    dtExplicit = min(dtExplicit,dtMax);
    dt         = min(dt,        dtMax);
  }
  
  dtExplicit = ParallelUtility::getMinValue(dtExplicit);  // get min value over all processors 
  dt         = ParallelUtility::getMinValue(dt);          // get min value over all processors 


  if( timeSteppingMethod==implicitTimeStepping )
  {
    printF("cgWave::getTimeStep: explicit dt=%9.2e, implicit dt=%9.2e, ratio=%6.2f\n",dtExplicit,dt,dt/dtExplicit);
  }



  return 0;
}
