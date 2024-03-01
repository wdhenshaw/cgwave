#include "CgWave.h"
#include "CompositeGridOperators.h"   
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
  CompositeGridOperators & operators        = dbase.get<CompositeGridOperators>("operators");
  const real & cfl                          = dbase.get<real>("cfl");
  real & dt                                 = dbase.get<real>("dt");
  const real & ad4                          = dbase.get<real>("ad4"); // coeff of the artificial dissipation.
  const Real & dtMax                        = dbase.get<Real>("dtMax"); 
  const real & c                            = dbase.get<real>("c");
  const int & orderOfAccuracyInTime         = dbase.get<int>("orderOfAccuracyInTime");
  const int & orderOfAccuracy               = dbase.get<int>("orderOfAccuracy");
  IntegerArray & gridIsImplicit             = dbase.get<IntegerArray>("gridIsImplicit");
  const int & preComputeUpwindUt            = dbase.get<int>("preComputeUpwindUt");
  // const int & chooseImplicitTimeStepFromCFL = dbase.get<int>("chooseImplicitTimeStepFromCFL");
  const int & chooseTimeStepFromExplicitGrids= dbase.get<int>("chooseTimeStepFromExplicitGrids"); 

  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

  const int useUpwindDissipation = ad4>0.;  // ** FIX ME**

  printF("CgWave::getTimeStep: c=%g, cfl=%g, timeSteppingMethod=%d\n",c,cfl,(int)timeSteppingMethod);

  const bool allGridsAreImplicit = min(gridIsImplicit);

  Real stabilityBound=1.; 
  if( orderOfAccuracyInTime==2 )
  {

    // if( orderOfAccuracy==2 && useUpwindDissipation )
    // {
    //   stabilityBound=.7; // *************** TEMPORARY: FD22s needs cfl=.7 for some reason ... FIX ME 
    // }

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
  }

  dt              = REAL_MAX;   // actual time-step
  Real dtExplicit = REAL_MAX;   // dtExplicit : time-step when all grids are treated explicitly
  
  int numberOfDimensions = cg.numberOfDimensions();

  // Holds min/max grid spacing for each grid 
  RealArray & dxMinMax = dbase.get<RealArray>("dxMinMax");
  dxMinMax.redim(cg.numberOfComponentGrids(),2);

  // -------------------------------------------------------------------
  // gridCFL(grid) =  "c/dx" = effective CFL/dt number for each grid 
  // NOTE: We store this variable since "dt" may be changed later
  // -------------------------------------------------------------------

  if( !dbase.has_key("gridCFL") )
  {
    RealArray & gridCFL = dbase.put<RealArray>("gridCFL");
    gridCFL.redim(cg.numberOfComponentGrids());
    gridCFL=-1.;

  }
  RealArray & gridCFL = dbase.get<RealArray>("gridCFL");

  // gridSize(grid) = "dx" used in computing dt *new* Jan 4, 2024
  if( !dbase.has_key("gridSize") )
  {
    RealArray & gridSize = dbase.put<RealArray>("gridSize");
    gridSize.redim(cg.numberOfComponentGrids());
    gridSize=-1.;
  }
  RealArray & gridSize = dbase.get<RealArray>("gridSize");

  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];

    Real dtGrid = REAL_MAX;  // max dt for this grid 
    
    real dx[3]={1.,1.,1.};
    if( mg.isRectangular() )
    {
      mg.getDeltaX(dx);
      

      if( mg.numberOfDimensions()==2 )
      {
        gridSize(grid) = 1./( sqrt( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) ) );

        dtGrid = cfl*1./( c*sqrt( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) ) );  
      }
      else
      {
        gridSize(grid) = 1./( sqrt( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) + 1./(dx[2]*dx[2]) ) );

        dtGrid = cfl*1./( c*sqrt( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) + 1./(dx[2]*dx[2]) ) );
      }

      printF("getTimeStep: grid=%d : dx=%9.2e, dy=%9.2e, gridSize=%9.2e, dt=%9.3e\n",grid,dx[0],dx[1],gridSize(grid),dtGrid);

      dxMinMax(grid,0)= numberOfDimensions == 2? min(dx[0],dx[1]) : min(dx[0],dx[1],dx[2]);
      dxMinMax(grid,1)= numberOfDimensions == 2? min(dx[0],dx[1]) : max(dx[0],dx[1],dx[2]);

      gridCFL(grid)= 1./(dtGrid/cfl);  // = "c/dx"
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

        gridSize(grid) = dxMin;        

        dtGrid = (cfl/c) * dxMin;

        gridCFL(grid)= 1./(dtGrid/cfl); // "c/dx"

        printF("getTimeStep: grid=%d, dxMin=%9.2e, dxMax=%9.2e, gridSize=%9.2e, dt=%9.3e\n",grid,dxMin,dxMax,gridSize(grid),dtGrid);        
      }
        
    } // end curvilinear grid

    if( orderOfAccuracyInTime==2 )
    {
      // ----- order in time =2 may have a smaller stability bound ----
      dtGrid *= stabilityBound; 
      // gridCFL(grid) *= stabilityBound;
    }    

    dtExplicit = min(dtExplicit,dtGrid);

    if( !gridIsImplicit(grid) )
      dt = min(dt,dtGrid); 
    
  } // end for grid 

  if( true )
  {
    // ********** NEW WAY : Jan 4, 2024 **************

    // Given gridSize(grid) compute the time-step dt

    Real dtExplicitGrids = dtMax>0 ? dtMax : REAL_MAX;  
    Real dtImplicitGrids = dtMax>0 ? dtMax : REAL_MAX;

    // ----- order-in-time=2 + order-in-space > 2 may have a smaller stability bound ----
    const Real scaledCFL = orderOfAccuracyInTime==2 ? cfl*stabilityBound : cfl;

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      // c*dt/dx = cfl : 
      Real dtGrid = (scaledCFL/c) * gridSize(grid);
      if( gridIsImplicit(grid) )
        dtImplicitGrids = min( dtImplicitGrids, dtGrid );
      else
        dtExplicitGrids = min( dtExplicitGrids, dtGrid );
    }
    dtExplicitGrids = ParallelUtility::getMinValue(dtExplicitGrids);  // get min value over all processors 
    dtImplicitGrids = ParallelUtility::getMinValue(dtImplicitGrids);  // get min value over all processors 

    if( timeSteppingMethod==explicitTimeStepping )
    {
      dt = dtExplicitGrids;
    }
    else 
    { // -- implicit time-stepping ---

      dt = dtExplicitGrids; // by default choose dt from any explicit grids (or dtMax if set)

      if( chooseTimeStepFromExplicitGrids ) 
      {
        // --- choose the implicit dt based on the CFL parameter (which can be large) ----
        if( allGridsAreImplicit )
        {
          dt=dtImplicitGrids;
          printF("\n >>>>chooseTimeStep: implicit-time-stepping: all grids implicit: choose c*dt/dxMin = cfl, or dt=dtMax: setting dt=%9.2e, dtMax=%9.2e, cfl=%g\n\n",dt,dtMax,cfl);         
        }
        else
        {
          const Real dtRatio = dtExplicitGrids/min(dtExplicitGrids,dtImplicitGrids);
          printF("\n >>>>chooseTimeStep: implicit-time-stepping: choose c*dt/dxMin from explicit grids only, or dtMax: setting dt=%9.2e, dtMax=%9.2e cfl=%g, "
                 " dtRatio=dtExplicit/dtAllGrids=%5.2f\n\n",dt,dtMax,cfl,dtRatio);
        }
      }

    }
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ ) 
    {
      gridCFL(grid)= c/gridSize(grid);
      printF("getTimeStep: grid=%d: gridSize=%10.2e, gridCFL*dt=c*dt/gridSize=%10.2e\n",grid,gridSize(grid),gridCFL(grid)*dt);
    }
      

  }
  else
  {
    // ******** OLD WAY 

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
      printF("cgWave::getTimeStep: cfl=%g, explicit dt=%9.2e, implicit dt=%9.2e, ratio=%6.2f\n",
           cfl,dtExplicit,dt,dt/dtExplicit);

      if( chooseTimeStepFromExplicitGrids ) // && !allGridsAreImplicit )
      {
        // --- choose the implicit dt based on the CFL parameter (which can be large) ----
        if( allGridsAreImplicit )
        {
           dt=dtExplicit;
          printF("\n >>>>cgWave:chooseTimeStep: choose dt from CFL and explicit time-step: setting dt=%9.2e, cfl=%g\n\n",dt,cfl);         
        }
        else
        {
          const Real dtRatio = dt/dtExplicit;
          dt=dtExplicit;

          printF("\n >>>>cgWave:chooseTimeStep: choose dt from CFL and explicit grids only: setting dt=%9.2e, cfl=%g, dtExplicit/dtAllGrids=%5.2f\n\n",
            dt,cfl,dtRatio);
        }
      }
    }

    // Save the grid CFL number (used for upwind dissipation coefficient in advWave.bf90 )
    // gridCFL = cfl*( dt/gridCFL );   //  **check me : should be adjust this if dtMax is chosen ?? ++++++++++++++++++++=
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      printF("getTimeStep: grid=%d: gridCFL = c/dx =%10.2e, gridCFL*dt=%10.2e\n",grid,gridCFL(grid),gridCFL(grid)*dt);
    }

  }

  return 0;
}
