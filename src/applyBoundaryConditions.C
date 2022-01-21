// This file automatically generated from applyBoundaryConditions.bC with bpp.
#include "CgWave.h"
#include "CompositeGridOperators.h";    
#include "PlotStuff.h"
#include "display.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"

// -------- function prototypes for Fortran routines --------
#define bcOptWave EXTERN_C_NAME(bcoptwave)
extern "C"
{

void bcOptWave( const int&nd, 
                                const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                                const int&gridIndexRange, const int& dimRange, const int &isPeriodic, real&u, const int&mask,
                                const real&rsxy, const real&xy, const int&boundaryCondition, const real & frequencyArray,
                                const DataBase *pdb, const int&ipar, const real&rpar, int&ierr );

}

#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )

// The getBcOptParameters macro is defined here:
// =======================================================================================
// Macro: get the parameters for calling the optimized fortran BC routine
// =======================================================================================


// ======================================================================================================
/// \brief apply boundary conditions
// ======================================================================================================
int CgWave::
applyBoundaryConditions( realCompositeGridFunction & u, real t )
{
    real cpu0=getCPU();

    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::numberOfProcessors());

    const int & debug = dbase.get<int>("debug");

    if( debug & 8 )
        printF("applyBoundaryConditions at t=%9.3e\n",t);

  // realCompositeGridFunction *& ucg = dbase.get<realCompositeGridFunction*>("ucg");
  // realCompositeGridFunction & u = ucg[current];

    const int & orderOfAccuracy   = dbase.get<int>("orderOfAccuracy");
    const Real & c                = dbase.get<real>("c");
    const real & dt               = dbase.get<real>("dt");
  // const real & ad4              = dbase.get<real>("ad4"); // coeff of the artificial dissipation.
  // bool useUpwindDissipation     = ad4  > 0.;
    const int & upwind                = dbase.get<int>("upwind");
    bool useUpwindDissipation     = upwind!=0;
    const int & solveHelmholtz    = dbase.get<int>("solveHelmholtz");
    IntegerArray & gridIsImplicit = dbase.get<IntegerArray>("gridIsImplicit");

    const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption"); 

    const int & applyKnownSolutionAtBoundaries = dbase.get<int>("applyKnownSolutionAtBoundaries"); // by default, do NOT apply known solution at boundaries

    const int & addForcing                  = dbase.get<int>("addForcing");
    const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const bool twilightZone = forcingOption==twilightZoneForcing; 

    const int & numberOfFrequencies     = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray    = dbase.get<RealArray>("frequencyArray");  

    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");

    BoundaryConditionParameters bcParams;
    bcParams.orderOfExtrapolation=orderOfAccuracy+1;

  // if( false )
  // {
  //   // Is this needed ?
  //   if( cg.numberOfComponentGrids()>1 )
  //   {
  //     for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  //     {
  //    u[grid].updateGhostBoundaries();
  //     }
  //   } 
  // }
    

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];

    cpu0 = getCPU();

  // *** Note: interpolate also does a periodic update **** 
    u.interpolate();

    real cpu1 = getCPU();
    timing(timeForInterpolate)+= cpu1-cpu0;
    
    cpu0=cpu1;

  // ---- check for valid boundary conditions ---
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        ForBoundary(side,axis) 
        {
            int bc = mg.boundaryCondition(side,axis);
            if( !( bc<=0 || bc==dirichlet || bc==evenSymmetry || bc==neumann  ) )
            {
                printF("CgWave:applyBoundaryConditions:ERROR: grid=%i side=%i axis=%i bc=%i not implemented\n",grid,side,axis,bc);
                printF(" You should set boundary conditions in the interactiveUpdate \n");
                OV_ABORT("ERROR");

                printF(" Setting bc=dirichlet\n");
                mg.setBoundaryCondition(side,axis,dirichlet);
            }
        }
    }

    int numGhost = orderOfAccuracy/2;
    if( useUpwindDissipation ) numGhost++;
    const int assignBCForImplicit=0; 

    bool useOpt= true; // twilightZone && true;
    
    if( debug & 4 && t<=2*dt )
        printF("CgWave: applyBC: useOpt=%d\n",(int)useOpt);
    if( useOpt )
    {
    // ----- Optimized Boundary Conditions ------
        bool isImplicit    = false; // at least one grid is implicit
        bool isAllImplicit = true;  // all grids are implicit
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            isImplicit    = isImplicit || gridIsImplicit(grid); 
            isAllImplicit = isAllImplicit && gridIsImplicit(grid); 
      // +++ NOTE: BCs for implicit are done in takeImplicitTimeStep +++
            if( !gridIsImplicit(grid) )
            {
                MappedGrid & mg = cg[grid];
        // const IntegerArray & gid = mg.gridIndexRange();
                
                OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

        // get parameters for calling fortran
                    IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
                    ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( u[grid],indexRangeLocal,dimLocal,bcLocal );
                    const bool isRectangular=mg.isRectangular();
                    real dx[3]={1.,1.,1.};
                    if( isRectangular )
                        mg.getDeltaX(dx);
          // int assignKnownSolutionAtBoundaries = 0;  // changed below 
                    DataBase *pdb = &dbase;
          // Real cfl1 = pdb->get<real>("cfl");
          // printF(" CFL from pdb: cfl1=%g\n",cfl1);
                    int knownSolutionOption=0; // no known solution
                    if( dbase.has_key("userDefinedKnownSolutionData") )
                    {
                        DataBase & db =  dbase.get<DataBase>("userDefinedKnownSolutionData");
                        const aString & userKnownSolution = db.get<aString>("userKnownSolution");
                        if( userKnownSolution=="planeWave"  )
                        {
                            knownSolutionOption=1;                   // this number must match in bcOptWave.bf90
              // assignKnownSolutionAtBoundaries=1;
                        }
                        else if( userKnownSolution=="gaussianPlaneWave"  ) 
                        {
                            knownSolutionOption=2;                   // this number must match in bcOptWave.bf90
              // assignKnownSolutionAtBoundaries=1;  
                        }    
                        else if( userKnownSolution=="boxHelmholtz"  ) 
                        {
                            knownSolutionOption=3;                   // this number must match in bcOptWave.bf90
              // assignKnownSolutionAtBoundaries=1;  // not needed for square or box but is needed for cic **fix me**
                        }
                        else if( userKnownSolution=="polyPeriodic"  ) 
                        {
                            knownSolutionOption=4;                   // this number must match in bcOptWave.bf90
              // assignKnownSolutionAtBoundaries=1;  
                        }    
                    }
                    int gridType = isRectangular ? 0 : 1;
                    int gridIsImplicit = 0; 
                    int numberOfComponents = 1;          // for now CgWave only has a single component 
                    int uc = 0;                          // first component
                    int ipar[] = {
                        uc                  ,            // ipar( 0)
                        numberOfComponents  ,            // ipar( 1)
                        grid                ,            // ipar( 2)
                        gridType            ,            // ipar( 3)
                        orderOfAccuracy     ,            // ipar( 4)
                        gridIsImplicit      ,            // ipar( 5)
                        twilightZone        ,            // ipar( 6)
                        np                  ,            // ipar( 7)
                        debug               ,            // ipar( 8)
                        myid                ,            // ipar( 9)
                        applyKnownSolutionAtBoundaries,  // ipar(10)
                        knownSolutionOption,             // ipar(11)
                        addForcing,                      // ipar(12)
                        forcingOption,                   // ipar(13)
                        useUpwindDissipation,            // ipar(14)
                        numGhost,                        // ipar(15)
                        assignBCForImplicit,             // ipar(16)
                        bcApproach,                      // ipar(17)
                        numberOfFrequencies              // ipar(18)
                                              };
                    real rpar[] = {
                        t                , //  rpar( 0)
                        dt               , //  rpar( 1)
                        dx[0]            , //  rpar( 2)
                        dx[1]            , //  rpar( 3)
                        dx[2]            , //  rpar( 4)
                        mg.gridSpacing(0), //  rpar( 5)
                        mg.gridSpacing(1), //  rpar( 6)
                        mg.gridSpacing(2), //  rpar( 7)
                        (real &)(dbase.get<OGFunction* >("tz")) ,        //  rpar( 8) ! pointer for exact solution -- new : 110311 
                        REAL_MIN,         //  rpar( 9)
                        c                 //  rpar(10)
                                                };
                    real *pu = uLocal.getDataPointer();
                    int *pmask = maskLocal.getDataPointer();
                    real temp, *pxy=&temp, *prsxy=&temp;
                    if( !isRectangular )
                    {
                        mg.update(MappedGrid::THEinverseVertexDerivative);
                        #ifdef USE_PPP
                          prsxy=mg.inverseVertexDerivative().getLocalArray().getDataPointer();
                        #else
                          prsxy=mg.inverseVertexDerivative().getDataPointer();
                        #endif    
                    }
                    bool vertexNeeded = twilightZone || knownSolutionOption!=0;
                    if( vertexNeeded )
                    {
                        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
                        #ifdef USE_PPP
                          pxy=mg.vertex().getLocalArray().getDataPointer();
                        #else
                          pxy=mg.vertex().getDataPointer();
                        #endif    
                    }

                int ierr=0;
                bcOptWave(mg.numberOfDimensions(),
                                    uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
                                    uLocal.getBase(2),uLocal.getBound(2),
                                    indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),
                                    *pu, *pmask, *prsxy, *pxy,  bcLocal(0,0), frequencyArray(0),
                                    pdb, ipar[0],rpar[0], ierr );

        // ...swap periodic edges 
                u[grid].periodicUpdate();
                u[grid].updateGhostBoundaries();
            }

        } // end for grid 

        if( useUpwindDissipation )
        {
      // For upwind dissipation we should assign another line of points next to interpolation points
      // to support the wider upwind stencil

            if( !isAllImplicit )
                u.applyBoundaryCondition(0,BCTypes::extrapolateInterpolationNeighbours,BCTypes::allBoundaries,0.,t,bcParams);
            else
            {
                if( isImplicit && !isAllImplicit )
                {
                    printF("applyBoundaryConditions: t=%9.3e: FIX ME: extrapolateInterpolationNeighbours for PARTIALY IMPLICIT time-stepping\n"); 
                    OV_ABORT("error");
                }
            }
        }

    // ** TESTING ***
        if( false )
        {
            printF("applyBC: TESTING: SET UNUSED GHOST TO ZERO\n");
    
            for( int ghost=2; ghost<=2; ghost++ )
            {
                bcParams.ghostLineToAssign=ghost;
                u.applyBoundaryCondition(0,BCTypes::extrapolate,dirichlet,0.,t,bcParams);
            }
            bcParams.ghostLineToAssign=1; // reset 
        }


        return 0; 

    } // end for useOpt

    else
    {
    // ==== OLD WAY ====

        if( knownSolutionOption=="userDefinedKnownSolution" ) 
      //  !solveHelmholtz )  // when solving Helmholtz do not use exact solution for BC's
        {
            
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const IntegerArray & gid = mg.gridIndexRange();
                
                realArray & ug = u[grid];
                ForBoundary(side,axis) 
                {
                    if( mg.boundaryCondition(side,axis)>0 )
                    {
                        getBoundaryIndex(gid,side,axis,I1,I2,I3);
            // int numGhost=orderOfAccuracy/2;
            // if( useUpwindDissipation ) numGhost++;
            // if( side==0 )
            //   Iv[axis] = Range(gid(0,axis)-numGhost,gid(0,axis));
            // else 
            //   Iv[axis] = Range(gid(1,axis),gid(1,axis)+numGhost);

                        getUserDefinedKnownSolution( t, grid, ug, I1,I2,I3 );
                    }
                }
            }
                          
        }
        else  
        {


      // -- dirichlet --
            u.applyBoundaryCondition(0,BCTypes::dirichlet,dirichlet,0.,t);

        }
        
    // -- Extrapolate ghost outside Dirichlet boundaries
    //  -- COULD DO BETTER WITH COMPATIBILITY ---  *** FIX ME ***
    // Compatibility:
    //      u = g  at x=a 
    //      c^2( u_xx + g_yy ) + f = g_tt 
    // 
        int numGhost = orderOfAccuracy/2;
        if( useUpwindDissipation ) numGhost++;
        for( int ghost=1; ghost<=numGhost; ghost++ )
        {
            bcParams.ghostLineToAssign=ghost;
            u.applyBoundaryCondition(0,BCTypes::extrapolate,dirichlet,0.,t,bcParams);
        }
        bcParams.ghostLineToAssign=1; // reset 
            
            
    // bcParams.lineToAssign=1;  // set first ghost line 
    // bcParams.extraInTangentialDirections=1;
    // u.applyBoundaryCondition(0,BCTypes::dirichlet,dirichlet,0.,t,bcParams);
        if( false )
        {
            if(  orderOfAccuracy>=4 || (orderOfAccuracy>=2 && useUpwindDissipation) )
            {
                bcParams.lineToAssign=2;  // set second ghost line 
                bcParams.extraInTangentialDirections=2;
                u.applyBoundaryCondition(0,BCTypes::dirichlet,dirichlet,0.,t,bcParams);
            }
            if( orderOfAccuracy>=4 && useUpwindDissipation )
            {
                bcParams.lineToAssign=3;  // set second ghost line 
                bcParams.extraInTangentialDirections=2;
                u.applyBoundaryCondition(0,BCTypes::dirichlet,dirichlet,0.,t,bcParams);
            }
        }
            

    // -- even symmetry ---

        u.applyBoundaryCondition(0,BCTypes::evenSymmetry,evenSymmetry,0.,t);
        bcParams.ghostLineToAssign=2;
        u.applyBoundaryCondition(0,BCTypes::evenSymmetry,evenSymmetry,0.,t,bcParams);
        if( orderOfAccuracy>=4 && useUpwindDissipation )
        {
            bcParams.ghostLineToAssign=3;
            u.applyBoundaryCondition(0,BCTypes::evenSymmetry,evenSymmetry,0.,t,bcParams);
        }
            
        
    // BoundaryConditionParameters extrapParams;
    // extrapParams.orderOfExtrapolation=orderOfAccuracy+1

        bcParams.ghostLineToAssign=1; // reset 

    } // end old way


    if( useUpwindDissipation )
    {
    // For upwind dissipation we should assign another line of points next to interpolation points
    // to support the wider upwind stencil

          // extrapParams.orderOfExtrapolation=dbase.get<int>("orderOfExtrapolationForInterpolationNeighbours");
          // if( debug & 4 )
          //   printf("***advanceStructured: orderOfExtrapolationForInterpolationNeighbours=%i\n",
          //          dbase.get<int>("orderOfExtrapolationForInterpolationNeighbours"));
        
          // // MappedGridOperators & mgop = mgp!=NULL ? *op : (*cgop)[grid];
          // // fieldNext.setOperators(mgop);
          // fieldNext.applyBoundaryCondition(Ca,BCTypes::extrapolateInterpolationNeighbours,BCTypes::allBoundaries,0.,t+dt,
          //                                  extrapParams,grid);

        u.applyBoundaryCondition(0,BCTypes::extrapolateInterpolationNeighbours,BCTypes::allBoundaries,0.,t,bcParams);
    }
    
/* -----
      if( boundaryCondition!=BCTypes::evenSymmetry && (orderOfAccuracy==4 || ad4!=0.) )
      { // extrapolate 2nd ghostline for 4th order when the grid only supports second order
      if( secondOrderGrid )
      {
      bcParams.ghostLineToAssign=2;
      bcParams.orderOfExtrapolation=4;
      u.applyBoundaryCondition(0,BCTypes::extrapolate,BCTypes::allBoundaries,0.,t,bcParams);
   // also extrapolate unused points next to interpolation points -- this allows us to
   // avoid making a grid with 2 lines of interpolation.
      u.applyBoundaryCondition(0,BCTypes::extrapolateInterpolationNeighbours);
      }
      }
      ---- */
        
  // ------ Zero out unused ghost points -----
  // May be needed for CgWaveHoltz and PETSc solver 
  // Do we need to use Overture's sero put unused to get other unused points ??
    if( false )
    {
    // ** BROKEN IN PARALLEL ***
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();
            const IntegerArray & dim = mg.dimension();
        
            ForBoundary(side,axis) 
            {
                getBoundaryIndex(dim,side,axis,I1,I2,I3);

                int ia=side==0 ? dim(side,axis)          : gid(side,axis)+numGhost;
                int ib=side==0 ? gid(side,axis)-numGhost : dim(side,axis);
                if( ib >= ia )
                {
                    Iv[axis] = Range(ia,ib);
                    OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                    bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
                    if( ok )
                    {
                        uLocal(I1,I2,I3)=0.;
                    }
                }
            }
        }
    }
    
    if( debug & 8 )
        printF("...BC's, now finishBoundaryConditions...\n");


  // *** check -- does finishBoundaryConditions update parallel ghost ??
    BoundaryConditionParameters extrapParams;
    extrapParams.orderOfExtrapolation=orderOfAccuracy+1;
    u.finishBoundaryConditions(extrapParams);

  // // check this: 
  // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  //   u[grid].updateGhostBoundaries();

    timing(timeForBoundaryConditions) += getCPU()-cpu0;
        
    return 0;
}

