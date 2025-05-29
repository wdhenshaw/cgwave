// This file automatically generated from applyBoundaryConditions.bC with bpp.
#include "CgWave.h"
#include "CompositeGridOperators.h"    
#include "PlotStuff.h"
#include "display.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "AssignInterpNeighbours.h"

// -------- function prototypes for Fortran routines --------
#define bcOptWave EXTERN_C_NAME(bcoptwave)
#define abcWave EXTERN_C_NAME(abcwave)
extern "C"
{

void bcOptWave( const int&nd, 
                                const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                                const int&gridIndexRange, const int& dimRange, const int &isPeriodic, real&u, const real&un, const int&mask,
                                const real&rsxy, const real&xy, real &uTemp1, real & uTemp2, 
                                const int&boundaryCondition, const real & frequencyArray,
                                const DataBase *pdb, const int&ipar, const real&rpar, int&ierr );

void abcWave(const int&nd,
            const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
            const int&ndf1a,const int&ndf1b,const int&ndf2a,const int&ndf2b,const int&ndf3a,const int&ndf3b,
            const int & gid,
            const real&u, const real&un, const real&f, const int&mask, const real&rsxy, const real&xy,
            const int&bc, const int&boundaryCondition, const int&ipar, const real&rpar, int&ierr );

}

#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )

// The getBcOptParameters macro is defined here:
// =======================================================================================
// Macro: get the parameters for calling the optimized fortran BC routine
// =======================================================================================


// ============================================================================
// For upwind dissipation we assign another line of points next to interpolation points
// to support the wider upwind stencil
// ============================================================================


// ** FOR OLD WAY **


// ======================================================================================================
/// \brief Apply boundary conditions to a composite grid function.
/// 
/// \param u (input/output) : apply BCs to this grid function
/// \param un (input/output) : solution at the previous time
/// \param t (input) : current time
/// \param applyExplicitBoundaryConditions (input) : if true, apply explicit boundary conditions to all grids, even if
///    a grid is advanced implicitly. This can be used to set the BC's for initial conditions, for example.
/// \param fillImplicitBoundaryConditions (input) : if true, fill the right-hand-side for implicit time-stepping
///  into u.
// ======================================================================================================
int CgWave::
applyBoundaryConditions( realCompositeGridFunction & u, realCompositeGridFunction & un, real t,
                                                  bool applyExplicitBoundaryConditions /* = false */,
                                                  bool fillImplicitBoundaryConditions  /* = false */ )
{
    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::numberOfProcessors());

    const int & debug = dbase.get<int>("debug");

    dbase.get<int>("bcCount")++;    // count the number of times applyBC is called

  // if( applyExplicitBoundaryConditions )
  // {
  //   printF("**** DO NOT APPLY EXPLICIT BCS ----------- TEMP -------------\n");
  //   return 0;
  // }

    if( debug & 8 )
        printF("applyBoundaryConditions at t=%9.3e\n",t);

  // realCompositeGridFunction *& ucg = dbase.get<realCompositeGridFunction*>("ucg");
  // realCompositeGridFunction & u = ucg[current];

    const int & orderOfAccuracy      = dbase.get<int>("orderOfAccuracy");
    const Real & c                   = dbase.get<real>("c");
    const real & dt                  = dbase.get<real>("dt");
    const int & upwind               = dbase.get<int>("upwind");
    const int & implicitUpwind       = dbase.get<int>("implicitUpwind");
    const int & orderOfExtrapolation = dbase.get<int>("orderOfExtrapolation");

    bool useUpwindDissipation        = upwind!=0;

    const int & solveHelmholtz       = dbase.get<int>("solveHelmholtz");
    IntegerArray & gridIsImplicit    = dbase.get<IntegerArray>("gridIsImplicit");

    const aString & knownSolutionOption       = dbase.get<aString>("knownSolutionOption"); 
    const int & applyKnownSolutionAtBoundaries = dbase.get<int>("applyKnownSolutionAtBoundaries"); // by default, do NOT apply known solution at boundaries

    const int & addForcing                  = dbase.get<int>("addForcing");
    const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const bool twilightZone                 = forcingOption==twilightZoneForcing; 

    const int & numberOfFrequencies         = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");  
    const RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");  

    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");

    const AssignInterpolationNeighboursEnum & assignInterpNeighbours = 
                                                          dbase.get<AssignInterpolationNeighboursEnum>("assignInterpNeighbours");


    BoundaryConditionParameters bcParams;
  // bcParams.orderOfExtrapolation=orderOfAccuracy+1; *wdh* July 21, 2024
    bcParams.orderOfExtrapolation=orderOfExtrapolation;


    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];


    real cpu0 = getCPU();
    if( !fillImplicitBoundaryConditions )
    {
    // printF("applyBC: Call interpolate...\n");
    // *** Note: interpolate also does a periodic update **** 

        u.interpolate();
    }

    real cpu1 = getCPU();
    timing(timeForInterpolate)+= cpu1-cpu0;

  // printF("Done interpolate.\n");
    
    cpu0=cpu1;

  // ---- check for valid boundary conditions ---
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        ForBoundary(side,axis) 
        {
            int bc = mg.boundaryCondition(side,axis);
            if( !( bc<=0 || bc==dirichlet || bc==evenSymmetry || bc==neumann || bc==exactBC || bc==abcEM2 || bc==characteristic || bc==absorbing ) )
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

    const int assignBCForImplicit = fillImplicitBoundaryConditions;   // used in call to bcOptWave
    if( fillImplicitBoundaryConditions && t<=3.*dt && debug>3 )
    {
        printF("applyBoundaryConditions: fill RHS for implicit time-stepping, t=%9.2e\n",t);
    }

  // assignCornerGhostPoints = true if bcOptWave -> cornerWave will assign corner ghost points
    const int assignCornerGhostPoints = bcApproach==useCompatibilityBoundaryConditions ||
                                                                            bcApproach==useLocalCompatibilityBoundaryConditions;

    bool useOpt= true; // twilightZone && true;

    bool isAllImplicit = true;  // all grids are implicit
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
    // BCs for implicit time-stepping are now filled in here *wdh* Dec 13, 2023.
    // OLD: +++ NOTE: BCs for implicit are done in takeImplicitTimeStep +++
        isAllImplicit = isAllImplicit && gridIsImplicit(grid);
    }
    if( applyExplicitBoundaryConditions )  
        isAllImplicit=false; 


  // --- apply known solutions on "exact" boundaries ----
    if( knownSolutionOption=="userDefinedKnownSolution" ) 
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();
            
            realArray & ug = u[grid];
            ForBoundary(side,axis) 
            {
                if( mg.boundaryCondition(side,axis)==exactBC )
                {
          // -- assign boundary and ghost values to be exact ---
                    if( t<=2*dt )
                        printF("applyBC: Set knownSolution at exact boundary (side,axis,grid)=(%d,%d,%d), numGhost=%d\n",side,axis,grid,numGhost);

                    getBoundaryIndex(gid,side,axis,I1,I2,I3,numGhost);
                    if( side==0 )
                        Iv[axis] = Range(gid(0,axis)-numGhost,gid(0,axis));
                    else 
                        Iv[axis] = Range(gid(1,axis),gid(1,axis)+numGhost);

                    getUserDefinedKnownSolution( t, grid, ug, I1,I2,I3 );
                }
            }
        }

    }

    if( debug & 4 && t<=2*dt )
        printF("CgWave: applyBC: useOpt=%d\n",(int)useOpt);

    if( useOpt )
    {
    // ----- Optimized Boundary Conditions ------
        bool isImplicit    = false; // at least one grid is implicit
    // bool isAllImplicit = true;  // all grids are implicit
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            realMappedGridFunction & ug = u[grid];

            isImplicit  = (isImplicit || gridIsImplicit(grid) ) && !applyExplicitBoundaryConditions; 
      // isAllImplicit = isAllImplicit && gridIsImplicit(grid);

      // +++ NOTE: BCs for implicit are done in takeImplicitTimeStep +++

            
            if( !gridIsImplicit(grid) || applyExplicitBoundaryConditions || fillImplicitBoundaryConditions )
            {
        // --------- EXPLICIT GRID ------------

                if( bcApproach==useLocalCompatibilityBoundaryConditions )
                {
                    if( fillImplicitBoundaryConditions )
                    {
                        printF("applyBoundaryConditions:ERROR: bcApproach==useLocalCompatibilityBoundaryConditions not implemented for implicit time-stepping\n");
                        OV_ABORT("error");
                    }

          // -- Apply local compatibility boundary conditions ---
                    ug.applyBoundaryCondition(0,BCTypes::dirichlet,dirichlet,0.,t);   // *** FIX ME **** THIS SHOULD BE u[grid] **********
          // u.applyBoundaryCondition(0,BCTypes::dirichlet,dirichlet,0.,t);   // *** FIX ME **** THIS SHOULD BE u[grid] **********

                    assignLCBC( u[grid], t, dt, grid );

          // OV_ABORT("applyBC: FINISH ME: for bcApproach==useLocalCompatibilityBoundaryConditions");
                    if( useUpwindDissipation )
                    {
                        const int ghost = orderOfAccuracy/2 + 1; // extrapolate an extra ghost
                        bcParams.ghostLineToAssign=ghost;
                        ug.applyBoundaryCondition(0,BCTypes::extrapolate,dirichlet,0.,t,bcParams);  // *** FIX ME **** THIS SHOULD BE u[grid] **********
                        ug.applyBoundaryCondition(0,BCTypes::extrapolate,neumann  ,0.,t,bcParams);  // *** FIX ME **** THIS SHOULD BE u[grid] **********
                        bcParams.ghostLineToAssign=1; // reset  
                    }         

                }
                else
                {     
                    MappedGrid & mg = cg[grid];
          // const IntegerArray & gid = mg.gridIndexRange();
                    
                    OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
                    OV_GET_SERIAL_ARRAY(Real,un[grid],unLocal);
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
                            else
                            {
                                knownSolutionOption=1000;  // all other user defined known solutions
                            }   
                        }
                        int addForcingBC= forcingOption==noForcing ? 0 : 1;
                        if( applyKnownSolutionAtBoundaries )
                            addForcingBC=1; 
                        int gridType = isRectangular ? 0 : 1;
            // int gridIsImplicit = 0; 
                        int numberOfComponents = 1;          // for now CgWave only has a single component 
                        int uc = 0;                          // first component
                        int ipar[] = {
                            uc                  ,            // ipar( 0)
                            numberOfComponents  ,            // ipar( 1)
                            grid                ,            // ipar( 2)
                            gridType            ,            // ipar( 3)
                            orderOfAccuracy     ,            // ipar( 4)
                            gridIsImplicit(grid),            // ipar( 5)  // added Nov 22, 2023
                            twilightZone        ,            // ipar( 6)
                            np                  ,            // ipar( 7)
                            debug               ,            // ipar( 8)
                            myid                ,            // ipar( 9)
                            applyKnownSolutionAtBoundaries,  // ipar(10)
                            knownSolutionOption,             // ipar(11)
                            addForcingBC,                    // ipar(12)
                            forcingOption,                   // ipar(13)
                            useUpwindDissipation,            // ipar(14)
                            numGhost,                        // ipar(15)
                            assignBCForImplicit,             // ipar(16)
                            bcApproach,                      // ipar(17)
                            numberOfFrequencies,             // ipar(18)
                            assignCornerGhostPoints,         // ipar(19)
                            orderOfExtrapolation             // ipar(20)
                                                  };
                        Real cEM2 = c;
                        if( solveHelmholtz ) 
                        {
              // Adjust c for the EM2 absorbing BC to account for time-discretization errors
              //   D+t (Dx ) w + A+( ... )
                            if( frequencyArraySave(0)*dt>0 )
                            {
                                cEM2 = c*tan(frequencyArray(0)*dt/2.)/(frequencyArraySave(0)*dt/2.);
                // printF("\n XXXXXXX cEM2=%e XXXXXXX\n\n",cEM2);
                            }
                            else
                            {
                                if( debug>3 )
                                {
                                    if( frequencyArraySave(0)==0 )
                                        printF("WARNING:bcMacros frequencyArraySave(0)=%12.4e. NOT ADJUSTING c for EM2 absorbing BC\n",frequencyArraySave(0));
                                    if( dt<=0 )
                                        printF("WARNING:bcMacros dt<= 0 ! dt=%12.4e. NOT ADJUSTING c for EM2 absorbing BC\n",dt);
                                }
                            }
                        }           
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
                            c,                //  rpar(10)
                            cEM2              //  rpar(11)
                                                    };
                        real *pu = uLocal.getDataPointer();
                        real *pun = unLocal.getDataPointer();
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
                        real *puTemp1=&temp, *puTemp2=&temp;
                        if( orderOfAccuracy==4 && bcApproach==useCompatibilityBoundaryConditions )
                        {
              // WE currently need work space for bcOptWave and standard CBC at order four ******************FIX ME: just make a stencil ****
                            if( !dbase.has_key("uTemp1") )
                            {
                                RealArray & uTemp1 = dbase.put<RealArray>("uTemp1");
                                RealArray & uTemp2 = dbase.put<RealArray>("uTemp2");
                // -- find the grid with physical boundaries with the most points ---
                                int numElements=0;
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    MappedGrid & mg = cg[grid];
                                    const IntegerArray & boundaryCondition = mg.boundaryCondition();
                                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    if( max(boundaryCondition)>0 )
                                        numElements = max( numElements, maskLocal.elementCount() );
                                }
                                printF(">>> INFO CgWave::applyBC: allocate uTemp1 and utemp2 for order 4 CBC: numElements=%d\n",numElements);
                                if( numElements>0 )
                                {
                                    uTemp1.redim(numElements);
                                    uTemp2.redim(numElements);
                                }
                            }
                            RealArray & uTemp1 = dbase.get<RealArray>("uTemp1");
                            RealArray & uTemp2 = dbase.get<RealArray>("uTemp2");
                            puTemp1 = uTemp1.getDataPointer();
                            puTemp2 = uTemp2.getDataPointer();
                        }

          // if( true )
          //   ::display(bcLocal,"applyBC: bcLocal","%3i ");
                    int ierr=0;
                    bcOptWave(mg.numberOfDimensions(),
                                        uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
                                        uLocal.getBase(2),uLocal.getBound(2),
                                        indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),
                                        *pu, *pun, *pmask, *prsxy, *pxy, *puTemp1, *puTemp2,
                                        bcLocal(0,0), frequencyArray(0),
                                        pdb, ipar[0],rpar[0], ierr );

          // NOTES:
          //   DO NOT CALL ABC FOR EM order 4 implicit -- FIX ME FOR UPWINDING 
          //   DO CALL ABC FOR EM ORDER 2 implicit and upwinding since we need to assign ABCs after upwind step

                    bool applyABCs=  !gridIsImplicit(grid); // no need to apply ABC's if we are implicit
                    if( useUpwindDissipation && !implicitUpwind )
                        applyABCs=true;               // do call ABCs if we are upwinding 

                    if( !fillImplicitBoundaryConditions && applyABCs )
          // if( !fillImplicitBoundaryConditions && (orderOfAccuracy==2 || !gridIsImplicit(grid)) )
          // if( !fillImplicitBoundaryConditions )
                    {
            // Non-reflecting and Absorbing boundary conditions
            // ***NOTE*** symmetry corners and edges are assigned in this next routine *fix me*
            // OV_GET_SERIAL_ARRAY(Real,un[grid],unLocal);

                        if( gridIsImplicit(grid) && orderOfAccuracy==4 )
                        {
                            printF("applyBC:ERROR: implicit time-stepping, order=4, EM RBCs -- cannot use abcWave since order=4 Em BC's use extrapolation instead to CBCs\n");
                            OV_ABORT("FINISH ME");
                        }

                        Real *uptr = uLocal.getDataPointer();
                        Real *uOldptr = unLocal.getDataPointer();

                        RealArray ff(2,2);
                        Real *fptr = ff.getDataPointer();

                        abcWave( mg.numberOfDimensions(), 
                                                uLocal.getBase(0),uLocal.getBound(0),
                                                uLocal.getBase(1),uLocal.getBound(1),
                                                uLocal.getBase(2),uLocal.getBound(2),
                                                ff.getBase(0),ff.getBound(0),
                                                ff.getBase(1),ff.getBound(1),
                                                ff.getBase(2),ff.getBound(2),
                                                indexRangeLocal(0,0),
                                                *uOldptr, *uptr, *fptr, *pmask, *prsxy, *pxy,
                                                bcLocal(0,0), bcLocal(0,0), ipar[0], rpar[0], ierr );
                    }
                }

        // ...swap periodic edges ....
                if( !fillImplicitBoundaryConditions)
                {
                    ug.periodicUpdate();
                    ug.updateGhostBoundaries();
                }

            }
            else
            {
        // -------- IMPLICIT GRID (but we are not filling the RHS for the implicit solve) ------

                if( useUpwindDissipation && !implicitUpwind )
                {
          // IPCU -- upwinding added in an explicit step 
                    if( debug & 4 && t<=2*dt )
                    {
                        printF("+++ extrap extra ghost for IPCU scheme\n");
                    }

                    const int ghost = orderOfAccuracy/2 + 1; // extrapolate an extra ghost
                    bcParams.ghostLineToAssign=ghost;
                    ug.applyBoundaryCondition(0,BCTypes::extrapolate,dirichlet,0.,t,bcParams);
                    ug.applyBoundaryCondition(0,BCTypes::extrapolate,neumann  ,0.,t,bcParams);
                    bcParams.ghostLineToAssign=1; // reset  
                } 

            }



        } // end for grid 


    // For upwind dissipation we should assign another line of points next to interpolation points
    // to support the wider upwind stencil
        if( !fillImplicitBoundaryConditions)
        {
                if( useUpwindDissipation )
                {
                    if( assignInterpNeighbours==interpolateInterpNeighbours ||
                            assignInterpNeighbours==defaultAssignInterpNeighbours  )
                    {
            // --- Interpolate interp neighbours *new* June 10, 2022
            // printF("Call assignInterpolationNeighbours...\n");
                        if( !dbase.has_key("interpNeighbours") )
                        {
                            AssignInterpNeighbours & interpNeighbours = dbase.put<AssignInterpNeighbours>("interpNeighbours");
                            interpNeighbours.setAssignmentType( AssignInterpNeighbours::interpolateInterpolationNeighbours );
              // Interpolation width: This could potentially be 1 less than the normal interp width: 
                            const int interpolationWidth = orderOfAccuracy+1;  
                            interpNeighbours.setInterpolationWidth( interpolationWidth );
                        }   
                        AssignInterpNeighbours & interpNeighbours = dbase.get<AssignInterpNeighbours>("interpNeighbours");
                        const int numberOfComponents=1;  // ** FIX ME **
                        Range C = numberOfComponents;
            // We could pass the TZ function for checking errors 
                        OGFunction* pExact = NULL; // twilightZone ? dbase.get<OGFunction* >("tz") : NULL;
                        if( debug & 8 )
                            printf("++++ CgWave::applyBC: interpolate Interpolation Neighbours at t=%9.3e\n",t);
                        interpNeighbours.assignInterpolationNeighbours( u, C, pExact, t  );
            // printF("Done call assignInterpolationNeighbours.\n");
                    }
                    else
                    {
                        assert( assignInterpNeighbours==extrapolateInterpNeighbours );
                        if( !isAllImplicit )
                        {
                            printF("Call extrapolateInterpolationNeighbours...\n");
                            u.applyBoundaryCondition(0,BCTypes::extrapolateInterpolationNeighbours,BCTypes::allBoundaries,0.,t,bcParams);
                            printF("Done call extrapolateInterpolationNeighbours.\n");
                        }
                        else
                        {
                            if( isImplicit && !isAllImplicit )
                            {
                                printF("applyBoundaryConditions: t=%9.3e: FIX ME: extrapolateInterpolationNeighbours for PARTIALY IMPLICIT time-stepping\n"); 
                                OV_ABORT("error");
                            }
                        }      
                    }
                }
        }

        timing(timeForBoundaryConditions) += getCPU()-cpu0;

    // ====================== RETURN ===========
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

  // For upwind dissipation we should assign another line of points next to interpolation points
  // to support the wider upwind stencil


    if( useUpwindDissipation )
    {
    // For upwind dissipation we assign another line of points next to interpolation points
    // to support the wider upwind stencil    
      // assign another line of points next to interpolation points
            if( assignInterpNeighbours==interpolateInterpNeighbours ||
                    assignInterpNeighbours==defaultAssignInterpNeighbours )
            {
        // --- Interpolate interp neighbours *new* June 10, 2022
                if( !dbase.has_key("interpNeighbours") )
                {
                    AssignInterpNeighbours & interpNeighbours = dbase.put<AssignInterpNeighbours>("interpNeighbours");
                    interpNeighbours.setAssignmentType( AssignInterpNeighbours::interpolateInterpolationNeighbours );
          // Interpolation width: This could potentially be 1 less than the normal interp width: 
                    const int interpolationWidth = orderOfAccuracy+1;  
                    interpNeighbours.setInterpolationWidth( interpolationWidth );
                }   
                AssignInterpNeighbours & interpNeighbours = dbase.get<AssignInterpNeighbours>("interpNeighbours");
                const int numberOfComponents=1;  // ** FIX ME **
                Range C = numberOfComponents;
        // We could pass the TZ function for checking errors 
                OGFunction* pExact = NULL; // twilightZone ? dbase.get<OGFunction* >("tz") : NULL;
        // printf("++++ CgWave::applyBC: interpolate Interpolation Neighbours at t=%9.3\n",t);
                interpNeighbours.assignInterpolationNeighbours( u, C, pExact, t  );
            }
            else
            {
        // -- extrapolate interpolation neighbours
                assert( assignInterpNeighbours==extrapolateInterpNeighbours );
                u.applyBoundaryCondition(0,BCTypes::extrapolateInterpolationNeighbours,BCTypes::allBoundaries,0.,t,bcParams);
            }
    }



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


  // if( bcApproach != useLocalCompatibilityBoundaryConditions )  // *wdh* March 27, 2023 LCBC should do corners I think
    if( !assignCornerGhostPoints ) // *wdh* Nov 29, 2023 bcOptWave -> cornerWave will assign corner ghost in some cases
    {
        if( t<2*dt )
            printF("applyBC: call finishBoundaryConditions to assign points at corners\n");

        BoundaryConditionParameters extrapParams;
    // extrapParams.orderOfExtrapolation=orderOfAccuracy+1; // *wdh* July 21, 2024
        extrapParams.orderOfExtrapolation=orderOfExtrapolation;
    // this function sets corner ghost, does a periodic update and parallel ghost update: 
        u.finishBoundaryConditions(extrapParams);
    }
    else
    {
    // I think this is needed:
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            u[grid].periodicUpdate();
            u[grid].updateGhostBoundaries();
        }
    }

    timing(timeForBoundaryConditions) += getCPU()-cpu0;

  // printF("\n ++++ timing(timeForBoundaryConditions) = %9.3e\n",timing(timeForBoundaryConditions));

        
    return 0;
}


// ======================================================================================================
/// \brief Apply boundary conditions for an eigenfunction
/// 
/// \param u (input/output) : apply BCs to this grid function
// ======================================================================================================
int CgWave::
applyEigenFunctionBoundaryConditions( realCompositeGridFunction & u )
{

    const Real t=0; 

    const int & orderOfAccuracy      = dbase.get<int>("orderOfAccuracy");
    const int & orderOfExtrapolation = dbase.get<int>("orderOfExtrapolation");
  // const Real & c                   = dbase.get<real>("c");
  // const real & dt                  = dbase.get<real>("dt");
    const int & upwind               = dbase.get<int>("upwind");
    const int & implicitUpwind       = dbase.get<int>("implicitUpwind");

    bool useUpwindDissipation        = upwind!=0;

    const AssignInterpolationNeighboursEnum & assignInterpNeighbours = 
                                                          dbase.get<AssignInterpolationNeighboursEnum>("assignInterpNeighbours");


    BoundaryConditionParameters bcParams;
    bcParams.orderOfExtrapolation=orderOfAccuracy+1;


    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];

  // printF("Call interpolate...\n");
    real cpu0 = getCPU();

  // *** Note: interpolate also does a periodic update **** 

    u.interpolate();

    real cpu1 = getCPU();
    timing(timeForInterpolate)+= cpu1-cpu0;
  // printF("Done interpolate.\n");
    
    cpu0=cpu1;

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        const IntegerArray & boundaryCondition = mg.boundaryCondition();
        if( orderOfAccuracy>2 && min(abs(boundaryCondition-neumann))==0 )
        {
            printF("applyEigenFunctionBoundaryConditions:ERROR: finish Neumann BC for order=%d\n",orderOfAccuracy);
            OV_ABORT("finish me");
        }

        realMappedGridFunction & ug = u[grid];

        ug.applyBoundaryCondition(0,BCTypes::dirichlet,dirichlet,0.,t); 
        ug.applyBoundaryCondition(0,BCTypes::neumann,neumann,0.,t);  

    //       assignLCBC( u[grid], t, dt, grid );

    // -- ghost points -- do this for now
        int ghost=1;
        bcParams.ghostLineToAssign=ghost;
        ug.applyBoundaryCondition(0,BCTypes::extrapolate,dirichlet,0.,t,bcParams);  // *** FIX ME **** THIS SHOULD BE u[grid] **********
        if( useUpwindDissipation )
        {
            ghost = orderOfAccuracy/2 + 1; // extrapolate an extra ghost
            bcParams.ghostLineToAssign=ghost;
            ug.applyBoundaryCondition(0,BCTypes::extrapolate,dirichlet,0.,t,bcParams);  // *** FIX ME **** THIS SHOULD BE u[grid] **********
            ug.applyBoundaryCondition(0,BCTypes::extrapolate,neumann  ,0.,t,bcParams);  // *** FIX ME **** THIS SHOULD BE u[grid] **********
            bcParams.ghostLineToAssign=1; // reset  
        }           

    }

  // For upwind dissipation we should assign another line of points next to interpolation points
  // to support the wider upwind stencil
    bool isImplicit=false, isAllImplicit=false;
        if( useUpwindDissipation )
        {
            if( assignInterpNeighbours==interpolateInterpNeighbours ||
                    assignInterpNeighbours==defaultAssignInterpNeighbours  )
            {
        // --- Interpolate interp neighbours *new* June 10, 2022
        // printF("Call assignInterpolationNeighbours...\n");
                if( !dbase.has_key("interpNeighbours") )
                {
                    AssignInterpNeighbours & interpNeighbours = dbase.put<AssignInterpNeighbours>("interpNeighbours");
                    interpNeighbours.setAssignmentType( AssignInterpNeighbours::interpolateInterpolationNeighbours );
          // Interpolation width: This could potentially be 1 less than the normal interp width: 
                    const int interpolationWidth = orderOfAccuracy+1;  
                    interpNeighbours.setInterpolationWidth( interpolationWidth );
                }   
                AssignInterpNeighbours & interpNeighbours = dbase.get<AssignInterpNeighbours>("interpNeighbours");
                const int numberOfComponents=1;  // ** FIX ME **
                Range C = numberOfComponents;
        // We could pass the TZ function for checking errors 
                OGFunction* pExact = NULL; // twilightZone ? dbase.get<OGFunction* >("tz") : NULL;
                if( debug & 8 )
                    printf("++++ CgWave::applyBC: interpolate Interpolation Neighbours at t=%9.3e\n",t);
                interpNeighbours.assignInterpolationNeighbours( u, C, pExact, t  );
        // printF("Done call assignInterpolationNeighbours.\n");
            }
            else
            {
                assert( assignInterpNeighbours==extrapolateInterpNeighbours );
                if( !isAllImplicit )
                {
                    printF("Call extrapolateInterpolationNeighbours...\n");
                    u.applyBoundaryCondition(0,BCTypes::extrapolateInterpolationNeighbours,BCTypes::allBoundaries,0.,t,bcParams);
                    printF("Done call extrapolateInterpolationNeighbours.\n");
                }
                else
                {
                    if( isImplicit && !isAllImplicit )
                    {
                        printF("applyBoundaryConditions: t=%9.3e: FIX ME: extrapolateInterpolationNeighbours for PARTIALY IMPLICIT time-stepping\n"); 
                        OV_ABORT("error");
                    }
                }      
            }
        }

  // WARNING -- do NOT do this for LCBC !! **********************
    BoundaryConditionParameters extrapParams;
  // extrapParams.orderOfExtrapolation=orderOfAccuracy+1;  // *wdh* Jul 21, 2024
    extrapParams.orderOfExtrapolation=orderOfExtrapolation; 
    u.finishBoundaryConditions(extrapParams);

    timing(timeForBoundaryConditions) += getCPU()-cpu0;

    return 0;
}




