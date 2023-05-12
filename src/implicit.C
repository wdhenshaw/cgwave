// This file automatically generated from implicit.bC with bpp.
#include "CgWave.h"
#include "display.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "Oges.h"
#include "CompositeGridOperators.h"
#include "SparseRep.h"

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


// The getBcOptParameters macro is defined here:
// =======================================================================================
// Macro: get the parameters for calling the optimized fortran BC routine
// =======================================================================================



#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )


// ============================================================================================
/// \brief Take an implicit time step to new time t
/// \param t (input) : new time
// ============================================================================================
int CgWave::takeImplicitStep( Real t )
{
    real cpu0=getCPU();

    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::numberOfProcessors());

    if( !dbase.has_key("impSolver") )
    {
        formImplicitTimeSteppingMatrix();
    }

    const int & debug                = dbase.get<int>("debug");
    FILE *& debugFile                = dbase.get<FILE*>("debugFile");
    FILE *& pDebugFile               = dbase.get<FILE*>("pDebugFile");

    const Real & c                   = dbase.get<real>("c");
    const real & dt                  = dbase.get<real>("dt");
    const int & orderOfAccuracy      = dbase.get<int>("orderOfAccuracy");

    const int & numberOfFrequencies  = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray = dbase.get<RealArray>("frequencyArray");
    const int & solveHelmholtz       = dbase.get<int>("solveHelmholtz");
    const int & computeEigenmodes    = dbase.get<int>("computeEigenmodes");

    const int & upwind               = dbase.get<int>("upwind");
    const int & implicitUpwind       = dbase.get<int>("implicitUpwind");
    bool useUpwindDissipation        = upwind;
    int & totalImplicitIterations    = dbase.get<int>("totalImplicitIterations");  
    int & totalImplicitSolves        = dbase.get<int>("totalImplicitSolves");  

    const int & addForcing                  = dbase.get<int>("addForcing");
    const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const bool twilightZone = forcingOption==twilightZoneForcing; 

    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");
    const int & applyKnownSolutionAtBoundaries = dbase.get<int>("applyKnownSolutionAtBoundaries"); // by default, do NOT apply known solution at boundaries

    const int & current                  = dbase.get<int>("current"); // hold the current solution index
    const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");  

    const int cur = current;   // current time level
    const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
    const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

    realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
    realCompositeGridFunction & up = u[prev];       // previous time 
    realCompositeGridFunction & uc = u[cur];        // current time 
    realCompositeGridFunction & un = u[next];       // new time

    Oges & impSolver = dbase.get<Oges>("impSolver");
    if( impSolver.isSolverIterative() )
    {
    // Iterative solvers need a separate RHS 
        if( !dbase.has_key("impSolverRHS") )
        {
            printF(">>>>cgWave::INFO: impSolver is iterative, we need a separate RHS.\n");
            dbase.put<realCompositeGridFunction>("impSolverRHS");
            realCompositeGridFunction & rhs = dbase.get<realCompositeGridFunction>("impSolverRHS");
            rhs.updateToMatchGrid(cg);
        }
        realCompositeGridFunction & rhs = dbase.get<realCompositeGridFunction>("impSolverRHS");
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,up[grid],upLocal);
            OV_GET_SERIAL_ARRAY(Real,uc[grid],ucLocal);
            OV_GET_SERIAL_ARRAY(Real,un[grid],unLocal);
            OV_GET_SERIAL_ARRAY(Real,rhs[grid],rhsLocal);

      // --- For implicit time-stepping, unLocal holds the RHS to the implicit equations: 
            rhsLocal=unLocal;

      // ---- initial guess for iterative solver ----
      //  *to do* For periodic solutions use periodic guess

            if( solveHelmholtz && numberOfFrequencies==1 && !computeEigenmodes )
            {
        // Assume: 
        //    uc = v(x,t)    = u(x)*cos(omega*t)
        //    un = v(x,t-dt) = u(x)*cos(omega*(t-dt))
        // Thus
        //    v(x,t) = v(x,t-dt)* cos(omega*t)/cos(omega*(t-dt))
                Index I1,I2,I3;
                getIndex(cg[grid].dimension(),I1,I2,I3);
                bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3);        

                if( ok )
                {
                    if( grid==0 && ( debug & 4 ) )
                        printF("CgWave:impSolve: choose initial guess for implicit solver assuming a periodic solution.\n");
           // Real diff = max(fabs(ucLocal(I1,I2,I3)-vLocal(I1,I2,I3,0)));
           // printF(" >> diff | uc - v |=%9.2e\n",diff);

                      Real cosFreqt  = cos( frequencyArray(0)*(t-dt) );
                      if( fabs(cosFreqt)<REAL_EPSILON*100. ) cosFreqt=1.;
                      Real cosFactor = cos( frequencyArray(0)*(t) )/cosFreqt;

                      unLocal(I1,I2,I3) = ucLocal(I1,I2,I3)*cosFactor;

           // ****** WRONG: un and uc ONLY HAVE ONE COMPONENT ******************* FIX ME 
                      for( int freq=1; freq<numberOfFrequencies; freq++ )
                      {
                          cosFreqt  = cos( frequencyArray(freq)*(t-dt) );
                          if( fabs(cosFreqt)<REAL_EPSILON*100. ) cosFreqt=1.;
                          cosFactor = cos( frequencyArray(freq)*(t) )/cosFreqt;

                          unLocal(I1,I2,I3) += ucLocal(I1,I2,I3,freq)*cos( frequencyArray(freq)*dt );
                      }
                }
            }
            else
            {
        // unLocal = ucLocal;           // initial guess *** FIX ME ***
                unLocal = 2.*ucLocal - upLocal; // initial guess *** FIX ME ***
            }


        }
    }
    realCompositeGridFunction & rhs = impSolver.isSolverIterative() ? dbase.get<realCompositeGridFunction>("impSolverRHS") : un;


  // *wdh* May 2, 2023
  // CompositeGridOperators & op = dbase.get<CompositeGridOperators>("operators");
  // un.setOperators( op );
  // rhs.setOperators( op );

  // --- Fill in the RHS for implicit boundary conditions ----

    int numGhost = orderOfAccuracy/2;
    if( useUpwindDissipation && implicitUpwind ) numGhost++; 

    if( false )
    {
        rhs.display("RHS before filling BCs");
    }

    const int assignBCForImplicit = 1;  // used in call to bcOptWave

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
    // const IntegerArray & gid = mg.gridIndexRange();
        
        OV_GET_SERIAL_ARRAY(Real,rhs[grid],rhsLocal);
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

    // get parameters for calling fortran
            IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
            ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( rhs[grid],indexRangeLocal,dimLocal,bcLocal );
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
            real *pu = rhsLocal.getDataPointer();
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
    // if( true )
    //   ::display(bcLocal,"implicit: bcLocal","%3i ");

        bcOptWave(mg.numberOfDimensions(),
                            rhsLocal.getBase(0),rhsLocal.getBound(0),rhsLocal.getBase(1),rhsLocal.getBound(1),
                            rhsLocal.getBase(2),rhsLocal.getBound(2),
                            indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),
                            *rhsLocal.getDataPointer(), *pmask, *prsxy, *pxy,  bcLocal(0,0), frequencyArray(0),
                            pdb, ipar[0],rpar[0], ierr );

    } // end for grid 



  // ------- SOLVE THE IMPLICIT EQUATIONS -----

    bool outputMatrix=false; // true;  // for debugging 

    if( outputMatrix )
    {
        Oges::debug=63;
        impSolver.set(OgesParameters::THEkeepSparseMatrix,true);
    }

  // int solverType;
  // impSolver.get(OgesParameters::THEsolverType,solverType);
    if( impSolver.isSolverIterative() )
    {
        if( false )
        {
            rhs.display("RHS before implicit solve");
            OV_ABORT("stop here for now");
        }

    // **  for iterative solvers we want un to be the good guess for u(t+dt) **************

        impSolver.solve( un,rhs );  

        int numIterations = impSolver.getNumberOfIterations();
        totalImplicitIterations += numIterations;
        totalImplicitSolves++;
        if( true )
            printF("CgWave::takeImplicitStep: implicit solve: max-res= %8.2e (iterations=%i) ***\n",
                                  impSolver.getMaximumResidual(),numIterations);
    }
    else
    {
        impSolver.solve( un,un );   
    }
    if( outputMatrix )
    {
        printF("cgWaveINFO: save the implicit matrix to file cgWaveMatrix.out (using writeMatrixToFile). \n");
        impSolver.writeMatrixToFile("cgWaveMatrix.out");

        aString fileName = "cgWaveSparseMatrix.dat";
        printF("cgWave:INFO: save the implicit matrix to file %s (using outputSparseMatrix)\n",(const char*)fileName);
        impSolver.outputSparseMatrix( fileName );

        OV_ABORT("stop here for now"); 
    }

    timing(timeForImplicitSolve) += getCPU()-cpu0;

    return 0;
}





#define ForStencil(m1,m2,m3)   for( m3=-halfWidth3; m3<=halfWidth3; m3++) for( m2=-halfWidth2; m2<=halfWidth2; m2++) for( m1=-halfWidth1; m1<=halfWidth1; m1++)

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)

// =======================================================================
// indexToEquation( n,i1,i2,i3 ) : defines the global index for each unknown in the system
//     n=component number (uc,vc,...)
//    (i1,i2,i3) = grid point
// =======================================================================
#define indexToEquation( n,i1,i2,i3 ) (n+1+ numberOfComponentsForCoefficients*(i1-equationNumberBase1+equationNumberLength1*(i2-equationNumberBase2+equationNumberLength2*(i3-equationNumberBase3))) + equationOffset)

// =======================================================================
// =======================================================================
#define setEquationNumber(m, ni,i1,i2,i3,  nj,j1,j2,j3 )equationNumber(m,i1,i2,i3)=indexToEquation( nj,j1,j2,j3)

// =======================================================================
// =======================================================================
#define setClassify(n,i1,i2,i3, type) classify(i1,i2,i3,n)=type

// ==========================================================================
// Macro: setup the variables needed to fill a sparse matrix on a mappedGrid
// ==========================================================================

// ==========================================================================
// Macro: fill the matrix with extrapolation for a given ghost=1,2,3,...
// ==========================================================================



// ============================================================================================
/// \brief Form and the matrix for implicit time-stepping
// ============================================================================================
int CgWave::formImplicitTimeSteppingMatrix()
{
    real cpu0=getCPU();

    int & debug                          = dbase.get<int>("debug");
    const Real & c                       = dbase.get<real>("c");
    const real & dt                      = dbase.get<real>("dt");
    const int & orderOfAccuracy          = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime    = dbase.get<int>("orderOfAccuracyInTime");
    const IntegerArray & gridIsImplicit  = dbase.get<IntegerArray>("gridIsImplicit");

    const int & upwind                   = dbase.get<int>("upwind");
    const int & implicitUpwind           = dbase.get<int>("implicitUpwind");

    const bool addUpwinding = upwind && implicitUpwind;

    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");

  // addUpwinding=false; // *********** TURN OFF FOR NOW ***************


    printF("\n ==================== FORM MATRIX FOR IMPLICIT TIME-STEPPING ===================\n");
    printF("   c=%.4g, dt=%9.3e, orderOfAccuracy=%d, orderOfAccuracyInTime=%d addUpwinding=%d, bcApproach=%d\n",
              c,dt,orderOfAccuracy,orderOfAccuracyInTime,addUpwinding,(int)bcApproach);
    printF(" ================================================================================\n");

  

    const int & numberOfComponentGrids = cg.numberOfComponentGrids(); 
    const int & numberOfDimensions = cg.numberOfDimensions(); 

  // coefficients in implicit time-stepping  
  //  D+t D-t u = c^2 Delta( cImp(1,0) *u^{n+1} + cImp(0,0) *u^n + cImp(-1,0)* u^{n-1} )
    RealArray & cImp              = dbase.get<RealArray>("cImp");  

    if( !dbase.has_key("impSolver") )
    {
        Oges & impSolver =  dbase.put<Oges>("impSolver");
        impSolver.setGridName( dbase.get<aString>("nameOfGridFile") );
    }
    Oges & impSolver = dbase.get<Oges>("impSolver");

  // impSolver.updateToMatchGrid( cg );                     

    if( dbase.has_key("implicitSolverParameters") )
    {  
        printf("CgWave::formImplicitTimeSteppingMatrix: Changing the implicit solver parameters. \n");
        OgesParameters & par = dbase.get<OgesParameters>("implicitSolverParameters");
        impSolver.setOgesParameters(par);
    }
    else
    {
        int solverType=OgesParameters::yale; 

    // solverType=OgesParameters::PETSc;
    // solverType=OgesParameters::PETScNew; // parallel

        impSolver.set(OgesParameters::THEsolverType,solverType); 

        if( solverType==OgesParameters::PETSc )
          impSolver.set(OgesParameters::THEsolverMethod,OgesParameters::biConjugateGradientStabilized);

    // impSolver.set(OgesParameters::THEparallelSolverMethod,OgesParameters::gmres);
    // impSolver.set(OgesParameters::THErelativeTolerance,max(tol,REAL_EPSILON*10.));
    // impSolver.set(OgesParameters::THEmaximumNumberOfIterations,10000);
    // if( iluLevels>=0 )
    //   impSolver.set(OgesParameters::THEnumberOfIncompleteLULevels,iluLevels);
    }

    int solverType;
    impSolver.get(OgesParameters::THEsolverType,solverType);

  // Oges::debug=7;
  // printF("Calling impSolver.setGrid()\n");
    impSolver.setGrid( cg );       // this will generate MG levels if using multigrid
  // printF("Calling impSolver.updateToMatchGrid()\n");
  // impSolver.updateToMatchGrid( cg ); 

  // --- save the name of the sparse solver we use ---
    if( !dbase.has_key("implicitSolverName") )
        dbase.put<aString>("implicitSolverName");
    aString & implicitSolverName =  dbase.get<aString>("implicitSolverName");
    implicitSolverName = impSolver.parameters.getSolverName();
    printF("\n === Implicit Time-Stepping Solver:\n %s\n =====\n",(const char*)implicitSolverName);

    CompositeGridOperators & op = dbase.get<CompositeGridOperators>("operators");
    op.setOrderOfAccuracy(orderOfAccuracy);

  // -- Use predefined equations in Oges for MG  *wdh* Jan 8 2023 -- TRY THIS ---
    bool useMultigrid= solverType==OgesParameters::multigrid;

  // bool usePredefined = useMultigrid && !addUpwinding && orderOfAccuracy==2; // true = old way
  // We can use predefined equations for
  //  Oges: 
  //    - orderInTime = 2 
  //    - no implicit upwinding
  //    - no CBC's ??
  // Ogmg:
  //    - order in time 2
  //    - no implicit upwinding
  //    - CBC's 
    bool usePredefined = useMultigrid && !addUpwinding && orderOfAccuracyInTime==2; // ** FIX ME for no-MG

    if( orderOfAccuracyInTime != 2  )
    {
        printF("CgWave::formImplicitTimeSteppingMatrix: Error orderOfAccuracyInTime != 2 not implemented yet.\n");
        OV_ABORT("ERROR");    
    }

  // if( !usePredefined &&  bcApproach==useCompatibilityBoundaryConditions )
  // {
  //   printF("CgWave::formImplicitTimeSteppingMatrix: Error useCompatibilityBoundaryConditions not implemented yet\n");
  //   OV_ABORT("ERROR");
  // }

    if( bcApproach==useLocalCompatibilityBoundaryConditions )
    {
        printF("CgWave::formImplicitTimeSteppingMatrix: Error useLocalCompatibilityBoundaryConditions not implemented yet\n");
        OV_ABORT("ERROR");
    }  
  


    if( usePredefined )
    {
    // ***** USE PREDEFINED EQUATIONS IN SOME CASES ******
        printF("\n >>>>>>> CgWave::implicitSolver: use predfined equations <<<<<<< \n\n");

    // ---- use Oges predefined equations ***OLD WAY*** ----
    
        IntegerArray boundaryConditions(2,3,numberOfComponentGrids);
        boundaryConditions = 0;

        RealArray bcData(2,2,3,numberOfComponentGrids);
        bcData=0.;

        Range all; 

    // Solve constCoeff(0,grid)*I +constCoeff(1,grid)*Laplacian 
    // We solve:  I - alpha*(c^2*dt^2)* Delta = ...
        RealArray constantCoeff(2,numberOfComponentGrids);


    // --- SET COEFF IN IMPLICIT MATRIX ----
        printF("CgWave::implicitSolver: A = I -alpha (c*dt)^2 Delta, c=%g, dt=%9.2e\n",c,dt);

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            Real alpha=cImp(-1,0);
            if( gridIsImplicit(grid)==0 )
                alpha=0.; // this grid is explicit

            constantCoeff(0,grid) = 1.;
            constantCoeff(1,grid) = - alpha*SQR(c*dt);
        }

    // Assign boundary conditions for Oges
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            ForBoundary(side,axis)
            {
                  if( mg.boundaryCondition(side,axis)==CgWave::dirichlet )
                  {
                      boundaryConditions(side,axis,grid) = OgesParameters::dirichlet;
                  }
                  else if( mg.boundaryCondition(side,axis)==CgWave::neumann )
                  { 
                      boundaryConditions(side,axis,grid) = OgesParameters::neumann;
                  }
                  else if( mg.boundaryCondition(side,axis) <= 0 )
                  { 
                      boundaryConditions(side,axis,grid) = mg.boundaryCondition(side,axis);
                  }       
                  else if( mg.boundaryCondition(side,axis) > 0 )
                  {
                      printF("CgWave::formImplicitTimeSteppingMatrix:ERROR: unknown boundaryCondition=%d for (side,axis,grid)=(%d,%d,%d)\n", mg.boundaryCondition(side,axis),side,axis,grid);
                      OV_ABORT("ERROR");
                  }

            }

        }


        impSolver.setEquationAndBoundaryConditions( OgesParameters::heatEquationOperator ,op,boundaryConditions,bcData,constantCoeff );

    }
    else
    {
    // --------------------------------------------------------------------------------------
    // ---- Fill in the implicit matrix : more general equations allowed than predefined ----
    // --------------------------------------------------------------------------------------

    // Here are the coefficients in the upwind dissipation operator (D+D-)^p 
    // There are 4 cases depending on whether the full wider stencil is available

    // fourth-order dissipation for 2nd-order scheme:
        Real upwindCoeff4[4][5] = { 1.,-4.,6.,-4.,1.,
                                                                1.,-3.,3.,-1.,0.,   // extrap right-most point D-^3 u(2)
                                                                0.,-1.,3.,-3.,1.,   // extrap left -most point D+^3 u(-2)
                                                                0.,-1.,2.,-1.,0.
                                                            };

    // sixth-order dissipation for 4th-order scheme
        Real upwindCoeff6[4][7] = {1.,-6.,15.,-20.,15.,-6.,1.,
                                                              1.,-5.,10.,-10., 5.,-1.,0.,  // extrap right-most point D-^5 u(3)
                                                              0.,-1., 5.,-10.,10.,-5.,1.,  // extrap left -most point D+^5 u(-3)
                                                              0.,-1., 4., -6., 4.,-1.,0.
              
                                                            };
    //  --- Coefficients in the sosup dissipation from Jordan Angel ---
    // These must match the values in advWave.bf90                          
        Real *upwindCoeff[4];

    // const int upwindHalfStencilWidth = orderOfAccuracy; 
        const int upwindHalfStencilWidth = (orderOfAccuracy+2)/2; 

    // *new* way to define coefficients: (see cgWave.pdf)
        Real adSosup = c*dt/( sqrt(1.*numberOfDimensions) * pow(2.,(orderOfAccuracy+1)) );
        if( (orderOfAccuracy/2) % 2 == 1 )
        {
            adSosup = - adSosup;  // make negative for orders 2,6,10, ...
        }

        if( orderOfAccuracy==2 )
        {
      // adSosup=-c*dt*1./8.;
            for( int m=0; m<4; m++ )
                upwindCoeff[m] =upwindCoeff4[m];
        }
        else if( orderOfAccuracy==4 )
        {
      // adSosup=c*dt*5./288.;
            for( int m=0; m<4; m++ )
                upwindCoeff[m] =upwindCoeff6[m];      
        }
        else if( orderOfAccuracy==6 )
        {
      // adSosup=-c*dt*31./8640.;
      // upwindDissCoeff=upwindDissCoeff8;
            OV_ABORT("ERROR orderOfAccuracy");
        }
        else
        {
            OV_ABORT("ERROR orderOfAccuracy");
        }


        const int cc=0; // component number

        Range all;
        int stencilWidth = orderOfAccuracy + 1;
        int numberOfGhostLines= orderOfAccuracy/2;  
    
    // if( TRUE )
    //     numberOfGhostLines++; // *********************************************** TEST ***********************

        int extraEntries = 1;  // we add 1 extra entry for interpolation equations

        if( addUpwinding )
        {
      // -- Note: we do not always have to add extra entries for upwinding
      //          e.g. on Cartesian grids there are zeros in the existing stencil that could be used. 
            extraEntries = 2*numberOfDimensions; // for upwinding equations
      // stencilWidth += 2;
            numberOfGhostLines += 1;   // for extrapolate an extra ghost line when upwinding
        }

        const int baseStencilSize = pow(stencilWidth,cg.numberOfDimensions());   // number of entries in default stencil 
        const int stencilSize=int( baseStencilSize + extraEntries );             // add extra for interpolation and upwind equations

        const int numberOfComponentsForCoefficients=1;
        const int stencilDimension=stencilSize*SQR(numberOfComponentsForCoefficients);

        const int baseStencilDimension=baseStencilSize*SQR(numberOfComponentsForCoefficients);



        printF(">>>> stenciWidth=%d, stencilSize=%d, numberOfGhostLines=%d, addUpwinding=%d\n",stencilWidth,stencilSize,numberOfGhostLines,addUpwinding);

    // use this coeff matrix: **FIX ME**
        if( !dbase.has_key("impCoeff") )
        {
            dbase.put<realCompositeGridFunction>("impCoeff");
        }
        realCompositeGridFunction & impCoeff = dbase.get<realCompositeGridFunction>("impCoeff"); 

        impCoeff.updateToMatchGrid(cg,stencilDimension,all,all,all); 
    // impCoeff.setIsACoefficientMatrix(true,baseStencilSize,numberOfGhostLines,numberOfComponentsForCoefficients);
        impCoeff.setIsACoefficientMatrix(true,stencilSize,numberOfGhostLines,numberOfComponentsForCoefficients);

    // TROUBLE FOR UPWIND CASE -- Operators probably base matrix side on orderOfAccuracy !!  *** FIX ME ***
        op.setStencilSize(stencilSize);
        op.setNumberOfComponentsForCoefficients(numberOfComponentsForCoefficients);
        impCoeff.setOperators(op);

    // Use these for indexing into coefficient matrices representing systems of equations
    // #define CE(c,e) (baseStencilSize*((c)+numberOfComponentsForCoefficients*(e)))
        #define M123(m1,m2,m3) (m1+halfWidth1+width*(m2+halfWidth2+width*(m3+halfWidth3)))
    // #define M123CE(m1,m2,m3,c,e) (M123(m1,m2,m3)+CE(c,e))    

        Index I1,I2,I3;
        Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
        Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];
        int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];
        int m1,m2,m3; 

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid &mg = cg[grid];
            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

            realMappedGridFunction & coeff = impCoeff[grid];
            MappedGridOperators & mgop = op[grid];

            OV_GET_SERIAL_ARRAY(real,coeff,coeffLocal);
            coeffLocal=0.; 

      // set up some variables we need to index into sparse coefficient matrices
                assert( coeff.sparse!=NULL );
                SparseRepForMGF & sparse = *coeff.sparse;
                int numberOfComponentsForCoefficients = sparse.numberOfComponents;  // size of the system of equations
                int numberOfGhostLines = sparse.numberOfGhostLines;
                int stencilSize = sparse.stencilSize;
                int stencilDim=stencilSize*numberOfComponentsForCoefficients; // number of coefficients per equation
                const int equationOffset=sparse.equationOffset;
                IntegerArray & equationNumber = sparse.equationNumber;
                IntegerArray & classify = sparse.classify;
                const int equationNumberBase1  =equationNumber.getBase(1);
                const int equationNumberLength1=equationNumber.getLength(1);
                const int equationNumberBase2  =equationNumber.getBase(2);
                const int equationNumberLength2=equationNumber.getLength(2);
                const int equationNumberBase3  =equationNumber.getBase(3);
                const int orderOfAccuracy=mgop.getOrderOfAccuracy(); 
        // stencil width's and half-width's :
                const int width = orderOfAccuracy+1;
        // const int width      = stencilWidth;
                const int halfWidth1 = (width-1)/2;
                const int halfWidth2 = numberOfDimensions>1 ? halfWidth1 : 0;
                const int halfWidth3 = numberOfDimensions>2 ? halfWidth1 : 0;
                Range M0 = baseStencilSize;    // ***** May 15, 2021 -> is this right
                Range M = coeff.dimension(0);


      // --- FILL INTERIOR EQUATIONS ----
      // Solve constCoeff(0,grid)*I +constCoeff(1,grid)*Laplacian 
      // We solve:  I - cImp(-1,0) * (c^2*dt^2)* Delta = ...

            const int mDiag = M123(0,0,0);              // index of diagonal entry

            if( gridIsImplicit(grid)==1 )
            {
        // ----- this grid is adavnced with IMPLICIT time-stepping ----
                getIndex(mg.gridIndexRange(),I1,I2,I3);
                RealArray lapCoeff(M0,I1,I2,I3);
                mgop.assignCoefficients(MappedGridOperators::laplacianOperator,lapCoeff,I1,I2,I3,0,0); // 
                

                Real ccLap = - cImp(-1,0)*SQR(c*dt);         // note minus
                coeffLocal(M0,I1,I2,I3)  = ccLap*lapCoeff;

        // set diagonal entry
                coeffLocal(mDiag,I1,I2,I3) += 1.0;

            }
            else
            {
        // ----- this grid is advanced with EXPLICIT time-stepping ----
        // set the matrix the IDENTITY
                printF("+++++ IMPLICIT: grid=%d (%s) IS TREATED EXPLICITLY\n",grid,(const char*)mg.getName());        

        // set diagonal entry
                coeffLocal(mDiag,I1,I2,I3) = 1.0;


            }


            if( addUpwinding )
            {
        // -------------------------------
        // --- ADD UPWIND DISSIPATION ----
        // -------------------------------

                const bool isRectangular = mg.isRectangular();
                Real dx[3]={1.,1.,1.};
                Real dr[3]={1.,1.,1.};
                if( isRectangular )
                    { // rectangular grid grid-spacings: 
                    mg.getDeltaX(dx);
                }
                else
                {
                    mg.update(MappedGrid::THEinverseVertexDerivative );
          // unit square grid spacings: 
                    for( int dir=0; dir<3; dir++ )
                        dr[dir]=mg.gridSpacing(dir);           
                }

                OV_GET_SERIAL_ARRAY_CONDITIONAL(real,mg.inverseVertexDerivative(),rxLocal,!isRectangular);
        // macro to make the rxLocal array look 5-dimensional 
                #define RX(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2))     

                const Real uDotFactor=.5; // from D0t 
                Real adxSosup[3];
        // upwind diss coeff for Cartesian grids: 
                adxSosup[0] = uDotFactor*adSosup/dx[0];
                adxSosup[1] = uDotFactor*adSosup/dx[1];
                adxSosup[2] = uDotFactor*adSosup/dx[2]; 

        // assert( isRectangular );

                FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the domain
                {
                    if( maskLocal(i1,i2,i3)>0 )
                    {
            // --- fill in the coefficients of the upwind dissipation formula ---
            // NOTE: this stencil is wider than the usual
            // NOTE: This formula must agree with the RHS computation in advWave.bf90

                        bool testUPW=false; // true;
                        if( testUPW )
                            coeffLocal(M,i1,i2,i3)=0.; 

                        int extraStencilLocation = baseStencilDimension; // start adding any extra equations here
                        for( int dir=0; dir<numberOfDimensions; dir++ )  // add dissipation along axis "dir"
                        {
                            int idv[3]={0,0,0};
                            idv[dir]=1; // active direction
              // check if left-most and right-most entries in the upwind stencil are valid 
                            const int i1l = i1-upwindHalfStencilWidth*idv[0], i1r = i1+upwindHalfStencilWidth*idv[0];
                            const int i2l = i2-upwindHalfStencilWidth*idv[1], i2r = i2+upwindHalfStencilWidth*idv[1];
                            const int i3l = i3-upwindHalfStencilWidth*idv[2], i3r = i3+upwindHalfStencilWidth*idv[2];

              // Note: there are at most four cases at any order, since we have order/2 layers of interpolation points 
              //  Example, order=2, upwind-order=4
              //    X---X---C---X---X           C = center point = valid discretization point
              //    X---X---C---X               missing right-most 
              //        X---C---X---X           missing left-most
              //        X---C---X               missing left and right-most
                            int upwCase=0; 
                            if( maskLocal(i1l,i2l,i3l)!=0 && maskLocal(i1r,i2r,i3r)!=0 )
                                upwCase=0; // centred, full-width stencil
                            else if( maskLocal(i1l,i2l,i3l)!=0 )
                                upwCase=1; // left biased stencil
                            else if( maskLocal(i1r,i2r,i3r)!=0 )
                                upwCase=2; // right biased stencil   
                            else  
                                upwCase=3; // centred smaller stencil     

                            if( !isRectangular )
                            {
                 // ---Upwind coefficients for a curvilinear grid ---

                 // diss-coeff ~= 1/(change in x along direction r(dir) )
                 // Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                                  if( numberOfDimensions==2 )
                                    adxSosup[dir] = adSosup*uDotFactor*sqrt( SQR(RX(i1,i2,i3,dir,0)) + SQR(RX(i1,i2,i3,dir,1)) )/dr[dir]; 
                                  else
                                    adxSosup[dir] = adSosup*uDotFactor*sqrt( SQR(RX(i1,i2,i3,dir,0)) + SQR(RX(i1,i2,i3,dir,1))  + SQR(RX(i1,i2,i3,dir,2)) )/dr[dir];                    
                            }

              // --- Put the upwind coefficient in the matrix --- 
                            for( int iStencil=-upwindHalfStencilWidth; iStencil<=upwindHalfStencilWidth; iStencil++ )
                            {
                                const int i1s = i1 + iStencil*idv[0], i2s = i2 + iStencil*idv[1],  i3s = i3 + iStencil*idv[2]; // (i1s,i2s,i3s) : stencil index 
                                Real upwStencilValue = upwindCoeff[upwCase][iStencil+upwindHalfStencilWidth]; 
                                if( upwStencilValue != 0. )
                                {
                                    int m; // put into coeff at this index: coeff(m,....) = ...
                                    if( iStencil>= -halfWidth1 && iStencil<=halfWidth1  )
                                    { 
                    // --- point fits in the existing stencil ---

                    //     2  X---X---U---X---X
                    //        |   |   |   |   |
                    //     1  X---X---U---X---X
                    //        |   |   |   |   |
                    //  m2 0  U---U---U---U---U    U = UPWIND Dissipation stencil
                    //        |   |   |   |   |
                    //    -1  X---X---U---X---X
                    //        |   |   |   |   |
                    //    -2  X---X---U---X---X
                    //        -2  -1  0   1   2  
                    //                m1 
                    // 
                    // -- first compute (m1,m2,m3) from iStencil : 
                                        int m1 = iStencil*idv[0], m2=iStencil*idv[1], m3=iStencil*idv[2];
                                        m = M123(m1,m2,m3); 

                                    }
                                    else
                                    { // point is outside the exisiting stencil, use an extra entry in the coefficient matrix 
                                        if( extraStencilLocation >= stencilDimension )
                                        {
                                            printF("[i1,i2]=[%d,%d], dir=%d, iStencil=%d, baseStencilDimension=%d, stencilDimension=%d, extraStencilLocation=%d, halfWidth1=%d upwindHalfStencilWidth=%d\n",
                                                          i1,i2,dir, iStencil,baseStencilDimension, stencilDimension, extraStencilLocation, halfWidth1, upwindHalfStencilWidth);
                                        }
                                        assert( extraStencilLocation<stencilDimension );

                                        m = extraStencilLocation;
                                        extraStencilLocation++;

                                    }

                  // int ii = indexToEquation( cc,i1,i2,i3 ); 
                  // printF("UPWIND: (i1,i2)=(%3d,%3d) dir=%d ->  ADD entry m=%3d, (i1s,i2s)=(%3d,%3d), indexToEqn=%d, adxSosup=%9.3e, upwStencilValue = %9.3e\n",
                  //           i1,i2,dir,m,i1s,i2s,ii,adxSosup[dir],upwStencilValue);

                  // ***TEST*** Just add difference coefficients 
                  // if( true ) upwStencilValue=0.;
                  // Note: adxSosup is NEGATIVE : we subtract upwind coefficient since it has been moved to the LHS

                                    coeffLocal(m,i1,i2,i3) -= adxSosup[dir] * upwStencilValue;
                  // TEST : coeffLocal(m,i1,i2,i3) += adxSosup[dir] * upwStencilValue;


                                    setEquationNumber(m, e,i1,i2,i3,  cc,i1s,i2s,i3s );      // macro to set equationNumber    

                  // equationNumber(m,i1,i2,i3)-=1;  // TEST 
                                }
                            }
                        }
                    } // end if mask 
                }

            }



      // --- FILL BOUNDARY CONDITIONS ----

            const int extrapOrder = orderOfAccuracy+1;

            const Real extrapCoeff3[] = {1.,-3.,3.,-1.};
            const Real extrapCoeff4[] = {1.,-4.,6.,-4.,1.};
            const Real extrapCoeff5[] = {1.,-5.,10.,-10.,5.,-1.};
            const Real *extrapCoeff;
            if( extrapOrder==3 )
                extrapCoeff = extrapCoeff3;
            else if( extrapOrder==4 )
                extrapCoeff = extrapCoeff4;
            else if( extrapOrder==5 )
                extrapCoeff = extrapCoeff5;
            else
              {
                printF("CgWave::formImplicitTimeSteppingMatrix unexpected extrapOrder=%d\n",extrapOrder);
                OV_ABORT("ERROR");
              }

            const int e=0, c=0; // equation number and component number 
            ForBoundary(side,axis)
            {

                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);

        // Set the index-shift for this side
                is1=is2=is3=0;
                isv[axis]=1-2*side;   // +1 on left and -1 on right      

                if( mg.boundaryCondition(side,axis)==dirichlet )
                {
          // ------------ FILL DIRICHLET BC ------------

                    printF("+++++ IMPLICIT BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) DIRICHLET\n",grid,side,axis);
                    coeffLocal(    M,Ib1,Ib2,Ib3) = 0.0;  // zero out any existing equations
                    coeffLocal(mDiag,Ib1,Ib2,Ib3) = 1.0;

          // --- EXTRAPOLATE GHOST LINES ---

                    for( int ghost=1; ghost<=numberOfGhostLines; ghost++ )
                    {
                            getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                            coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0;
                            FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                            {
                                int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                // --- fill in the coefficients of the extrapolation formula ---
                                for( int m=0; m<=extrapOrder; m++ )
                                {
                                    coeffLocal(m,i1m,i2m,i3m) = extrapCoeff[m];
                                    int j1=i1m + m*is1, j2=i2m + m*is2, j3=i3m + m*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                    setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                }
                            } // end FOR_3D
                    } // end for ghost

                }
                else if( mg.boundaryCondition(side,axis)==neumann )
                {
          // ------------ FILL NEUMANN BC ------------

                    printF("+++++ IMPLICIT BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) NEUMANN\n",grid,side,axis);

                    mg.update(MappedGrid::THEvertexBoundaryNormal);
                    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal); 

                    realSerialArray xCoeff(M0,Ib1,Ib2,Ib3), yCoeff(M0,Ib1,Ib2,Ib3), zCoeff; 
                    mgop.assignCoefficients(MappedGridOperators::xDerivative ,xCoeff, Ib1,Ib2,Ib3,0,0);
                    mgop.assignCoefficients(MappedGridOperators::yDerivative ,yCoeff, Ib1,Ib2,Ib3,0,0);
                    if( numberOfDimensions==3 )
                    {
                        zCoeff.redim(M0,Ib1,Ib2,Ib3);
                        mgop.assignCoefficients(MappedGridOperators::zDerivative ,zCoeff, Ib1,Ib2,Ib3,0,0);
                    }

                    FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                    {
                            
                        int i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)

            // Specify that this a "real" equation on the first ghost line: 
            // (A "real" equation has a possible non-zero right-hand-side)
                        setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost1);              

                        ForStencil(m1,m2,m3)
                        {
                            int m  = M123(m1,m2,m3);        // the single-component coeff-index
                            
                            coeffLocal(m,i1m,i2m,i3m) = normal(i1,i2,i3,0)*xCoeff(m,i1,i2,i3) + normal(i1,i2,i3,1)*yCoeff(m,i1,i2,i3);
                            if( numberOfDimensions==3 )
                                coeffLocal(m,i1m,i2m,i3m) += normal(i1,i2,i3,2)*zCoeff(m,i1,i2,i3);

              // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
                            int j1=i1+m1, j2=i2+m2, j3=i3+m3;                       // the stencil is centred on the boundary pt (i1,i2,i3)
                            setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                      }

                    } // end FOR_3D

          // fill ghost 2 with extrapolation
                    for( int ghost=2; ghost<=numberOfGhostLines; ghost++ )
                    {
                            getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                            coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0;
                            FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                            {
                                int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                // --- fill in the coefficients of the extrapolation formula ---
                                for( int m=0; m<=extrapOrder; m++ )
                                {
                                    coeffLocal(m,i1m,i2m,i3m) = extrapCoeff[m];
                                    int j1=i1m + m*is1, j2=i2m + m*is2, j3=i3m + m*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                    setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                }
                            } // end FOR_3D
                    } // end for ghost


                }
                else if(  mg.boundaryCondition(side,axis)> 0 )
                {
                    printF("fill implicit matrix:ERROR: unknown boundaryCondition=%d \n",mg.boundaryCondition(side,axis));
                    OV_ABORT("error");

                }
            }

          // // Evaluate the (single component) Laplace operator for points on the boundary
          // realSerialArray xxCoeff(M0,Ib1,Ib2,Ib3), yyCoeff(M0,Ib1,Ib2,Ib3), xCoeff(M0,Ib1,Ib2,Ib3), yCoeff(M0,Ib1,Ib2,Ib3), idCoeff(M0,Ib1,Ib2,Ib3);
          // mgop.assignCoefficients(MappedGridOperators::xxDerivative,xxCoeff,Ib1,Ib2,Ib3,0,0); //
          // mgop.assignCoefficients(MappedGridOperators::yyDerivative,yyCoeff,Ib1,Ib2,Ib3,0,0); //
          // mgop.assignCoefficients(MappedGridOperators::xDerivative ,xCoeff, Ib1,Ib2,Ib3,0,0);
          // mgop.assignCoefficients(MappedGridOperators::yDerivative ,yCoeff, Ib1,Ib2,Ib3,0,0);
          // mgop.assignCoefficients(MappedGridOperators::identityOperator,idCoeff,Ib1,Ib2,Ib3,0,0);

            if( addUpwinding && debug>3  )
            {
        // ::display(impCoeff[grid],sPrintF("implicit time-stepping matrix on grid=%d",grid));
                displayCoeff(impCoeff[grid],sPrintF("implicit time-stepping matrix on grid=%d",grid));
        // OV_ABORT("stop here for now");
            }

        }

        impCoeff.finishBoundaryConditions(); 

        impSolver.setCoefficientArray( impCoeff );   // supply coefficients to Oges

    }
    timing(timeForInitialize) += getCPU()-cpu0;

  
    return 0;
}
