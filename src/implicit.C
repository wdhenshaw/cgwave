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
                                const int&gridIndexRange, const int& dimRange, const int &isPeriodic, real&u, const real&un, const int&mask,
                                const real&rsxy, const real&xy, real &uTemp1, real & uTemp2, 
                                const int&boundaryCondition, const real & frequencyArray,
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

    const int & debug                    = dbase.get<int>("debug");
    FILE *& debugFile                    = dbase.get<FILE*>("debugFile");
    FILE *& pDebugFile                   = dbase.get<FILE*>("pDebugFile");

    const Real & c                       = dbase.get<real>("c");
    const real & dt                      = dbase.get<real>("dt");
    const int & orderOfAccuracy          = dbase.get<int>("orderOfAccuracy");

    const int & numberOfFrequencies      = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray     = dbase.get<RealArray>("frequencyArray");
    const RealArray & frequencyArraySave = dbase.get<RealArray>("frequencyArraySave");  
    const int & solveHelmholtz           = dbase.get<int>("solveHelmholtz");
    const int & computeEigenmodes        = dbase.get<int>("computeEigenmodes");

    const int & upwind                   = dbase.get<int>("upwind");
    const int & implicitUpwind           = dbase.get<int>("implicitUpwind");
    bool useUpwindDissipation            = upwind;
    int & totalImplicitIterations        = dbase.get<int>("totalImplicitIterations");  
    int & totalImplicitSolves            = dbase.get<int>("totalImplicitSolves");  

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


  // uc = 0.; // ** TEST **********************************************************

    const int assignBCForImplicit = 1;  // used in call to bcOptWave

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        const IntegerArray & gid = mg.gridIndexRange();
        
        OV_GET_SERIAL_ARRAY(Real,rhs[grid],rhsLocal);
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);


        ForBoundary(side,axis) 
        {
            if( mg.boundaryCondition(side,axis)==exactBC )
            {
        // -- assign boundary and ghost values to be exact ---
                if( t<=2*dt )
                    printF("takeImplicitStep: Set knownSolution at exact boundary (side,axis,grid)=(%d,%d,%d)\n",side,axis,grid);

                getBoundaryIndex(gid,side,axis,I1,I2,I3);
        // int extra=0;  // = numGhost 
        // if( side==0 )
        //   Iv[axis] = Range(gid(0,axis)-extra,gid(0,axis));
        // else 
        //   Iv[axis] = Range(gid(1,axis),gid(1,axis)+extra);

                getUserDefinedKnownSolution( t, grid, rhs[grid], I1,I2,I3 );

        // should we set values at ghost to zero??

            }
        }

    // get parameters for calling fortran
        OV_GET_SERIAL_ARRAY(Real,uc[grid],ucLocal); // *CHECK ME*
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
                else
                {
                    knownSolutionOption=1000;  // all other user defined known solutions
                }   
            }
            int addForcingBC= forcingOption==noForcing ? 0 : 1;
            if( applyKnownSolutionAtBoundaries )
                addForcingBC=1; 
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
                addForcingBC,                    // ipar(12)
                forcingOption,                   // ipar(13)
                useUpwindDissipation,            // ipar(14)
                numGhost,                        // ipar(15)
                assignBCForImplicit,             // ipar(16)
                bcApproach,                      // ipar(17)
                numberOfFrequencies              // ipar(18)
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
                    if( frequencyArraySave(0)==0 )
                        printF("WARNING: bcOpt: frequencyArraySave(0)=%12.4e. NOT ADJUSTING c for EM2 absorbing BC\n",frequencyArraySave(0));
                    if( dt<=0 )
                        printF("WARNING: bcOpt: dt<= 0 ! dt=%12.4e. NOT ADJUSTING c for EM2 absorbing BC\n",dt);
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
            real *pu = rhsLocal.getDataPointer();
            real *pun = ucLocal.getDataPointer();
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

        int ierr=0;
    // if( true )
    //   ::display(bcLocal,"implicit: bcLocal","%3i ");

        bcOptWave(mg.numberOfDimensions(),
                            rhsLocal.getBase(0),rhsLocal.getBound(0),rhsLocal.getBase(1),rhsLocal.getBound(1),
                            rhsLocal.getBase(2),rhsLocal.getBound(2),
                            indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),
                            *rhsLocal.getDataPointer(), *pun, *pmask, *prsxy, *pxy, *puTemp1, *puTemp2,
                            bcLocal(0,0), frequencyArray(0),
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

        if( debug & 2 || (t <= 10.*dt) )
            printF("CgWave::takeImplicitStep: implicit solve: max-res= %8.2e (iterations=%i) ***\n",
                                  impSolver.getMaximumResidual(),numIterations);
    }
    else
    {
        if( false )
        {
            rhs.display("RHS before implicit solve");
        }
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

    bool checkResiduals=false;
    if( checkResiduals )
    {
    // ---- CHECK SOME RESIDUALS ---
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const bool isRectangular = mg.isRectangular();

            OV_GET_SERIAL_ARRAY(Real,uc[grid],ucLocal);
            OV_GET_SERIAL_ARRAY(Real,un[grid],unLocal);

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
                {
                    dr[dir]=mg.gridSpacing(dir);   
                    dx[dir] = dr[dir];       
                }
            }      

            
            Real ca = c; // c*0;

            ForBoundary(side,axis)
            {
                int is = 1-2*side;
                if( mg.boundaryCondition(side,axis)==absorbing )
                {
                    getBoundaryIndex(mg.gridIndexRange(),side,axis,I1,I2,I3);
          // Engquist-Majda order2 scheme
          //  D+t (-D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : left 
          //  D+t ( D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : right 
  
                      RealArray unx(I1,I2,I3), unxx(I1,I2,I3), unyy(I1,I2,I3);
                      unxx(I1,I2,I3) = (unLocal(I1+1,I2,I3)-2.*unLocal(I1,I2,I3)+unLocal(I1-1,I2,I3))/(dx[0]*dx[0]);
                      unyy(I1,I2,I3) = (unLocal(I1,I2+1,I3)-2.*unLocal(I1,I2,I3)+unLocal(I1,I2-1,I3))/(dx[1]*dx[1]);

                      RealArray ucx(I1,I2,I3), ucxx(I1,I2,I3), ucyy(I1,I2,I3);
                      ucxx(I1,I2,I3) = (ucLocal(I1+1,I2,I3)-2.*ucLocal(I1,I2,I3)+ucLocal(I1-1,I2,I3))/(dx[0]*dx[0]);
                      ucyy(I1,I2,I3) = (ucLocal(I1,I2+1,I3)-2.*ucLocal(I1,I2,I3)+ucLocal(I1,I2-1,I3))/(dx[1]*dx[1]); 
                      
                      RealArray res(I1,I2,I3);
  
             // .5*(-2.*c/dxa2)
                      if( axis==0 ) 
                      {             
                          unx(I1,I2,I3) = (unLocal(I1+1,I2,I3) -unLocal(I1-1,I2,I3))/(2.*dx[0]);
                          ucx(I1,I2,I3) = (ucLocal(I1+1,I2,I3) -ucLocal(I1-1,I2,I3))/(2.*dx[0]);
                          res = -is*(unx-ucx)/dt  + (.5*ca)*( unxx + ucxx ) + (.25*ca)*( unyy + ucyy );
                      }
                      else
                      {
                          unx(I1,I2,I3) = (unLocal(I1,I2+1,I3) -unLocal(I1,I2-1,I3))/(2.*dx[1]);
                          ucx(I1,I2,I3) = (ucLocal(I1,I2+1,I3) -ucLocal(I1,I2-1,I3))/(2.*dx[1]);
                          res = -is*(unx-ucx)/dt  + (.5*ca)*( unyy + ucyy ) + (.25*ca)*( unxx + ucxx );
                      }

                    Real maxRes = max(fabs(res));
                    printF("AFTER IMP solve: (side,axis,grid)=(%d,%d,%d) c=%12.4e, bc=absorbing maxRes=%9.2e\n",side,axis,grid,c,maxRes); 
                    ::display(res,"res","%9.2e ");
                        
                }
            }
        }
        OV_ABORT("STOP here for now");
    } // end if checkResiduals

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
// Declare variables needed for Lap^2 in curvilinear coordinates
// Input: 
//   R = r,s,t
//   N = 0,1,2  
// ============================================================================================  


// ============================================================================================
// Fill the upwind dissipation equations into the matrix
//             *OLD WAY* 
//   ====  Use smaller stencils where the wide stenci does not fit ====
///  
// ============================================================================================  


// ==========================================================================
// Macro: fill the matrix with extrapolation formula
//
// Input: eq,uc : equation and component number 
// Input: i1g,i2g,i3g : point to extrapolate (unused point)
// Input: is1,is2,is3 : direction to extrpolate 
// 
// ==========================================================================


// ============================================================================================
// Macro: Fill the upwind dissipation equations into the matrix
//             *NEW WAY* 
//   ====  Extrapolate unused points in the stencil ====
///  
// ============================================================================================  


// ============================================================================================
/// \brief Form and the matrix for implicit time-stepping
// ============================================================================================
int CgWave::formImplicitTimeSteppingMatrix()
{
    real cpu0=getCPU();

    int & debug                          = dbase.get<int>("debug");
    const Real & c                       = dbase.get<real>("c");
    const real & dt                      = dbase.get<real>("dt");
    const Real & damp                    = dbase.get<Real>("damp");
    const Real & dampSave                = dbase.get<Real>("dampSave");

    const int & solveHelmholtz           = dbase.get<int>("solveHelmholtz");
    const RealArray & frequencyArray     = dbase.get<RealArray>("frequencyArray");
    const RealArray & frequencyArraySave = dbase.get<RealArray>("frequencyArraySave");    

    const int & useSuperGrid             = dbase.get<int>("useSuperGrid");
    const IntegerArray & superGrid       = dbase.get<IntegerArray>("superGrid");  

    const int & orderOfAccuracy          = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime    = dbase.get<int>("orderOfAccuracyInTime");
    const IntegerArray & gridIsImplicit  = dbase.get<IntegerArray>("gridIsImplicit");

    const int & upwind                   = dbase.get<int>("upwind");
    const int & implicitUpwind           = dbase.get<int>("implicitUpwind");

    const bool addUpwinding = upwind && implicitUpwind;

    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");

  // addUpwinding=false; // *********** TURN OFF FOR NOW ***************


    printF("\n ==================== FORM MATRIX FOR IMPLICIT TIME-STEPPING ===================\n");
    printF("   c=%.4g, dt=%9.3e, orderOfAccuracy=%d, orderOfAccuracyInTime=%d  bcApproach=%d\n",
              c,dt,orderOfAccuracy,orderOfAccuracyInTime,(int)bcApproach);
    printF(" upwind=%d, add upwinding to implicit matrix=%d, useSuperGrid=%d\n",upwind,addUpwinding,useSuperGrid);
    if( damp!=0 )
        printF(" Linear damping is on: damp=%14.6e, dampSave=%14.6e\n",damp,dampSave);
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


        if( numberOfDimensions==3 )
        {
      // *wdh* July 5, 2023

            Real fillinRatio=0;
      //solver.get(OgesParameters::THEfillinRatio,fillinRatio );
      // printF("Current fillinRatio = %g\n",fillinRatio);
            int stencilWidth = orderOfAccuracy+1;
            int stencilSize=int( pow(stencilWidth,cg.numberOfDimensions())+1 );  // add 1 for interpolation equations

            fillinRatio = (stencilSize+2)*2 + 20.; // what should this be ?
            printF("formImplicitTimeSteppingMatrix: Set new fillinRatio = %g\n",fillinRatio);
            impSolver.set(OgesParameters::THEfillinRatio,fillinRatio );
        }
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
  // bool usePredefined = useMultigrid && !addUpwinding && orderOfAccuracyInTime==2; // ** FIX ME for no-MG
    bool usePredefined = !addUpwinding && orderOfAccuracyInTime==2; // ** FIX ME for no-MG


  // We cannot use predefined if there are Cartesian grids with superGrid
    bool isAllRectangular=true;
    bool isAllCurvilinear=true;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        bool isRectangular = mg.isRectangular();
        isAllRectangular = isAllRectangular && isRectangular;
        isAllCurvilinear  = isAllCurvilinear &&  !isRectangular;

        usePredefined = usePredefined && ( !superGrid(grid) || !isRectangular );

        ForBoundary(side,axis)
        {
      // Check for special BCs
            if( mg.boundaryCondition(side,axis)==CgWave::absorbing ||  
                    mg.boundaryCondition(side,axis)==CgWave::abcEM2 )
            {
              usePredefined=false;
            }

        }    
    } 

  // We cannot use predefined if there are Cartesian grids with superGrid

    

  // if( orderOfAccuracyInTime != 2  )
  // {
  //   printF("CgWave::formImplicitTimeSteppingMatrix: Error orderOfAccuracyInTime != 2 not implemented yet.\n");
  //   OV_ABORT("ERROR");    
  // }

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

            constantCoeff(0,grid) = 1. + damp*dt*.5;
            constantCoeff(1,grid) = - alpha*SQR(c*dt);
        }

    // Assign boundary conditions for Oges
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            ForBoundary(side,axis)
            {
                  if( mg.boundaryCondition(side,axis)==CgWave::dirichlet ||
             // mg.boundaryCondition(side,axis)==CgWave::absorbing ||    // ** DO THIS FOR NOW : absorbing terminated with Dirichlet
                          mg.boundaryCondition(side,axis)==CgWave::exactBC )
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


    // if( useSuperGrid )
    // {
    //   OV_ABORT("implicit: FINISH ME FOR SuperGrid");
    // }

    // Here are the coefficients in the upwind dissipation operator (D+D-)^p 
    // There are 4 cases depending on whether the full wider stencil is available

    // fourth-order dissipation for 2nd-order scheme:
    // Real upwindCoeff4[4][5] = { 1.,-4.,6.,-4.,1.,
    //                             1.,-3.,3.,-1.,0.,   // extrap right-most point D-^3 u(2)
    //                             0.,-1.,3.,-3.,1.,   // extrap left -most point D+^3 u(-2)
    //                             0.,-1.,2.,-1.,0.
    //                           };

    // *wdh* July 3, 2023 change sign of upwindCoeff so we don't need to flip sign of upwind-coeff
        Real upwindCoeff4[4][5] = { -1.,4.,-6.,+4.,-1.,
                                                                -1.,3.,-3.,+1., 0.,   // extrap right-most point D-^3 u(2)
                                                                  0.,1.,-3.,+3.,-1.,   // extrap left -most point D+^3 u(-2)
                                                                  0.,1.,-2.,+1., 0.
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

    // Real adSosup = c*dt/( sqrt(1.*numberOfDimensions) * pow(2.,(orderOfAccuracy+1)) );

        if( orderOfAccuracy==2 )
        {
            for( int m=0; m<4; m++ )
                upwindCoeff[m] =upwindCoeff4[m];
        }
        else if( orderOfAccuracy==4 )
        {
            for( int m=0; m<4; m++ )
                upwindCoeff[m] =upwindCoeff6[m];      
        }
        else if( orderOfAccuracy==6 )
        {
            OV_ABORT("ERROR orderOfAccuracy");
        }
        else
        {
            OV_ABORT("ERROR orderOfAccuracy");
        }

        if( addUpwinding )
            printF("\n >>>>>>> CgWave::implicitSolver: USE IMPLICIT UPWIND DISS <<<<<<< \n\n");

    // const int cc=0; // component number
        const int e=0, cc=0; // equation number and component number 

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

    // TROUBLE FOR UPWIND CASE -- Operators probably base matrix size on orderOfAccuracy !!  *** FIX ME ***
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
        int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];
        int jv[3], &j1=jv[0], &j2=jv[1], &j3=jv[2];
        int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];
        int m1,m2,m3; 

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid &mg = cg[grid];
            const bool isRectangular = mg.isRectangular();

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
                intArray & equationNumber = sparse.equationNumber;
                intArray & classify = sparse.classify;
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

      // const bool isRectangular = mg.isRectangular();
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
                {
                    dr[dir]=mg.gridSpacing(dir);   
                    dx[dir] = dr[dir];       
                }
            }      


      // --- FILL INTERIOR EQUATIONS ----
      // Solve constCoeff(0,grid)*I +constCoeff(1,grid)*Laplacian 
      // We solve:  I - cImp(-1,0) * (c^2*dt^2)* Delta = ...

            const int mDiag = M123(0,0,0);              // index of diagonal entry

            if( gridIsImplicit(grid)==1 )
            {
        // ----- this grid is adavnced with IMPLICIT time-stepping ----
  

                getIndex(mg.gridIndexRange(),I1,I2,I3);
                RealArray lapCoeff(M0,I1,I2,I3);   

                if( !isRectangular || !superGrid(grid) )
                {
                    mgop.assignCoefficients(MappedGridOperators::laplacianOperator,lapCoeff,I1,I2,I3,0,0); // 
                }
                else
                {
          // ------ rectangular and SuperGrid ----------
                    printF("&&&&&& formImplicitTimeSteppingMatrix: Adjust Cartesian grid=%d for superGrid, addUpwinding=%d &&&&&&&\n",grid,(int)addUpwinding);

                    IntegerArray & useAbsorbingLayer = dbase.get<IntegerArray>("useAbsorbingLayer");

          // --- Coefficients for xx and x derivatives ---
                    RealArray ddCoeff(M0,I1,I2,I3), dCoeff(M0,I1,I2,I3);
                    mgop.assignCoefficients( MappedGridOperators::xxDerivative,ddCoeff,I1,I2,I3,0,0);
                    mgop.assignCoefficients( MappedGridOperators::xDerivative,  dCoeff,I1,I2,I3,0,0);

                    if( useAbsorbingLayer(0,grid) )
                    {
            // -- scale coefficients using superGrid functions --
                        RealArray *& etaxSuperGrid = dbase.get<RealArray*>("etaxSuperGrid" );
                        RealArray & etax = etaxSuperGrid[grid];
                        for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )
                        {
                            for( int m=M0.getBase(); m<=M0.getBound(); m++ )
              // for( int m=0; m<stencilSize; m++ )
                            {
                                ddCoeff(m,i1,I2,I3) *= etax(i1,0);  // scale by "(r.x)^2"
                                dCoeff(m,i1,I2,I3)  *= etax(i1,1);  // scale by "r.xx"
                            }
                        }
                    }

                    lapCoeff(M0,I1,I2,I3) = ddCoeff(M0,I1,I2,I3)+ dCoeff(M0,I1,I2,I3);    // transformed .xx derivative

          // --- Coefficients for yy and y derivatives ---
                    mgop.assignCoefficients( MappedGridOperators::yyDerivative,ddCoeff,I1,I2,I3,0,0);
                    mgop.assignCoefficients( MappedGridOperators::yDerivative,  dCoeff,I1,I2,I3,0,0);

                    if( useAbsorbingLayer(1,grid) )
                    {
            // -- scale coefficients using superGrid functions --  
                        RealArray *& etaySuperGrid = dbase.get<RealArray*>("etaySuperGrid" );
                        RealArray & etay = etaySuperGrid[grid];
                        for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )
                        {
                            for( int m=M0.getBase(); m<=M0.getBound(); m++ )
              // for( int m=0; m<stencilSize; m++ )
                            {
                                ddCoeff(m,I1,i2,I3) *= etay(i2,0);
                                dCoeff(m,I1,i2,I3)  *= etay(i2,1);
                            }
                        }
                    }

                    lapCoeff(M0,I1,I2,I3) += ddCoeff(M0,I1,I2,I3)+ dCoeff(M0,I1,I2,I3);    // transformed .yy derivative

                }
                
                Real ccLap = - cImp(-1,0)*SQR(c*dt);         // note minus
                coeffLocal(M0,I1,I2,I3)  = ccLap*lapCoeff;

        // set diagonal entry
                coeffLocal(mDiag,I1,I2,I3) += 1.0 + damp*.5*dt;


                if( orderOfAccuracy==4 )
                {
           // Add modified equation term 
           //     (L_2h)^2

                    printF("\n ***ADD Fourth-order in time coefficients to the Implicit Matrix ****\n\n");

                    if( isRectangular )
                    {
            // 2D: 
            // (Delta_2h)^2  = (D+xD-x)^2 + (D+yD-y)^2 + 2 (D+xD-x)(D+yD-y)
            // 3D 
            // (Delta_2h)^2  = (D+xD-x)^2 + (D+yD-y)^2 + (D+zD-z)^2 + 2 (D+xD-x)(D+yD-y) + 2 (D+xD-x)(D+zD-z)+ 2 (D+yD-y)(D+zD-z)

                        Range R5 = Range(-2,2);
                        RealArray lapSq(R5,R5,R5);
                        lapSq = 0.;

                        const Real dx4 = dx[0]*dx[0]*dx[0]*dx[0];
                        const Real dy4 = dx[1]*dx[1]*dx[1]*dx[1];
                        lapSq(-2,0,0) +=  1./dx4;  lapSq(0,-2,0) +=  1./dy4;
                        lapSq(-1,0,0) += -4./dx4;  lapSq(0,-1,0) += -4./dy4;
                        lapSq( 0,0,0) +=  6./dx4;  lapSq(0, 0,0) +=  6./dy4;
                        lapSq( 1,0,0) += -4./dx4;  lapSq(0, 1,0) += -4./dy4;
                        lapSq( 2,0,0) +=  1./dx4;  lapSq(0, 2,0) +=  1./dy4;

                        const Real dxy2 = dx[0]*dx[0]*dx[1]*dx[1];
                        lapSq(-1, 1,0) += 2./dxy2; lapSq(0, 1,0) +=-4./dxy2; lapSq(1, 1,0) += 2./dxy2;
                        lapSq(-1, 0,0) +=-4./dxy2; lapSq(0, 0,0) += 8./dxy2; lapSq(1, 0,0) +=-4./dxy2;
                        lapSq(-1,-1,0) += 2./dxy2; lapSq(0,-1,0) +=-4./dxy2; lapSq(1,-1,0) += 2./dxy2;

                        if( numberOfDimensions==3 )
                        {
                            const Real dz4 = dx[2]*dx[2]*dx[2]*dx[2];
                            lapSq(0,0,-2) +=  1./dz4;  
                            lapSq(0,0,-1) += -4./dz4;  
                            lapSq(0,0, 0) +=  6./dz4;  
                            lapSq(0,0, 1) += -4./dz4;  
                            lapSq(0,0, 2) +=  1./dz4; 

                            const Real dxz2 = dx[0]*dx[0]*dx[2]*dx[2];
                            lapSq(-1,0, 1) += 2./dxz2; lapSq(0,0, 1) +=-4./dxz2; lapSq(1,0, 1) += 2./dxz2;
                            lapSq(-1,0, 0) +=-4./dxz2; lapSq(0,0, 0) += 8./dxz2; lapSq(1,0, 0) +=-4./dxz2;
                            lapSq(-1,0,-1) += 2./dxz2; lapSq(0,0,-1) +=-4./dxz2; lapSq(1,0,-1) += 2./dxz2;

                            Real dyz2 = dx[1]*dx[1]*dx[2]*dx[2];
                            lapSq(0,-1, 1) += 2./dyz2; lapSq(0,0, 1) +=-4./dyz2; lapSq(0,1, 1) += 2./dyz2;
                            lapSq(0,-1, 0) +=-4./dyz2; lapSq(0,0, 0) += 8./dyz2; lapSq(0,1, 0) +=-4./dyz2;
                            lapSq(0,-1,-1) += 2./dyz2; lapSq(0,0,-1) +=-4./dyz2; lapSq(0,1,-1) += 2./dyz2;
                        }

            // coefficients in implicit time-stepping  
            //  D+t D-t u =              c^2 Delta( cImp(1,0) *u^{n+1} + cImp(0,0) *u^n + cImp(-1,0)* u^{n-1} )   :  second-order coeff cImp(-1:1,0)
            //              -(c^4*dt^2/12) Delta^2( cImp(1,1) *u^{n+1} + cImp(0,1) *u^n + cImp(-1,1)* u^{n-1}  )  :  foruth-order ceoff cImp(-1:1,1) 
            // For accuracy the weights depend on one parameter beta2 for second-order,
            // and a second parameter beta4 for fourth-order: (See notes in research/timeStepping/implicitTaylorSchemes.pdf)

                        const Real cLapSq = cImp(-1,1)*(c*dt)*(c*dt)*(c*dt)*(c*dt);  // CHECK ME 
                        ForStencil(m1,m2,m3)
                        {
                            int m  = M123(m1,m2,m3);   
                            coeffLocal(m,I1,I2,I3) += cLapSq*lapSq(m1,m2,m3);
                        }

                    }
                    else
                    {
            // OV_ABORT("implicit matrix: finish me for order=4 CURVLINEAR");

                        printF("implicit matrix: order=4 CURVLINEAR : **CHECK ME**\n");

                        OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rxLocal);
            // macro to make the rxLocal array look 5-dimensional 
                        #define DD(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2))     

            // ----- COMPUTE DERIVATIVES OF METRICS -----
                        Range Rd2=numberOfDimensions*numberOfDimensions;
                        RealArray ddx(I1,I2,I3,Rd2), ddy(I1,I2,I3,Rd2);
                        mgop.derivative( MappedGridOperators::xDerivative,rxLocal,ddx,I1,I2,I3,Rd2);            
                        mgop.derivative( MappedGridOperators::yDerivative,rxLocal,ddy,I1,I2,I3,Rd2);

                        const int extra=1; 
                        Index J1,J2,J3;
                        getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
                        RealArray ddxx(J1,J2,J3,Rd2), ddxy(J1,J2,J3,Rd2), ddyy(J1,J2,J3,Rd2);          // COMPUTE THESE AT MORE POINTS FOR BELOW
                        mgop.derivative( MappedGridOperators::xxDerivative,rxLocal,ddxx,J1,J2,J3,Rd2);            
                        mgop.derivative( MappedGridOperators::xyDerivative,rxLocal,ddxy,J1,J2,J3,Rd2);            
                        mgop.derivative( MappedGridOperators::yyDerivative,rxLocal,ddyy,J1,J2,J3,Rd2);             


                        #define DDX(i1,i2,i3,m1,m2) ddx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                        #define DDY(i1,i2,i3,m1,m2) ddy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 

                        #define DDXX(i1,i2,i3,m1,m2) ddxx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                        #define DDXY(i1,i2,i3,m1,m2) ddxy(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                        #define DDYY(i1,i2,i3,m1,m2) ddyy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 




            // Define stencils for parametric derivatives 
                        Range R5 = Range(-2,2);
                        RealArray rrrrCoeff(R5), ssssCoeff(R5);
                        rrrrCoeff = 0.;  ssssCoeff=0.;
                        const Real dr4 = dr[0]*dr[0]*dr[0]*dr[0];
                        const Real ds4 = dr[1]*dr[1]*dr[1]*dr[1];
            // (D+D-)^2 stencil 
                        rrrrCoeff(-2) =  1./dr4;  ssssCoeff(-2) =  1./ds4;
                        rrrrCoeff(-1) = -4./dr4;  ssssCoeff(-1) = -4./ds4;
                        rrrrCoeff( 0) =  6./dr4;  ssssCoeff( 0) =  6./ds4;
                        rrrrCoeff( 1) = -4./dr4;  ssssCoeff( 1) = -4./ds4;
                        rrrrCoeff( 2) =  1./dr4;  ssssCoeff( 2) =  1./ds4;            
  
                        RealArray rrrCoeff(R5), sssCoeff(R5);
                        rrrCoeff=0.; sssCoeff=0.;
            // D0(D+D-) stencil 
                        rrrCoeff(-2) = -1./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff(-2) = -1./(2.*dr[1]*dr[1]*dr[1]);
                        rrrCoeff(-1) = +2./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff(-1) = +2./(2.*dr[1]*dr[1]*dr[1]);
                        rrrCoeff( 0) =  0.;                         sssCoeff( 0) =  0.;    
                        rrrCoeff( 1) = -2./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff( 1) = -2./(2.*dr[1]*dr[1]*dr[1]);
                        rrrCoeff( 2) =  1./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff( 2) =  1./(2.*dr[1]*dr[1]*dr[1]);

                        RealArray rrCoeff(R5), ssCoeff(R5);
                        rrCoeff=0.; ssCoeff=0.;
            // D+D-
                        rrCoeff(-1) = 1./(dr[0]*dr[0]); ssCoeff(-1) = 1./(dr[1]*dr[1]); 
                        rrCoeff( 0) =-2./(dr[0]*dr[0]); ssCoeff( 0) =-2./(dr[1]*dr[1]); 
                        rrCoeff( 1) = 1./(dr[0]*dr[0]); ssCoeff( 1) = 1./(dr[1]*dr[1]); 

                        RealArray rCoeff(R5), sCoeff(R5);
                        rCoeff=0.; sCoeff=0.;
            // Dz stencil 
                        rCoeff(-1) = -1./(2.*dr[0]); sCoeff(-1) = -1./(2.*dr[1]); 
                        rCoeff( 0) =  0.;            sCoeff( 0) =  0.;
                        rCoeff( 1) =  1./(2.*dr[0]); sCoeff( 1) =  1./(2.*dr[1]); 

            // Identity stencil 
                        RealArray iCoeff(R5);
                        iCoeff=0.;
                        iCoeff(0)=1.; 

                        const Real cLapSq = cImp(-1,1)*(c*dt)*(c*dt)*(c*dt)*(c*dt);
                        if( numberOfDimensions==2 )
                        {      
                            RealArray ddxxx(I1,I2,I3,Rd2), ddxyy(I1,I2,I3,Rd2), ddyyy(I1,I2,I3,Rd2);

              // Operators have no third x derivative :(
              // Take x derivative of xx derivative 
                            for( int dir=0; dir<=Rd2.getBound(); dir++ )
                            {
                // rxxx.x
                                ddxxx(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddxx(I1+1,I2,I3,dir) - ddxx(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,0)*( ddxx(I1,I2+1,I3,dir) - ddxx(I1,I2-1,I3,dir) )/(2.*dr[1]);

                // rxyy.x 
                                ddxyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,0)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1]);                                    

                // rxyy.y
                                ddyyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,1)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1]);  
                                                                                                                                                    
                            }  

                            #define DDXXX(i1,i2,i3,m1,m2) ddxxx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                            #define DDXYY(i1,i2,i3,m1,m2) ddxyy(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
                            #define DDYYY(i1,i2,i3,m1,m2) ddyyy(i1,i2,i3,(m1)+numberOfDimensions*(m2))  

                                                            
                            FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the domain
                            {
                                if( maskLocal(i1,i2,i3)>0 )
                                {
                                    Real rx    =    DD(i1,i2,i3,0,0);
                                    Real ry    =    DD(i1,i2,i3,0,1);
                                    Real rxx   =   DDX(i1,i2,i3,0,0);
                                    Real rxy   =   DDY(i1,i2,i3,0,0);
                                    Real ryy   =   DDY(i1,i2,i3,0,1);  // ry.y
                                    Real rxxx  =  DDXX(i1,i2,i3,0,0);
                                    Real rxxy  =  DDXY(i1,i2,i3,0,0);  // rx.xy
                                    Real rxyy  =  DDYY(i1,i2,i3,0,0);  // rx.yy
                                    Real ryyy  =  DDYY(i1,i2,i3,0,1);  // ry.yy
                                    Real rxxxx = DDXXX(i1,i2,i3,0,0);  // rx.xxx
                                    Real rxxyy = DDXYY(i1,i2,i3,0,1);  // ry.xyy
                                    Real ryyyy = DDYYY(i1,i2,i3,0,1);  // ry.yyy   

                                    Real sx    =    DD(i1,i2,i3,1,0);
                                    Real sy    =    DD(i1,i2,i3,1,1);
                                    Real sxx   =   DDX(i1,i2,i3,1,0);
                                    Real sxy   =   DDY(i1,i2,i3,1,0);
                                    Real syy   =   DDY(i1,i2,i3,1,1);  
                                    Real sxxx  =  DDXX(i1,i2,i3,1,0);
                                    Real sxxy  =  DDXY(i1,i2,i3,1,0);
                                    Real sxyy  =  DDYY(i1,i2,i3,1,0);  // sx.yy
                                    Real syyy  =  DDYY(i1,i2,i3,1,1);  
                                    Real sxxxx = DDXXX(i1,i2,i3,1,0);  // sx.xxx
                                    Real sxxyy = DDXYY(i1,i2,i3,1,1);  // sy.xyy
                                    Real syyyy = DDYYY(i1,i2,i3,1,1);  // sy.yyy                                                       

                  // ---- COEFFICIENTS OF LAPLACIAN SQUARED  : from laplacianCoefficients.mpl ----
                                    Real urrrr = pow(rx, 0.4e1) + 0.2e1 * ry * ry * rx * rx + pow(ry, 0.4e1);
                                    Real urrrs = 0.4e1 * pow(rx, 0.3e1) * sx + 0.4e1 * pow(ry, 0.3e1) * sy + 0.2e1 * ry * (sy * rx * rx + 0.2e1 * ry * sx * rx) + 0.2e1 * sy * ry * rx * rx;
                                    Real urrss = 6 * rx * rx * sx * sx + 6 * ry * ry * sy * sy + 2 * ry * (2 * sy * sx * rx + ry * sx * sx) + 2 * sy * (sy * rx * rx + 2 * ry * sx * rx);
                                    Real ursss = 0.4e1 * rx * pow(sx, 0.3e1) + 0.4e1 * ry * pow(sy, 0.3e1) + 0.2e1 * ry * sy * sx * sx + 0.2e1 * sy * (0.2e1 * sy * sx * rx + ry * sx * sx);
                                    Real ussss = pow(sx, 0.4e1) + 0.2e1 * sy * sy * sx * sx + pow(sy, 0.4e1);
                                    Real urrr = 6 * rx * rx * rxx + 6 * ry * ry * ryy + 2 * ry * (2 * rxy * rx + ry * rxx) + 2 * ryy * rx * rx + 4 * ry * rxy * rx;
                                    Real urrs = ry * (3 * syy * ry + 3 * sy * ryy) + 7 * sy * ry * ryy + ry * (2 * syy * ry + 2 * sy * ryy) + syy * ry * ry + 2 * ry * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * sy * (2 * rxy * rx + ry * rxx) + 4 * ryy * sx * rx + 2 * ry * (2 * sxy * rx + 2 * sx * rxy) + 2 * syy * rx * rx + 4 * sy * rxy * rx + rx * (3 * sxx * rx + 3 * sx * rxx) + 7 * sx * rx * rxx + rx * (2 * sxx * rx + 2 * sx * rxx) + sxx * rx * rx;
                                    Real urss = 7 * ry * sy * syy + sy * (3 * syy * ry + 3 * sy * ryy) + ryy * sy * sy + sy * (2 * syy * ry + 2 * sy * ryy) + 2 * ry * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * ryy * sx * sx + 4 * ry * sx * sxy + 4 * syy * sx * rx + 2 * sy * (2 * sxy * rx + 2 * sx * rxy) + 7 * rx * sx * sxx + sx * (3 * sxx * rx + 3 * sx * rxx) + rxx * sx * sx + sx * (2 * sxx * rx + 2 * sx * rxx);
                                    Real usss = 6 * sx * sx * sxx + 6 * sy * sy * syy + 2 * sy * (2 * sx * sxy + sxx * sy) + 2 * syy * sx * sx + 4 * sy * sx * sxy;
                                    Real urr = 4 * rx * rxxx + 4 * rxyy * rx + 3 * rxx * rxx + 2 * ryy * rxx + 4 * ry * rxxy + 4 * rxy * rxy + 4 * ry * ryyy + 3 * ryy * ryy;
                                    Real urs = 4 * sxxx * rx + 4 * sxyy * rx + 6 * sxx * rxx + 2 * syy * rxx + 4 * sx * rxxx + 4 * sy * rxxy + 8 * sxy * rxy + 4 * sx * rxyy + 4 * ry * sxxy + 4 * syyy * ry + 2 * ryy * sxx + 6 * syy * ryy + 4 * sy * ryyy;
                                    Real uss = 4 * sx * sxxx + 4 * sx * sxyy + 3 * sxx * sxx + 2 * sxx * syy + 4 * sxxy * sy + 4 * sxy * sxy + 4 * sy * syyy + 3 * syy * syy;
                                    Real ur = rxxxx + 2 * rxxyy + ryyyy;
                                    Real us = sxxxx + 2 * sxxyy + syyyy;

                  // printF(" (i1,i2)=(%d,%d) urrrr=%g, urrrs=%g, urrss=%g, ursss=%g, ussss=%g, urrr=%g urrs=%g urss=%g urr=%g urs=%g uss=%g ur=%g us=%g\n",
                  //   urrrr,urrrs,urrss,ursss,ussss,urrr,urrs,urss,urss,usss,urr,urs,uss,ur,us);

                                    ForStencil(m1,m2,m3)
                                    {
                                        int m  = M123(m1,m2,m3);   
                                        coeffLocal(m,i1,i2,i3) += 
                                                                        cLapSq*(  urrrr*rrrrCoeff(m1)*   iCoeff(m2) 
                                                                                        + urrrs* rrrCoeff(m1)*   sCoeff(m2)
                                                                                        + urrss*  rrCoeff(m1)*  ssCoeff(m2)
                                                                                        + ursss*   rCoeff(m1)* sssCoeff(m2)
                                                                                        + ussss*   iCoeff(m1)*ssssCoeff(m2)
                                                                                        + urrr * rrrCoeff(m1)*   iCoeff(m2)
                                                                                        + urrs *  rrCoeff(m1)*   sCoeff(m2)
                                                                                        + urss *   rCoeff(m1)*  ssCoeff(m2)
                                                                                        + usss *   iCoeff(m1)* sssCoeff(m2)
                                                                                        + urr  *  rrCoeff(m1)*   iCoeff(m2)
                                                                                        + urs  *   rCoeff(m1)*   sCoeff(m2)
                                                                                        + uss  *   iCoeff(m1)*  ssCoeff(m2)
                                                                                        + ur   *   rCoeff(m1)*   iCoeff(m2)
                                                                                        + us   *   iCoeff(m1)*   sCoeff(m2)
                                                                                        );
                                    }
                                }
                            }

                        }
                        else
                        {
              // --------------- THREE DIMENSIONS ----------------

                            RealArray rxza(I1,I2,I3,Rd2), ddy(I1,I2,I3,Rd2);
                            mgop.derivative( MappedGridOperators::zDerivative,rxLocal,rxza,I1,I2,I3,Rd2);            

                            #define DDZ(i1,i2,i3,m1,m2) rxza(i1,i2,i3,(m1)+numberOfDimensions*(m2))

                            const int extra=1; 
                            Index J1,J2,J3;
                            getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
                            RealArray ddxz(J1,J2,J3,Rd2), ddyz(J1,J2,J3,Rd2), ddzz(J1,J2,J3,Rd2);          // COMPUTE THESE AT MORE POINTS FOR BELOW
                            mgop.derivative( MappedGridOperators::xzDerivative,rxLocal,ddxz,J1,J2,J3,Rd2);            
                            mgop.derivative( MappedGridOperators::yzDerivative,rxLocal,ddyz,J1,J2,J3,Rd2);            
                            mgop.derivative( MappedGridOperators::zzDerivative,rxLocal,ddzz,J1,J2,J3,Rd2);             

                            #define DDXZ(i1,i2,i3,m1,m2) ddxz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                            #define DDYZ(i1,i2,i3,m1,m2) ddyz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                            #define DDZZ(i1,i2,i3,m1,m2) ddzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  

                            RealArray ddxxx(I1,I2,I3,Rd2), ddxyy(I1,I2,I3,Rd2), ddyyy(I1,I2,I3,Rd2);
                            RealArray ddxzz(I1,I2,I3,Rd2), ddzzz(I1,I2,I3,Rd2), ddyzz(I1,I2,I3,Rd2);

              // Operators have no third x derivative :(
              // Take x derivative of xx derivative 
                            for( int dir=0; dir<=Rd2.getBound(); dir++ )
                            {
                // rxxx.x
                                ddxxx(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddxx(I1+1,I2,I3,dir) - ddxx(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,0)*( ddxx(I1,I2+1,I3,dir) - ddxx(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,0)*( ddxx(I1,I2,I3+1,dir) - ddxx(I1,I2,I3-1,dir) )/(2.*dr[2]);

                // rxyy.x 
                                ddxyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,0)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,0)*( ddyy(I1,I2,I3+1,dir) - ddyy(I1,I2,I3-1,dir) )/(2.*dr[2]);


                // rxyy.y
                                ddyyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,1)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,1)*( ddyy(I1,I2,I3+1,dir) - ddyy(I1,I2,I3-1,dir) )/(2.*dr[2]);

                // rxzz.x 
                                ddxzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,0)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,0)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);

                // rxzz.y 
                                ddyzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,1)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,1)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);

                // rxzz.z
                                ddzzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,2)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                        +DD(I1,I2,I3,1,2)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                        +DD(I1,I2,I3,2,2)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);

                            } 

                            #define DDXXX(i1,i2,i3,m1,m2) ddxxx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                            #define DDXYY(i1,i2,i3,m1,m2) ddxyy(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
                            #define DDYYY(i1,i2,i3,m1,m2) ddyyy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 
                            #define DDXZZ(i1,i2,i3,m1,m2) ddxzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
                            #define DDZZZ(i1,i2,i3,m1,m2) ddzzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
                            #define DDYZZ(i1,i2,i3,m1,m2) ddyzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))   

                            RealArray ttttCoeff(R5);
                            ttttCoeff = 0.; 
                            const Real dt4 = dr[2]*dr[2]*dr[2]*dr[2];
              // (D+D-)^2 stencil 
                            ttttCoeff(-2) =  1./dt4; 
                            ttttCoeff(-1) = -4./dt4; 
                            ttttCoeff( 0) =  6./dt4; 
                            ttttCoeff( 1) = -4./dt4; 
                            ttttCoeff( 2) =  1./dt4;             
      
                            RealArray tttCoeff(R5);
                            tttCoeff=0.; 
              // D0(D+D-) stencil 
                            tttCoeff(-2) = -1./(2.*dr[2]*dr[2]*dr[2]); 
                            tttCoeff(-1) = +2./(2.*dr[2]*dr[2]*dr[2]); 
                            tttCoeff( 0) =  0.;                        
                            tttCoeff( 1) = -2./(2.*dr[2]*dr[2]*dr[2]); 
                            tttCoeff( 2) =  1./(2.*dr[2]*dr[2]*dr[2]); 

                            RealArray ttCoeff(R5);
                            ttCoeff=0.;
              // D+D-
                            ttCoeff(-1) = 1./(dr[2]*dr[2]);
                            ttCoeff( 0) =-2./(dr[2]*dr[2]);
                            ttCoeff( 1) = 1./(dr[2]*dr[2]);

                            RealArray tCoeff(R5);
                            tCoeff=0.; 
              // Dz stencil 
                            tCoeff(-1) = -1./(2.*dr[2]);
                            tCoeff( 0) =  0.;           
                            tCoeff( 1) =  1./(2.*dr[2]);            

                            FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the domain
                            {
                                if( maskLocal(i1,i2,i3)>0 )
                                {
                                        Real rx     =     DD(i1,i2,i3,0,0);
                                        Real ry     =     DD(i1,i2,i3,0,1);
                                        Real rz     =     DD(i1,i2,i3,0,2);
                                        Real rxx    =    DDX(i1,i2,i3,0,0);
                                        Real rxy    =    DDX(i1,i2,i3,0,1);  // .xy 
                                        Real rxz    =    DDX(i1,i2,i3,0,2);  // .xz
                                        Real ryy    =    DDY(i1,i2,i3,0,1);  // .yy
                                        Real ryz    =    DDY(i1,i2,i3,0,2);  // .yz                 
                                        Real rzz    =    DDZ(i1,i2,i3,0,2);  // .zz
                                        Real rxxx   =   DDXX(i1,i2,i3,0,0);  // .xxx
                                        Real rxxy   =   DDXX(i1,i2,i3,0,1);  // .xxy
                                        Real rxyy   =   DDYY(i1,i2,i3,0,0);  // .xyy
                                        Real ryyy   =   DDYY(i1,i2,i3,0,1);  // .yyy
                                        Real rxxz   =   DDXX(i1,i2,i3,0,2);  // .xxz
                                        Real rxzz   =   DDZZ(i1,i2,i3,0,1);  // .xzz
                                        Real rzzz   =   DDZZ(i1,i2,i3,0,2);  // .zzz
                                        Real ryyz   =   DDYY(i1,i2,i3,0,2);  // .yyz
                                        Real ryzz   =   DDZZ(i1,i2,i3,0,1);  // .yzz
                                        Real rxxxx  =  DDXXX(i1,i2,i3,0,0);  // .xxxx
                                        Real rxxyy  =  DDXYY(i1,i2,i3,0,0);  // .xxyy
                                        Real ryyyy  =  DDYYY(i1,i2,i3,0,1);  // .yyy
                                        Real rxxzz  =  DDXZZ(i1,i2,i3,0,0);  // .xxzz
                                        Real rzzzz  =  DDZZZ(i1,i2,i3,0,2);  // .zzzz
                                        Real ryyzz  =  DDYZZ(i1,i2,i3,0,1);  // .yyzz                        
                                        Real sx     =     DD(i1,i2,i3,1,0);
                                        Real sy     =     DD(i1,i2,i3,1,1);
                                        Real sz     =     DD(i1,i2,i3,1,2);
                                        Real sxx    =    DDX(i1,i2,i3,1,0);
                                        Real sxy    =    DDX(i1,i2,i3,1,1);  // .xy 
                                        Real sxz    =    DDX(i1,i2,i3,1,2);  // .xz
                                        Real syy    =    DDY(i1,i2,i3,1,1);  // .yy
                                        Real syz    =    DDY(i1,i2,i3,1,2);  // .yz                 
                                        Real szz    =    DDZ(i1,i2,i3,1,2);  // .zz
                                        Real sxxx   =   DDXX(i1,i2,i3,1,0);  // .xxx
                                        Real sxxy   =   DDXX(i1,i2,i3,1,1);  // .xxy
                                        Real sxyy   =   DDYY(i1,i2,i3,1,0);  // .xyy
                                        Real syyy   =   DDYY(i1,i2,i3,1,1);  // .yyy
                                        Real sxxz   =   DDXX(i1,i2,i3,1,2);  // .xxz
                                        Real sxzz   =   DDZZ(i1,i2,i3,1,1);  // .xzz
                                        Real szzz   =   DDZZ(i1,i2,i3,1,2);  // .zzz
                                        Real syyz   =   DDYY(i1,i2,i3,1,2);  // .yyz
                                        Real syzz   =   DDZZ(i1,i2,i3,1,1);  // .yzz
                                        Real sxxxx  =  DDXXX(i1,i2,i3,1,0);  // .xxxx
                                        Real sxxyy  =  DDXYY(i1,i2,i3,1,0);  // .xxyy
                                        Real syyyy  =  DDYYY(i1,i2,i3,1,1);  // .yyy
                                        Real sxxzz  =  DDXZZ(i1,i2,i3,1,0);  // .xxzz
                                        Real szzzz  =  DDZZZ(i1,i2,i3,1,2);  // .zzzz
                                        Real syyzz  =  DDYZZ(i1,i2,i3,1,1);  // .yyzz                        
                                        Real tx     =     DD(i1,i2,i3,2,0);
                                        Real ty     =     DD(i1,i2,i3,2,1);
                                        Real tz     =     DD(i1,i2,i3,2,2);
                                        Real txx    =    DDX(i1,i2,i3,2,0);
                                        Real txy    =    DDX(i1,i2,i3,2,1);  // .xy 
                                        Real txz    =    DDX(i1,i2,i3,2,2);  // .xz
                                        Real tyy    =    DDY(i1,i2,i3,2,1);  // .yy
                                        Real tyz    =    DDY(i1,i2,i3,2,2);  // .yz                 
                                        Real tzz    =    DDZ(i1,i2,i3,2,2);  // .zz
                                        Real txxx   =   DDXX(i1,i2,i3,2,0);  // .xxx
                                        Real txxy   =   DDXX(i1,i2,i3,2,1);  // .xxy
                                        Real txyy   =   DDYY(i1,i2,i3,2,0);  // .xyy
                                        Real tyyy   =   DDYY(i1,i2,i3,2,1);  // .yyy
                                        Real txxz   =   DDXX(i1,i2,i3,2,2);  // .xxz
                                        Real txzz   =   DDZZ(i1,i2,i3,2,1);  // .xzz
                                        Real tzzz   =   DDZZ(i1,i2,i3,2,2);  // .zzz
                                        Real tyyz   =   DDYY(i1,i2,i3,2,2);  // .yyz
                                        Real tyzz   =   DDZZ(i1,i2,i3,2,1);  // .yzz
                                        Real txxxx  =  DDXXX(i1,i2,i3,2,0);  // .xxxx
                                        Real txxyy  =  DDXYY(i1,i2,i3,2,0);  // .xxyy
                                        Real tyyyy  =  DDYYY(i1,i2,i3,2,1);  // .yyy
                                        Real txxzz  =  DDXZZ(i1,i2,i3,2,0);  // .xxzz
                                        Real tzzzz  =  DDZZZ(i1,i2,i3,2,2);  // .zzzz
                                        Real tyyzz  =  DDYZZ(i1,i2,i3,2,1);  // .yyzz                        

                  // ---- COEFFICIENTS OF 3D LAPLACIAN SQUARED : from laplacianCoefficients.mpl ----
                                    Real urrrr = pow(rx, 0.4e1) + 0.2e1 * ry * ry * rx * rx + 0.2e1 * rz * rz * rx * rx + pow(ry, 0.4e1) + 0.2e1 * rz * rz * ry * ry + pow(rz, 0.4e1);
                                    Real urrrs = 0.4e1 * pow(rx, 0.3e1) * sx + 0.4e1 * pow(rz, 0.3e1) * sz + 0.2e1 * rz * (sz * ry * ry + 0.2e1 * rz * sy * ry) + 0.2e1 * sz * rz * ry * ry + 0.2e1 * rz * (sz * rx * rx + 0.2e1 * rz * sx * rx) + 0.2e1 * sz * rz * rx * rx + 0.4e1 * pow(ry, 0.3e1) * sy + 0.2e1 * ry * (sy * rx * rx + 0.2e1 * ry * sx * rx) + 0.2e1 * sy * ry * rx * rx;
                                    Real urrss = 6 * rx * rx * sx * sx + 6 * rz * rz * sz * sz + 2 * rz * (2 * sz * sy * ry + rz * sy * sy) + 2 * sz * (sz * ry * ry + 2 * rz * sy * ry) + 2 * rz * (2 * sz * sx * rx + rz * sx * sx) + 2 * sz * (sz * rx * rx + 2 * rz * sx * rx) + 6 * ry * ry * sy * sy + 2 * ry * (2 * sy * sx * rx + ry * sx * sx) + 2 * sy * (sy * rx * rx + 2 * ry * sx * rx);
                                    Real ursss = 0.4e1 * rx * pow(sx, 0.3e1) + 0.4e1 * rz * pow(sz, 0.3e1) + 0.2e1 * rz * sz * sy * sy + 0.2e1 * sz * (0.2e1 * sz * sy * ry + rz * sy * sy) + 0.2e1 * rz * sz * sx * sx + 0.2e1 * sz * (0.2e1 * sz * sx * rx + rz * sx * sx) + 0.4e1 * ry * pow(sy, 0.3e1) + 0.2e1 * ry * sy * sx * sx + 0.2e1 * sy * (0.2e1 * sy * sx * rx + ry * sx * sx);
                                    Real ussss = pow(sx, 0.4e1) + 0.2e1 * sy * sy * sx * sx + 0.2e1 * sz * sz * sx * sx + pow(sy, 0.4e1) + 0.2e1 * sz * sz * sy * sy + pow(sz, 0.4e1);
                                    Real urrrt = 0.4e1 * pow(rz, 0.3e1) * tz + 0.2e1 * rz * (tz * ry * ry + 0.2e1 * rz * ty * ry) + 0.2e1 * tz * rz * ry * ry + 0.2e1 * rz * (tz * rx * rx + 0.2e1 * rz * tx * rx) + 0.2e1 * tz * rz * rx * rx + 0.4e1 * pow(ry, 0.3e1) * ty + 0.2e1 * ry * (ty * rx * rx + 0.2e1 * ry * tx * rx) + 0.2e1 * ty * ry * rx * rx + 0.4e1 * pow(rx, 0.3e1) * tx;
                                    Real urrtt = 6 * rz * rz * tz * tz + 2 * rz * (2 * tz * ty * ry + rz * ty * ty) + 2 * tz * (tz * ry * ry + 2 * rz * ty * ry) + 2 * rz * (2 * tz * tx * rx + rz * tx * tx) + 2 * tz * (tz * rx * rx + 2 * rz * tx * rx) + 2 * ry * (2 * ty * tx * rx + ry * tx * tx) + 2 * ty * (ty * rx * rx + 2 * ry * tx * rx) + 6 * ry * ry * ty * ty + 6 * rx * rx * tx * tx;
                                    Real urttt = 0.4e1 * rz * pow(tz, 0.3e1) + 0.2e1 * rz * tz * ty * ty + 0.2e1 * tz * (0.2e1 * tz * ty * ry + rz * ty * ty) + 0.2e1 * rz * tz * tx * tx + 0.2e1 * tz * (0.2e1 * tz * tx * rx + rz * tx * tx) + 0.2e1 * ry * ty * tx * tx + 0.2e1 * ty * (0.2e1 * ty * tx * rx + ry * tx * tx) + 0.4e1 * ry * pow(ty, 0.3e1) + 0.4e1 * rx * pow(tx, 0.3e1);
                                    Real utttt = pow(tx, 0.4e1) + 0.2e1 * ty * ty * tx * tx + 0.2e1 * tz * tz * tx * tx + pow(ty, 0.4e1) + 0.2e1 * tz * tz * ty * ty + pow(tz, 0.4e1);
                                    Real ussst = 0.4e1 * pow(sz, 0.3e1) * tz + 0.2e1 * sz * (tz * sy * sy + 0.2e1 * sz * ty * sy) + 0.2e1 * tz * sz * sy * sy + 0.2e1 * sz * (tz * sx * sx + 0.2e1 * sz * tx * sx) + 0.2e1 * tz * sz * sx * sx + 0.2e1 * sy * (ty * sx * sx + 0.2e1 * sy * tx * sx) + 0.2e1 * ty * sy * sx * sx + 0.4e1 * pow(sy, 0.3e1) * ty + 0.4e1 * pow(sx, 0.3e1) * tx;
                                    Real usstt = 6 * sz * sz * tz * tz + 2 * sz * (2 * tz * ty * sy + sz * ty * ty) + 2 * tz * (tz * sy * sy + 2 * sz * ty * sy) + 2 * sz * (2 * tz * tx * sx + sz * tx * tx) + 2 * tz * (tz * sx * sx + 2 * sz * tx * sx) + 2 * sy * (2 * ty * tx * sx + sy * tx * tx) + 2 * ty * (ty * sx * sx + 2 * sy * tx * sx) + 6 * sy * sy * ty * ty + 6 * sx * sx * tx * tx;
                                    Real usttt = 0.4e1 * sz * pow(tz, 0.3e1) + 0.2e1 * sz * tz * ty * ty + 0.2e1 * tz * (0.2e1 * tz * ty * sy + sz * ty * ty) + 0.2e1 * sz * tz * tx * tx + 0.2e1 * tz * (0.2e1 * tz * tx * sx + sz * tx * tx) + 0.2e1 * sy * ty * tx * tx + 0.2e1 * ty * (0.2e1 * ty * tx * sx + sy * tx * tx) + 0.4e1 * sy * pow(ty, 0.3e1) + 0.4e1 * sx * pow(tx, 0.3e1);
                                    Real urrr = 6 * rx * rx * rxx + 6 * rz * rz * rzz + 2 * rz * (2 * ryz * ry + rz * ryy) + 2 * rzz * ry * ry + 4 * rz * ryz * ry + 2 * rz * (2 * rxz * rx + rz * rxx) + 2 * rzz * rx * rx + 4 * rz * rxz * rx + 6 * ry * ry * ryy + 2 * ry * (2 * rxy * rx + ry * rxx) + 2 * ryy * rx * rx + 4 * ry * rxy * rx;
                                    Real urrs = 7 * sx * rx * rxx + rz * (3 * rz * szz + 3 * sz * rzz) + rz * (2 * rz * szz + 2 * sz * rzz) + szz * rz * rz + 7 * sz * rz * rzz + 2 * rz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * sz * (2 * ryz * ry + rz * ryy) + 2 * rz * (2 * ry * syz + 2 * ryz * sy) + 2 * szz * ry * ry + 4 * rzz * sy * ry + 4 * sz * ryz * ry + 2 * rz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * sz * (2 * rxz * rx + rz * rxx) + 2 * rz * (2 * rx * sxz + 2 * rxz * sx) + 2 * szz * rx * rx + 4 * rzz * sx * rx + 4 * sz * rxz * rx + ry * (3 * syy * ry + 3 * sy * ryy) + ry * (2 * syy * ry + 2 * sy * ryy) + syy * ry * ry + 7 * sy * ry * ryy + 4 * ryy * sx * rx + 4 * sy * rxy * rx + 2 * ry * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * sy * (2 * rxy * rx + ry * rxx) + 2 * ry * (2 * sxy * rx + 2 * sx * rxy) + 2 * syy * rx * rx + rx * (2 * sxx * rx + 2 * sx * rxx) + sxx * rx * rx + rx * (3 * sxx * rx + 3 * sx * rxx);
                                    Real urss = 7 * rx * sx * sxx + sz * (3 * rz * szz + 3 * sz * rzz) + rzz * sz * sz + sz * (2 * rz * szz + 2 * sz * rzz) + 7 * rz * sz * szz + 2 * sz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * rzz * sy * sy + 2 * sz * (2 * ry * syz + 2 * ryz * sy) + 2 * rz * (2 * syz * sy + sz * syy) + 4 * rz * syz * sy + 4 * szz * sy * ry + 2 * rz * (2 * sxz * sx + sz * sxx) + 2 * sz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * rzz * sx * sx + 2 * sz * (2 * rx * sxz + 2 * rxz * sx) + 4 * rz * sxz * sx + 4 * szz * sx * rx + sy * (3 * syy * ry + 3 * sy * ryy) + ryy * sy * sy + sy * (2 * syy * ry + 2 * sy * ryy) + 7 * ry * sy * syy + 4 * ry * sx * sxy + 4 * syy * sx * rx + 2 * ry * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * ryy * sx * sx + 2 * sy * (2 * sxy * rx + 2 * sx * rxy) + sx * (3 * sxx * rx + 3 * sx * rxx) + rxx * sx * sx + sx * (2 * sxx * rx + 2 * sx * rxx);
                                    Real usss = 6 * sx * sx * sxx + 6 * sz * sz * szz + 2 * sz * (2 * syz * sy + sz * syy) + 2 * szz * sy * sy + 4 * sz * syz * sy + 2 * sz * (2 * sxz * sx + sz * sxx) + 2 * szz * sx * sx + 4 * sz * sxz * sx + 6 * sy * sy * syy + 2 * sy * (2 * sx * sxy + sxx * sy) + 2 * syy * sx * sx + 4 * sy * sx * sxy;
                                    Real urrt = rz * (3 * tzz * rz + 3 * tz * rzz) + rz * (2 * tzz * rz + 2 * tz * rzz) + tzz * rz * rz + 7 * tz * rz * rzz + 2 * rz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * tz * (2 * ryz * ry + rz * ryy) + 2 * rz * (2 * tyz * ry + 2 * ty * ryz) + 2 * tzz * ry * ry + 4 * rzz * ty * ry + 4 * tz * ryz * ry + 2 * rz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * tz * (2 * rxz * rx + rz * rxx) + 2 * rz * (2 * txz * rx + 2 * tx * rxz) + 2 * tzz * rx * rx + 4 * rzz * tx * rx + 4 * tz * rxz * rx + ry * (2 * tyy * ry + 2 * ty * ryy) + tyy * ry * ry + ry * (3 * tyy * ry + 3 * ty * ryy) + 2 * ry * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ty * (2 * rxy * rx + ry * rxx) + 2 * ry * (2 * txy * rx + 2 * tx * rxy) + 2 * tyy * rx * rx + 7 * ty * ry * ryy + 4 * ryy * tx * rx + 4 * ty * rxy * rx + rx * (3 * txx * rx + 3 * tx * rxx) + rx * (2 * txx * rx + 2 * tx * rxx) + txx * rx * rx + 7 * tx * rx * rxx;
                                    Real urtt = tz * (3 * tzz * rz + 3 * tz * rzz) + rzz * tz * tz + tz * (2 * tzz * rz + 2 * tz * rzz) + 7 * rz * tz * tzz + 2 * rz * (2 * ty * tyz + tyy * tz) + 2 * tz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * rzz * ty * ty + 2 * tz * (2 * tyz * ry + 2 * ty * ryz) + 4 * rz * ty * tyz + 4 * tzz * ty * ry + 2 * rz * (2 * tx * txz + txx * tz) + 2 * tz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * rzz * tx * tx + 2 * tz * (2 * txz * rx + 2 * tx * rxz) + 4 * rz * tx * txz + 4 * tzz * tx * rx + ty * (3 * tyy * ry + 3 * ty * ryy) + ryy * ty * ty + ty * (2 * tyy * ry + 2 * ty * ryy) + 2 * ry * (2 * tx * txy + txx * ty) + 2 * ty * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ryy * tx * tx + 2 * ty * (2 * txy * rx + 2 * tx * rxy) + 7 * ry * ty * tyy + 4 * ry * tx * txy + 4 * tyy * tx * rx + tx * (3 * txx * rx + 3 * tx * rxx) + rxx * tx * tx + tx * (2 * txx * rx + 2 * tx * rxx) + 7 * rx * tx * txx;
                                    Real uttt = 6 * tz * tz * tzz + 2 * tz * (2 * ty * tyz + tyy * tz) + 2 * tzz * ty * ty + 4 * tz * ty * tyz + 2 * tz * (2 * tx * txz + txx * tz) + 2 * tzz * tx * tx + 4 * tz * tx * txz + 2 * ty * (2 * tx * txy + txx * ty) + 2 * tyy * tx * tx + 4 * ty * tx * txy + 6 * ty * ty * tyy + 6 * tx * tx * txx;
                                    Real usst = sz * (3 * tzz * sz + 3 * tz * szz) + sz * (2 * tzz * sz + 2 * tz * szz) + tzz * sz * sz + 7 * tz * sz * szz + 2 * sz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * tz * (2 * syz * sy + sz * syy) + 2 * sz * (2 * tyz * sy + 2 * ty * syz) + 2 * tzz * sy * sy + 4 * szz * ty * sy + 4 * tz * syz * sy + 2 * sz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * tz * (2 * sxz * sx + sz * sxx) + 2 * sz * (2 * txz * sx + 2 * tx * sxz) + 2 * tzz * sx * sx + 4 * szz * tx * sx + 4 * tz * sxz * sx + sy * (2 * tyy * sy + 2 * ty * syy) + tyy * sy * sy + sy * (3 * tyy * sy + 3 * ty * syy) + 2 * sy * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * ty * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * txy * sx + 2 * tx * sxy) + 2 * tyy * sx * sx + 7 * ty * sy * syy + 4 * syy * tx * sx + 4 * ty * sx * sxy + sx * (3 * txx * sx + 3 * tx * sxx) + sx * (2 * txx * sx + 2 * tx * sxx) + txx * sx * sx + 7 * tx * sx * sxx;
                                    Real ustt = tz * (2 * tzz * sz + 2 * tz * szz) + tz * (3 * tzz * sz + 3 * tz * szz) + szz * tz * tz + 7 * sz * tz * tzz + 2 * sz * (2 * ty * tyz + tyy * tz) + 2 * tz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * szz * ty * ty + 2 * tz * (2 * tyz * sy + 2 * ty * syz) + 4 * sz * ty * tyz + 4 * tzz * ty * sy + 2 * sz * (2 * tx * txz + txx * tz) + 2 * tz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * szz * tx * tx + 2 * tz * (2 * txz * sx + 2 * tx * sxz) + 4 * sz * tx * txz + 4 * tzz * tx * sx + ty * (3 * tyy * sy + 3 * ty * syy) + syy * ty * ty + ty * (2 * tyy * sy + 2 * ty * syy) + 2 * sy * (2 * tx * txy + txx * ty) + 2 * ty * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * syy * tx * tx + 2 * ty * (2 * txy * sx + 2 * tx * sxy) + 7 * sy * ty * tyy + 4 * sy * tx * txy + 4 * tyy * tx * sx + tx * (3 * txx * sx + 3 * tx * sxx) + sxx * tx * tx + tx * (2 * txx * sx + 2 * tx * sxx) + 7 * sx * tx * txx;
                                    Real urr = 4 * rx * rxxx + 4 * rxyy * rx + 4 * rxzz * rx + 3 * rxx * rxx + 2 * ryy * rxx + 2 * rzz * rxx + 4 * ry * rxxy + 4 * rz * rxxz + 4 * rxy * rxy + 4 * rxz * rxz + 4 * ry * ryyy + 4 * ryzz * ry + 3 * ryy * ryy + 2 * rzz * ryy + 4 * rz * ryyz + 4 * ryz * ryz + 4 * rz * rzzz + 3 * rzz * rzz;
                                    Real urs = 4 * rz * szzz + 4 * sz * rzzz + 6 * rzz * szz + 4 * rz * syyz + 4 * sz * ryyz + 2 * rzz * syy + 2 * szz * ryy + 4 * ryzz * sy + 8 * ryz * syz + 4 * ry * syzz + 4 * rz * sxxz + 4 * sz * rxxz + 2 * rzz * sxx + 2 * szz * rxx + 4 * rxzz * sx + 8 * rxz * sxz + 4 * rx * sxzz + 4 * sy * ryyy + 6 * syy * ryy + 4 * syyy * ry + 4 * ry * sxxy + 4 * sy * rxxy + 2 * ryy * sxx + 2 * syy * rxx + 4 * sxyy * rx + 8 * sxy * rxy + 4 * sx * rxyy + 4 * sx * rxxx + 6 * sxx * rxx + 4 * sxxx * rx;
                                    Real uss = 4 * sx * sxxx + 4 * sx * sxyy + 4 * sxzz * sx + 3 * sxx * sxx + 2 * sxx * syy + 2 * szz * sxx + 4 * sxxy * sy + 4 * sz * sxxz + 4 * sxy * sxy + 4 * sxz * sxz + 4 * sy * syyy + 4 * syzz * sy + 3 * syy * syy + 2 * szz * syy + 4 * sz * syyz + 4 * syz * syz + 4 * sz * szzz + 3 * szz * szz;
                                    Real urt = 4 * tz * rzzz + 6 * tzz * rzz + 4 * tzzz * rz + 4 * rz * tyyz + 4 * tz * ryyz + 2 * rzz * tyy + 2 * tzz * ryy + 4 * tyzz * ry + 8 * tyz * ryz + 4 * ty * ryzz + 4 * rz * txxz + 4 * tz * rxxz + 2 * rzz * txx + 2 * tzz * rxx + 4 * txzz * rx + 8 * txz * rxz + 4 * tx * rxzz + 4 * ty * ryyy + 6 * tyy * ryy + 4 * tyyy * ry + 2 * tyy * rxx + 4 * txyy * rx + 8 * txy * rxy + 4 * tx * rxyy + 4 * ry * txxy + 4 * ty * rxxy + 2 * ryy * txx + 4 * tx * rxxx + 6 * txx * rxx + 4 * txxx * rx;
                                    Real utt = 4 * tx * txxx + 4 * tx * txyy + 4 * tx * txzz + 3 * txx * txx + 2 * txx * tyy + 2 * txx * tzz + 4 * txxy * ty + 4 * txxz * tz + 4 * txy * txy + 4 * txz * txz + 4 * ty * tyyy + 4 * ty * tyzz + 3 * tyy * tyy + 2 * tyy * tzz + 4 * tyyz * tz + 4 * tyz * tyz + 4 * tz * tzzz + 3 * tzz * tzz;
                                    Real ust = 4 * tz * szzz + 6 * tzz * szz + 4 * tzzz * sz + 4 * sz * tyyz + 4 * tz * syyz + 2 * szz * tyy + 2 * tzz * syy + 4 * tyzz * sy + 8 * tyz * syz + 4 * ty * syzz + 4 * sz * txxz + 4 * tz * sxxz + 2 * szz * txx + 2 * tzz * sxx + 4 * txzz * sx + 8 * txz * sxz + 4 * tx * sxzz + 4 * ty * syyy + 6 * tyy * syy + 4 * tyyy * sy + 4 * sy * txxy + 4 * ty * sxxy + 2 * syy * txx + 2 * tyy * sxx + 4 * txyy * sx + 8 * txy * sxy + 4 * tx * sxyy + 4 * tx * sxxx + 6 * txx * sxx + 4 * txxx * sx;
                                    
                                    Real ur = rxxxx + 2 * rxxyy + 2 * rxxzz + ryyyy + 2 * ryyzz + rzzzz;
                                    Real us = sxxxx + 2 * sxxyy + 2 * sxxzz + syyyy + 2 * syyzz + szzzz;
                                    Real ut = txxxx + 2 * txxyy + 2 * txxzz + tyyyy + 2 * tyyzz + tzzzz;



                  // printF(" (i1,i2)=(%d,%d) urrrr=%g, urrrs=%g, urrss=%g, ursss=%g, ussss=%g, urrr=%g urrs=%g urss=%g urr=%g urs=%g uss=%g ur=%g us=%g\n",
                  //   urrrr,urrrs,urrss,ursss,ussss,urrr,urrs,urss,urss,usss,urr,urs,uss,ur,us);

                                    ForStencil(m1,m2,m3)
                                    {
                                        int m  = M123(m1,m2,m3);   
                                        coeffLocal(m,i1,i2,i3) += 
                                                                        cLapSq*(  urrrr*rrrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrrs* rrrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrss*  rrCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
                                                                                        + ursss*   rCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
                                                                                        + ussss*   iCoeff(m1)*ssssCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrrt* rrrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
                                                                                        + urrtt*  rrCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
                                                                                        + urttt*   rCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
                                                                                        + utttt*   iCoeff(m1)*   iCoeff(m2)*ttttCoeff(m3) 
                                                                                        + ussst*   iCoeff(m1)* sssCoeff(m2)*   tCoeff(m3) 
                                                                                        + usstt*   iCoeff(m1)*  ssCoeff(m2)*  ttCoeff(m3) 
                                                                                        + usttt*   iCoeff(m1)*   sCoeff(m2)* tttCoeff(m3) 

                                                                                        + urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrs *  rrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
                                                                                        + urss *   rCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
                                                                                        + usss *   iCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrt *  rrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
                                                                                        + urtt *   rCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
                                                                                        + uttt *   iCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
                                                                                        + usst *   iCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
                                                                                        + ustt *   iCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 


                                                                                        + urr  *  rrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                                                                                        + urs  *   rCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
                                                                                        + uss  *   iCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
                                                                                        + urt  *   rCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
                                                                                        + ust  *   iCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
                                                                                        + utt  *   iCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3)                                             

                                                                                        + ur   *   rCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                                                                                        + us   *   iCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
                                                                                        + ut   *   iCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
                                                                                        );
                                    }
                                }
                            }



                        }



                    } // end curvilinear


                }

            }
            else
            {
        // ----- this grid is advanced with EXPLICIT time-stepping ----
        // set the matrix the IDENTITY
                printF("+++++ IMPLICIT: grid=%d (%s) IS TREATED EXPLICITLY\n",grid,(const char*)mg.getName());        

        // set diagonal entry
                coeffLocal(mDiag,I1,I2,I3) = 1.0;

            }


            const Real extrapCoeff2[] = {1.,-2.,1.};
            const Real extrapCoeff3[] = {1.,-3.,3.,-1.};
            const Real extrapCoeff4[] = {1.,-4.,6.,-4.,1.};
            const Real extrapCoeff5[] = {1.,-5.,10.,-10.,5.,-1.};    

            if( addUpwinding )
            {
        // -------------------------------
        // --- ADD UPWIND DISSIPATION ----
        // -------------------------------

                Real upwindDissipationCoefficient = getUpwindDissipationCoefficient( grid ); 
        // Real adSosup = upwindDissipationCoefficient;
        // if( false )
        // { 
        //   Real adSosup = c*dt/( sqrt(1.*numberOfDimensions) * pow(2.,(orderOfAccuracy+1)) );

        //   printF(">>>>>>> adSosup=%12.4e, upwindDissipationCoefficient=%12.4e\n",adSosup, upwindDissipationCoefficient);
        //   upwindDissipationCoefficient = adSosup;
        // }


                OV_GET_SERIAL_ARRAY_CONDITIONAL(real,mg.inverseVertexDerivative(),rxLocal,!isRectangular);
        // macro to make the rxLocal array look 5-dimensional 
                #define DD(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2))     

                const Real uDotFactor=.5; // from D0t 
                Real adxSosup[3];
        // upwind diss coeff for Cartesian grids: 
                adxSosup[0] = uDotFactor*upwindDissipationCoefficient/dx[0];
                adxSosup[1] = uDotFactor*upwindDissipationCoefficient/dx[1];
                adxSosup[2] = uDotFactor*upwindDissipationCoefficient/dx[2]; 

                if( true )
                    printF("\nZZZ ADD UPWIND DISS TO IMPLICT MATRIX: grid=%d, dt=%14.6e, adxSosup=%14.6e, %14.6e, %14.6e ZZZ\n\n",grid, dt, adxSosup[0],adxSosup[1],adxSosup[2]);

        // getIndex(mg.gridIndexRange(),I1,I2,I3,-1);  // **** TESTING ***********

        

                if( true )
                {
          // ** NEW WAY : extrap interp neighbours as needed **
                        FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the domain
                        {
                            if( maskLocal(i1,i2,i3)>0 )
                            {
                // --- fill in the coefficients of the upwind dissipation formula ---
                // NOTE: this stencil is wider than the usual
                // NOTE: This formula must agree with the RHS computation in advWave.bf90
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
                                    if( !isRectangular && useSuperGrid==0 )
                                    {
                     // ---Upwind coefficients for a curvilinear grid ---
                     // diss-coeff ~= 1/(change in x along direction r(dir) )
                     // Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                                          if( numberOfDimensions==2 )
                                            adxSosup[dir] = upwindDissipationCoefficient*uDotFactor*sqrt( SQR(DD(i1,i2,i3,dir,0)) + SQR(DD(i1,i2,i3,dir,1)) )/dr[dir]; 
                                          else
                                            adxSosup[dir] = upwindDissipationCoefficient*uDotFactor*sqrt( SQR(DD(i1,i2,i3,dir,0)) + SQR(DD(i1,i2,i3,dir,1))  + SQR(DD(i1,i2,i3,dir,2)) )/dr[dir];                    
                                    }
                  // ------------------------------------------------
                  // --- Put the upwind coefficient in the matrix --- 
                  // ------------------------------------------------
                                    for( int iStencil=-upwindHalfStencilWidth; iStencil<=upwindHalfStencilWidth; iStencil++ )
                                    {
                                        const int i1s = i1 + iStencil*idv[0], i2s = i2 + iStencil*idv[1],  i3s = i3 + iStencil*idv[2]; // (i1s,i2s,i3s) : stencil index 
                    // Real upwStencilValue = upwindCoeff[upwCase][iStencil+upwindHalfStencilWidth];
                                        Real upwStencilValue = upwindCoeff[0][iStencil+upwindHalfStencilWidth]; // always use full stencil *new*
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
                                        coeffLocal(m,i1,i2,i3) -= adxSosup[dir] * upwStencilValue;
                                        setEquationNumber(m, e,i1,i2,i3,  cc,i1s,i2s,i3s );      // macro to set equationNumber   
                                    } // end for iStencil
                                    if( upwCase!=0 )
                                    {
                    // --- full upwind stencil involves unused points where the nask is ZERO ---
                    // Extrapolate the unused point(s) in the stencil
                                        iv[0]=i1;  iv[1]=i2; iv[2]=i3;  // centre point of the stencil 
                                        jv[0]=i1;  jv[1]=i2; jv[2]=i3;  // holds point to extrapolate
                                        isv[0]=0; isv[1]= 0; isv[2]=0;  // holds direction to extrapolate 
                                        if( upwCase==1 )
                                        {
                      // right-most point is missing   
                      //      O--O--C--O--E   E=extrapolate  
                      //   extrap4 = 1 -4 6 -4 1
                      //   extrap3 = 1 -3 3 -1 
                                            isv[dir] = -1;  jv[dir] = iv[dir] - upwindHalfStencilWidth*isv[dir];
                                            {
                                                if( false )
                                                    printF("implicit: EXTRAPOLATE upwind stencil point: grid=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d) extrap-order=%d\n",
                                                        grid,j1,j2,j3,is1,is2,is3,3);
                                                assert( maskLocal(j1,j2,j3)==0 );
                        // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                setClassify( e,j1,j2,j3, SparseRepForMGF::active );   
                                                Range all;
                                                coeffLocal(all,j1,j2,j3) = 0.0; // zero all coeff to start
                        // --- fill in the coefficients of the extrapolation formula ---
                                                for( int ie=0; ie<=3; ie++ )
                                                {
                                                    int m0 = ie + M0.getBase(); 
                                                    coeffLocal(m0,j1,j2,j3) = extrapCoeff3[ie];
                          // printF("  -- Fill m0=%d with value=%12.4e\n",m0,coeffLocal(m0,j1,j2,j3));
                                                    int j1e=j1 + ie*(is1), j2e=j2 + ie*(is2), j3e=j3 + ie*(is3);  // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                    setEquationNumber(m0, e,j1,j2,j3,  cc,j1e,j2e,j3e );         // macro to set equationNumber
                                                    if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                    {
                                                        printF("implicit: ERROR: extrapolation formula invovlves an unused point -- FIX ME: this should not happen?!\n");
                                                        printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                    j1,j2,j3,is1,is2,is3,grid,ie,j1e,j2e,j3e);
                                                        OV_ABORT("error");
                                                    }
                                                }
                                            }
                                        }
                                        else if( upwCase==2 )
                                        {
                      // left-most point is missing  
                      //     E--O--C--O--O   E=extrapolate  extrap3 = 1 -3 3 -1 
                                            isv[dir] = +1;  jv[dir] = iv[dir] - upwindHalfStencilWidth*isv[dir];
                                            {
                                                if( false )
                                                    printF("implicit: EXTRAPOLATE upwind stencil point: grid=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d) extrap-order=%d\n",
                                                        grid,j1,j2,j3,is1,is2,is3,3);
                                                assert( maskLocal(j1,j2,j3)==0 );
                        // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                setClassify( e,j1,j2,j3, SparseRepForMGF::active );   
                                                Range all;
                                                coeffLocal(all,j1,j2,j3) = 0.0; // zero all coeff to start
                        // --- fill in the coefficients of the extrapolation formula ---
                                                for( int ie=0; ie<=3; ie++ )
                                                {
                                                    int m0 = ie + M0.getBase(); 
                                                    coeffLocal(m0,j1,j2,j3) = extrapCoeff3[ie];
                          // printF("  -- Fill m0=%d with value=%12.4e\n",m0,coeffLocal(m0,j1,j2,j3));
                                                    int j1e=j1 + ie*(is1), j2e=j2 + ie*(is2), j3e=j3 + ie*(is3);  // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                    setEquationNumber(m0, e,j1,j2,j3,  cc,j1e,j2e,j3e );         // macro to set equationNumber
                                                    if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                    {
                                                        printF("implicit: ERROR: extrapolation formula invovlves an unused point -- FIX ME: this should not happen?!\n");
                                                        printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                    j1,j2,j3,is1,is2,is3,grid,ie,j1e,j2e,j3e);
                                                        OV_ABORT("error");
                                                    }
                                                }
                                            }
                                        }
                                        else if( upwCase==3 )
                                        {
                      // left and right most points are missing 
                      //     E--O--C--O--E   E=extrapolate  extrap2 = 1 -2 1 
                                            isv[dir] = -1;  jv[dir] = iv[dir] - upwindHalfStencilWidth*isv[dir];
                                            {
                                                if( false )
                                                    printF("implicit: EXTRAPOLATE upwind stencil point: grid=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d) extrap-order=%d\n",
                                                        grid,j1,j2,j3,is1,is2,is3,2);
                                                assert( maskLocal(j1,j2,j3)==0 );
                        // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                setClassify( e,j1,j2,j3, SparseRepForMGF::active );   
                                                Range all;
                                                coeffLocal(all,j1,j2,j3) = 0.0; // zero all coeff to start
                        // --- fill in the coefficients of the extrapolation formula ---
                                                for( int ie=0; ie<=2; ie++ )
                                                {
                                                    int m0 = ie + M0.getBase(); 
                                                    coeffLocal(m0,j1,j2,j3) = extrapCoeff2[ie];
                          // printF("  -- Fill m0=%d with value=%12.4e\n",m0,coeffLocal(m0,j1,j2,j3));
                                                    int j1e=j1 + ie*(is1), j2e=j2 + ie*(is2), j3e=j3 + ie*(is3);  // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                    setEquationNumber(m0, e,j1,j2,j3,  cc,j1e,j2e,j3e );         // macro to set equationNumber
                                                    if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                    {
                                                        printF("implicit: ERROR: extrapolation formula invovlves an unused point -- FIX ME: this should not happen?!\n");
                                                        printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                    j1,j2,j3,is1,is2,is3,grid,ie,j1e,j2e,j3e);
                                                        OV_ABORT("error");
                                                    }
                                                }
                                            }
                                            isv[dir] = +1;  jv[dir] = iv[dir] - upwindHalfStencilWidth*isv[dir];
                                            {
                                                if( false )
                                                    printF("implicit: EXTRAPOLATE upwind stencil point: grid=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d) extrap-order=%d\n",
                                                        grid,j1,j2,j3,is1,is2,is3,2);
                                                assert( maskLocal(j1,j2,j3)==0 );
                        // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                setClassify( e,j1,j2,j3, SparseRepForMGF::active );   
                                                Range all;
                                                coeffLocal(all,j1,j2,j3) = 0.0; // zero all coeff to start
                        // --- fill in the coefficients of the extrapolation formula ---
                                                for( int ie=0; ie<=2; ie++ )
                                                {
                                                    int m0 = ie + M0.getBase(); 
                                                    coeffLocal(m0,j1,j2,j3) = extrapCoeff2[ie];
                          // printF("  -- Fill m0=%d with value=%12.4e\n",m0,coeffLocal(m0,j1,j2,j3));
                                                    int j1e=j1 + ie*(is1), j2e=j2 + ie*(is2), j3e=j3 + ie*(is3);  // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                    setEquationNumber(m0, e,j1,j2,j3,  cc,j1e,j2e,j3e );         // macro to set equationNumber
                                                    if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                    {
                                                        printF("implicit: ERROR: extrapolation formula invovlves an unused point -- FIX ME: this should not happen?!\n");
                                                        printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                    j1,j2,j3,is1,is2,is3,grid,ie,j1e,j2e,j3e);
                                                        OV_ABORT("error");
                                                    }
                                                }
                                            }
                                        }
                                    } // end if upwCase != 0
                                }
                            } // end if mask 
                        }

                    if( false )
                        ::display(classify,"classify");
                }
                else
                {
          // ** OLD WAY **
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
                                    if( !isRectangular && useSuperGrid==0 )
                                    {
                     // ---Upwind coefficients for a curvilinear grid ---
                     // diss-coeff ~= 1/(change in x along direction r(dir) )
                     // Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                                          if( numberOfDimensions==2 )
                                            adxSosup[dir] = upwindDissipationCoefficient*uDotFactor*sqrt( SQR(DD(i1,i2,i3,dir,0)) + SQR(DD(i1,i2,i3,dir,1)) )/dr[dir]; 
                                          else
                                            adxSosup[dir] = upwindDissipationCoefficient*uDotFactor*sqrt( SQR(DD(i1,i2,i3,dir,0)) + SQR(DD(i1,i2,i3,dir,1))  + SQR(DD(i1,i2,i3,dir,2)) )/dr[dir];                    
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


                if( false )
                {
                    int grid=0;
                    displayCoeff(impCoeff[grid],sPrintF("AFTER FILL UPWINDING: implicit time-stepping matrix on grid=%d",grid));
                } 

            }



      // --- FILL BOUNDARY CONDITIONS ----

            const int extrapOrder = orderOfAccuracy+1;

      // const Real extrapCoeff3[] = {1.,-3.,3.,-1.};
      // const Real extrapCoeff4[] = {1.,-4.,6.,-4.,1.};
      // const Real extrapCoeff5[] = {1.,-5.,10.,-10.,5.,-1.};
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

            const int e=0, cc=0; // equation number and component number 
            ForBoundary(side,axis)
            {

                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);

        // Set the index-shift for this side
                is1=is2=is3=0;
                isv[axis]=1-2*side;   // +1 on left and -1 on right  

                const int axisp1 = (axis+1) % numberOfDimensions;    

                if( mg.boundaryCondition(side,axis)==dirichlet         ||
            // mg.boundaryCondition(side,axis)==CgWave::absorbing ||  // ** DO THIS FOR NOW : absorbing terminated with Dirichlet
                        mg.boundaryCondition(side,axis)==exactBC )
                {
          // ------------ FILL DIRICHLET BC ------------

                    printF("+++++ IMPLICIT BC: FILL MATRIX BC=%d FOR (grid,side,axis)=(%d,%d,%d) DIRICHLET/ABSORBING/EXACT\n",mg.boundaryCondition(side,axis),grid,side,axis);
                    FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3)
                    {
                        if( maskLocal(i1,i2,i3)!=0 ) // *wdh* avoid changing an "active" point 
                        {
                            coeffLocal(    M,i1,i2,i3) = 0.0;  // zero out any existing equations
                            coeffLocal(mDiag,i1,i2,i3) = 1.0;
                        }
                    }
          // coeffLocal(    M,Ib1,Ib2,Ib3) = 0.0;  // zero out any existing equations
          // coeffLocal(mDiag,Ib1,Ib2,Ib3) = 1.0;          

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

              // Specify that the above coeff value is the coefficient of component cc at the grid point (j1,j2,j3).
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
                else if( mg.boundaryCondition(side,axis)==CgWave::absorbing || 
                                  mg.boundaryCondition(side,axis)==CgWave::abcEM2 )
                {
                    printF("+++++ IMPLICIT BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) ABSORBING/EM2, c=%e\n",grid,side,axis,c);

          // OV_ABORT(" IMP BC ABORBING -- FINISH ME"); 

          // EM2 
          // Engquist-Majda order2 scheme
          //  D+t (-D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : left 
          //  D+t ( D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : right                
          //
          //   - D0x + c D+xD-x + .5*c D+yD-y  : left 
          //     D0x + c D+xD-x + .5*c D+yD-y  : right

          // FROM cg/mx/src/abc.bf: 
          // 
          // ! Here are 2nd-order in time approximations -- centered in space-time, solve for ghost at new time: 
          // !   D+t D0x ( u^n ) = A+t[ c1abcem2 * D+xD-x u^n + c2abcem2 D+yD-y u^n ] + f(t^n+dt/2)
          // !   Average in time operator:  A+t u^n = .5*( u^(n+1) + u^n )
          // #defineMacro ABCEM2X(i1,i2,i3,cc) ( (un(i1+is1,i2+is2,i3+is3,cc) //                     - (u(i1+is1,i2+is2,i3+is3,cc)-u(i1-is1,i2-is2,i3-is3,cc))//          - (dxa*dt)*( c1abcem2*uxx22r(i1,i2,i3,cc) + c2abcem2*uyy22r(i1,i2,i3,cc) //                      +c1abcem2*(                    -2.*un(i1,i2,i3,cc)+un(i1+is1,i2,i3,cc))/dxa**2   //                      +c2abcem2*( un(i1  ,i2-1,i3,cc)-2.*un(i1,i2,i3,cc)+un(i1  ,i2+1,i3,cc))/dx(1)**2 //                      + 2.*forcex(cc) )//                               )/(1.+c1abcem2*dt/dxa) )

          // #defineMacro ABCEM2Y(i1,i2,i3,cc) ( (un(i1+is1,i2+is2,i3+is3,cc) //                     - (u(i1+is1,i2+is2,i3+is3,cc)-u(i1-is1,i2-is2,i3-is3,cc))//          - (dya*dt)*( c1abcem2*uyy22r(i1,i2,i3,cc) + c2abcem2*uxx22r(i1,i2,i3,cc) //                      +c1abcem2*(                    -2.*un(i1,i2,i3,cc)+un(i1,i2+is2,i3,cc))/dya**2  //                      +c2abcem2*( un(i1-1,i2  ,i3,cc)-2.*un(i1,i2,i3,cc)+un(i1+1,i2  ,i3,cc))/dx(0)**2//                      + 2.*forcey(cc) )//                                   )/(1.+c1abcem2*dt/dya) )          

                    Real ca = c;
                    if( solveHelmholtz )
                    {
            // Adjust c for the EM2 absorbing BC to account for time-discretization errors
            //   D+t (Dx ) w + A+( ... )            
                        ca = c*tan(frequencyArray(0)*dt/2)/(frequencyArraySave(0)*dt/2); 
                    }

          // res = -is*(unx-ucx)/dt + (.5*c)*( unxx + ucxx) + (.25*c)*( unyy + ucyy );

                    Range Rw(-halfWidth1,halfWidth1);
                    RealArray abcCoeff(Rw,Rw,Rw);
                    abcCoeff=0.;

                    if( isRectangular )
                    {
                        const Real dxa2 = dx[axis]*dx[axis];
                        const Real dya2 = dx[axisp1]*dx[axisp1];
                        abcCoeff(-is1,-is2,0) = .5*(    ca/dxa2           )  + 1./(2.*dx[axis]*dt);    // ghost
                        abcCoeff(   0,   0,0) = .5*(-2.*ca/dxa2 - ca/dya2 );                           // boundary 
                        abcCoeff(+is1,+is2,0) = .5*(    ca/dxa2           )  - 1./(2.*dx[axis]*dt);    // first line in
                        
                        if( axis==0 )
                        {
                            abcCoeff(0,-1,0) = .5*( .5*ca/(dx[1]*dx[1]) ); 
                            abcCoeff(0,+1,0) = .5*( .5*ca/(dx[1]*dx[1]) ); 
                        }
                        else
                        {
                            abcCoeff(-1,0,0) = .5*( .5*ca/(dx[0]*dx[0]) ); 
                            abcCoeff(+1,0,0) = .5*( .5*ca/(dx[0]*dx[0]) );                   
                        }     
                    }
                    else
                    {
                          OV_ABORT(" IMP BC ABORBING -- FINISH ME: curvilinear"); 
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

                            coeffLocal(m,i1m,i2m,i3m) = abcCoeff(m1,m2,m3);

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

            if( (addUpwinding && debug>3)  )
            {
        // ::display(impCoeff[grid],sPrintF("implicit time-stepping matrix on grid=%d",grid));
                displayCoeff(impCoeff[grid],sPrintF("AFTER FILL BCS: implicit time-stepping matrix on grid=%d",grid));
                OV_ABORT("stop here for now");
            }


        }


        if( false )
        {
            int grid=0;
            displayCoeff(impCoeff[grid],sPrintF("BEFORE FINISH BC: implicit time-stepping matrix on grid=%d",grid));
        } 

        impCoeff.finishBoundaryConditions(); 

        if( false )
        {
            int grid=0;
            displayCoeff(impCoeff[grid],sPrintF("AFTER FINISH BC: implicit time-stepping matrix on grid=%d",grid));
            OV_ABORT("stop here for now");
        }    

        impSolver.setCoefficientArray( impCoeff );   // supply coefficients to Oges

    }
    timing(timeForInitialize) += getCPU()-cpu0;

  
    return 0;
}
