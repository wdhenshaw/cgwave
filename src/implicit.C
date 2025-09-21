// This file automatically generated from implicit.bC with bpp.
#include "CgWave.h"
#include "display.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "Oges.h"
#include "CompositeGridOperators.h"
#include "SparseRep.h"

// forward declaration
int coefficientsByDelta( CompositeGrid & cg, realArray & coeff, int grid, Index Rv[3], int coeffOption =4 );


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

    Real & timeForOgesSolve              = dbase.get<Real>("timeForOgesSolve");

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

    const int numberOfDimensions = cg.numberOfDimensions();

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

      // -- this is not really correct for Krylov 
            if(  solveHelmholtz && numberOfFrequencies==1 && !computeEigenmodes )
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
        // ========== Initial Guess for iterative solvers : linear extrapolation in time  ========

        // unLocal = ucLocal;           
                unLocal = 2.*ucLocal - upLocal; 
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

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];

  // ---- Fill boundary conditions for the implicit system ---
    if( true )
    {
        bool applyExplicitBoundaryConditions = false; // what should this be ? 
        bool fillImplicitBoundaryConditions  = true;
    // *fixed* Oct 1, 2024: need to pass uc for EM radiation BC
        applyBoundaryConditions( rhs, uc, t, applyExplicitBoundaryConditions, fillImplicitBoundaryConditions ); 

    // applyBoundaryConditions( rhs, rhs, t, applyExplicitBoundaryConditions, fillImplicitBoundaryConditions );
    }


    if( false || debug & 16 )
    {
        rhs.display("RHS AFTER filling BCs","%5.3f ");
    }

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

        Real cpuStart=getCPU();
        impSolver.solve( un,rhs ); 
        timeForOgesSolve += getCPU()-cpuStart; 

        int numIterations = impSolver.getNumberOfIterations();
        totalImplicitIterations += numIterations;
        totalImplicitSolves++;

        if( 1==0 || (debug & 2 && (t <= 10.*dt) ) )
            printF("CgWave::takeImplicitStep: implicit solve: t=%9.2e max-res= %8.2e (iterations=%i) ***\n",
                                  t,impSolver.getMaximumResidual(),numIterations);
    }
    else
    {
        if( false )
        {
            rhs.display("RHS before implicit solve");
        }
        Real cpuStart=getCPU();

        impSolver.solve( un,un ); 

        timeForOgesSolve += getCPU()-cpuStart; 

        if( false )
            ::display(un[0],"un after imp solve","%5.3f ");
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
    // ---- CHECK RESIDUALS IN RADIATION BOUNDARY CONDITIONS ---

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
            if( solveHelmholtz )
            {
        // Adjust c for the EM2 absorbing BC to account for time-discretization errors
        //   D+t (Dx ) w + A+( ... )            
                ca = c*tan(frequencyArray(0)*dt/2)/(frequencyArraySave(0)*dt/2); 
            }      

            ForBoundary(side,axis)
            {
                int is = 1-2*side;
                if( mg.boundaryCondition(side,axis)==absorbing || mg.boundaryCondition(side,axis)==abcEM2 )
                {
                    getBoundaryIndex(mg.gridIndexRange(),side,axis,I1,I2,I3);
          // Engquist-Majda order2 scheme
          //  D+t (-D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : left 
          //  D+t ( D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : right 
  
                    RealArray unx(I1,I2,I3), unxx(I1,I2,I3), unyy(I1,I2,I3);
                    RealArray ucx(I1,I2,I3), ucxx(I1,I2,I3), ucyy(I1,I2,I3);
                    if( orderOfAccuracy==2 )
                    {
                        unxx(I1,I2,I3) = (unLocal(I1+1,I2,I3)-2.*unLocal(I1,I2,I3)+unLocal(I1-1,I2,I3))/(dx[0]*dx[0]);
                        unyy(I1,I2,I3) = (unLocal(I1,I2+1,I3)-2.*unLocal(I1,I2,I3)+unLocal(I1,I2-1,I3))/(dx[1]*dx[1]);

                        ucxx(I1,I2,I3) = (ucLocal(I1+1,I2,I3)-2.*ucLocal(I1,I2,I3)+ucLocal(I1-1,I2,I3))/(dx[0]*dx[0]);
                        ucyy(I1,I2,I3) = (ucLocal(I1,I2+1,I3)-2.*ucLocal(I1,I2,I3)+ucLocal(I1,I2-1,I3))/(dx[1]*dx[1]); 
                    }
                    else if( orderOfAccuracy==4 )
                    {
                        unxx(I1,I2,I3) = (-unLocal(I1-2,I2,I3)+16.*unLocal(I1-1,I2,I3)-30.*unLocal(I1,I2,I3)+16.*unLocal(I1+1,I2,I3)-unLocal(I1+2,I2,I3))/(12.*dx[0]*dx[0]);
                        unyy(I1,I2,I3) = (-unLocal(I1,I2-2,I3)+16.*unLocal(I1,I2-1,I3)-30.*unLocal(I1,I2,I3)+16.*unLocal(I1,I2+1,I3)-unLocal(I1,I2+2,I3))/(12.*dx[1]*dx[1]);
                        ucxx(I1,I2,I3) = (-ucLocal(I1-2,I2,I3)+16.*ucLocal(I1-1,I2,I3)-30.*ucLocal(I1,I2,I3)+16.*ucLocal(I1+1,I2,I3)-ucLocal(I1+2,I2,I3))/(12.*dx[0]*dx[0]);
                        ucyy(I1,I2,I3) = (-ucLocal(I1,I2-2,I3)+16.*ucLocal(I1,I2-1,I3)-30.*ucLocal(I1,I2,I3)+16.*ucLocal(I1,I2+1,I3)-ucLocal(I1,I2+2,I3))/(12.*dx[1]*dx[1]);             
                    }
                    else
                    {
                        OV_ABORT("finish me");
                    }
                      
                    RealArray res(I1,I2,I3);
  
                    if( axis==0 ) 
                    {       
                        if( orderOfAccuracy==2 )
                        {      
                            unx(I1,I2,I3) = (unLocal(I1+1,I2,I3) -unLocal(I1-1,I2,I3))/(2.*dx[0]);
                            ucx(I1,I2,I3) = (ucLocal(I1+1,I2,I3) -ucLocal(I1-1,I2,I3))/(2.*dx[0]);
                        }
                        else
                        {
                            unx(I1,I2,I3) = (unLocal(I1-2,I2,I3)-8.*unLocal(I1-1,I2,I3)+8.*unLocal(I1+1,I2,I3)-unLocal(I1+2,I2,I3))/(12.*dx[0]);
                            ucx(I1,I2,I3) = (ucLocal(I1-2,I2,I3)-8.*ucLocal(I1-1,I2,I3)+8.*ucLocal(I1+1,I2,I3)-ucLocal(I1+2,I2,I3))/(12.*dx[0]);
                        }
                        res = -is*(unx-ucx)/dt  + (.5*ca)*( unxx + ucxx ) + (.25*ca)*( unyy + ucyy );
                    }
                    else
                    {
                        if( orderOfAccuracy==2 )
                        {
                            unx(I1,I2,I3) = (unLocal(I1,I2+1,I3) -unLocal(I1,I2-1,I3))/(2.*dx[1]);
                            ucx(I1,I2,I3) = (ucLocal(I1,I2+1,I3) -ucLocal(I1,I2-1,I3))/(2.*dx[1]);
                        }
                        else
                        {
                            unx(I1,I2,I3) = (unLocal(I1,I2-2,I3)-8.*unLocal(I1,I2-1,I3)+8.*unLocal(I1,I2+1,I3)-unLocal(I1,I2+2,I3))/(12.*dx[1]);
                            ucx(I1,I2,I3) = (ucLocal(I1,I2-2,I3)-8.*ucLocal(I1,I2-1,I3)+8.*ucLocal(I1,I2+1,I3)-ucLocal(I1,I2+2,I3))/(12.*dx[1]);              
                        }
                        res = -is*(unx-ucx)/dt  + (.5*ca)*( unyy + ucyy ) + (.25*ca)*( unxx + ucxx );
                    }

                    RealArray res2(I1,I2,I3);
                    if( orderOfAccuracy==4 )
                    {
             // check the CBC
                        RealArray unxxxx(I1,I2,I3), ucxxxx(I1,I2,I3);
                        RealArray unxxyy(I1,I2,I3), ucxxyy(I1,I2,I3);
                        RealArray unyyyy(I1,I2,I3), ucyyyy(I1,I2,I3);

                        unxxxx(I1,I2,I3) = (unLocal(I1-2,I2,I3)-4.*unLocal(I1-1,I2,I3)+6.*unLocal(I1,I2,I3)-4.*unLocal(I1+1,I2,I3)+unLocal(I1+2,I2,I3))/(dx[0]*dx[0]*dx[0]*dx[0]);
                        unyyyy(I1,I2,I3) = (unLocal(I1,I2-2,I3)-4.*unLocal(I1,I2-1,I3)+6.*unLocal(I1,I2,I3)-4.*unLocal(I1,I2+1,I3)+unLocal(I1,I2+2,I3))/(dx[1]*dx[1]*dx[1]*dx[1]);

                        unxxyy(I1,I2,I3) = (  4.*unLocal(I1,I2,I3)
                                                                  -2.*unLocal(I1-1,I2  ,I3)-2.*unLocal(I1+1,I2  ,I3) -2.*unLocal(I1  ,I2-1,I3)-2.*unLocal(I1  ,I2+1,I3)
                                                                        +unLocal(I1-1,I2-1,I3)   +unLocal(I1+1,I2-1,I3)    +unLocal(I1-1,I2+1,I3)   +unLocal(I1+1,I2+1,I3)
                                                              )/(dx[0]*dx[0]*dx[1]*dx[1]);
                        
                        ucxxxx(I1,I2,I3) = (ucLocal(I1-2,I2,I3)-4.*ucLocal(I1-1,I2,I3)+6.*ucLocal(I1,I2,I3)-4.*ucLocal(I1+1,I2,I3)+ucLocal(I1+2,I2,I3))/(dx[0]*dx[0]*dx[0]*dx[0]);
                        ucyyyy(I1,I2,I3) = (ucLocal(I1,I2-2,I3)-4.*ucLocal(I1,I2-1,I3)+6.*ucLocal(I1,I2,I3)-4.*ucLocal(I1,I2+1,I3)+ucLocal(I1,I2+2,I3))/(dx[1]*dx[1]*dx[1]*dx[1]);

                        ucxxyy(I1,I2,I3) = (  4.*ucLocal(I1,I2,I3)
                                                                  -2.*ucLocal(I1-1,I2,I3)-2.*ucLocal(I1+1,I2,I3) -2.*ucLocal(I1,I2-1,I3)-2.*ucLocal(I1,I2+1,I3)
                                                                  +ucLocal(I1-1,I2-1,I3)+ucLocal(I1+1,I2-1,I3)+ucLocal(I1-1,I2+1,I3)+ucLocal(I1+1,I2+1,I3)
                                                              )/(dx[0]*dx[0]*dx[1]*dx[1]);

                        if( axis==0 ) 
                        {       
                            RealArray unxxx(I1,I2,I3), ucxxx(I1,I2,I3);
                            unxxx(I1,I2,I3) = (-unLocal(I1-2,I2,I3)+2.*unLocal(I1-1,I2,I3)-2.*unLocal(I1+1,I2,I3)+unLocal(I1+2,I2,I3))/(2.*dx[0]*dx[0]*dx[0]);
                            ucxxx(I1,I2,I3) = (-ucLocal(I1-2,I2,I3)+2.*ucLocal(I1-1,I2,I3)-2.*ucLocal(I1+1,I2,I3)+ucLocal(I1+2,I2,I3))/(2.*dx[0]*dx[0]*dx[0]);
                            res2 = -is*(unxxx-ucxxx)/dt  + (.5*ca)*( unxxxx + ucxxxx ) + (.25*ca)*( unxxyy + ucxxyy );
                        }
                        else
                        {
                            RealArray unyyy(I1,I2,I3), ucyyy(I1,I2,I3);
                            unyyy(I1,I2,I3) = (-unLocal(I1,I2-2,I3)+2.*unLocal(I1,I2-1,I3)-2.*unLocal(I1,I2+1,I3)+unLocal(I1,I2+2,I3))/(2.*dx[1]*dx[1]*dx[1]);
                            ucyyy(I1,I2,I3) = (-ucLocal(I1,I2-2,I3)+2.*ucLocal(I1,I2-1,I3)-2.*ucLocal(I1,I2+1,I3)+ucLocal(I1,I2+2,I3))/(2.*dx[1]*dx[1]*dx[1]);              
                            res2 = -is*(unyyy-ucyyy)/dt  + (.5*ca)*( unyyyy + ucyyyy ) + (.25*ca)*( unxxyy + ucxxyy );
                        }
                    }

                    if( twilightZone )
                    {
                        assert( dbase.get<OGFunction*>("tz")!=NULL );
                        OGFunction & e = *dbase.get<OGFunction*>("tz"); 
                        OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
            // RealArray utx(I1,I2,I3);
                        const Real tc = t-dt;    // previous time
            // const Real th = t-.5*dt; // half time 
                        const int isRectangular=0;

                        e.gd( ucxx,xLocal,numberOfDimensions,isRectangular,0,2,0,0,I1,I2,I3,0,tc ); 
                        e.gd( ucyy,xLocal,numberOfDimensions,isRectangular,0,0,2,0,I1,I2,I3,0,tc ); 

                        e.gd( unxx,xLocal,numberOfDimensions,isRectangular,0,2,0,0,I1,I2,I3,0,t ); 
                        e.gd( unyy,xLocal,numberOfDimensions,isRectangular,0,0,2,0,I1,I2,I3,0,t ); 

                        if( axis==0 ) 
                        {
             //::display(ucx,"ucx (computed)","%9.2e ");                  
             //::display(unx,"unx (computed)","%9.2e ");                  
                          e.gd( ucx ,xLocal,numberOfDimensions,isRectangular,0,1,0,0,I1,I2,I3,0,tc );
                          e.gd( unx ,xLocal,numberOfDimensions,isRectangular,0,1,0,0,I1,I2,I3,0,t  );

              // ::display(ucx,"ucx (true)","%9.2e ");                  
              // ::display(unx,"unx (true)","%9.2e ");                  
                            res -= -is*(unx-ucx)/dt  + (.5*ca)*( unxx + ucxx ) + (.25*ca)*( unyy + ucyy );
                        }
                        else
                        {
                            e.gd( ucx ,xLocal,numberOfDimensions,isRectangular,0,0,1,0,I1,I2,I3,0,tc );
                            e.gd( unx ,xLocal,numberOfDimensions,isRectangular,0,0,1,0,I1,I2,I3,0,t  );

                            res -= -is*(unx-ucx)/dt  + (.5*ca)*( unyy + ucyy ) + (.25*ca)*( unxx + ucxx );
                        }  


                    }

                    Real maxRes = max(fabs(res));

                    if( orderOfAccuracy==4 )
                    {
                        Real maxRes2= max(fabs(res2));
                        if( 1==1 || maxRes>1.e10 || maxRes2>1.e-9 )
                        {
                            printF("AFTER IMP solve: (side,axis,grid)=(%d,%d,%d) c=%12.4e, cEM2=%12.4e, bc=EM2 maxRes=[%9.2e,%9.2e]\n",side,axis,grid,c,ca,maxRes,maxRes2);
                        }
                    }
                    else
                    {
                        if( maxRes > 1e-10 )
                        {
                            printF("WARNING : AFTER IMP solve: (side,axis,grid)=(%d,%d,%d) c=%12.4e, cEM2=%12.4e, bc=absorbing maxRes=%9.2e\n",side,axis,grid,c,ca,maxRes);
                            if( debug & 4 ) 
                                ::display(res,"res","%9.2e ");
                        }
                    }
                }
            }
        }
    // OV_ABORT("STOP here for now");
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
#define setEquationNumber(m, ni,i1,i2,i3,  nj,j1,j2,j3 )equationNumberLocal(m,i1,i2,i3)=indexToEquation( nj,j1,j2,j3)

// =======================================================================
// =======================================================================
#define setClassify(n,i1,i2,i3, type) classifyLocal(i1,i2,i3,n)=type

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

// ------------------------------------------------------
// Macro: FILL MATRIX WITH THE DIRICHLET OR EXACT BC
// ------------------------------------------------------

// ------------------------------------------------------
// Macro: FILL MATRIX WITH THE NEUMANN BC
// ------------------------------------------------------


// ============================================================================================
// Macro to set the extrapolation width
// ============================================================================================



// ============================================================================================
// Macro to add foruth order correction terms to the implicit matrix
//  
// Add modified equation term 
//        (L_2h)^2  
// ============================================================================================

// ============================================================================================
// Macro to add foruth order correction terms to the implicit matrix
//  
// Add modified equation term 
//        (L_2h)^2  
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
    const int & orderOfExtrapolation     = dbase.get<int>("orderOfExtrapolation");
    const IntegerArray & gridIsImplicit  = dbase.get<IntegerArray>("gridIsImplicit");

    const int & upwind                   = dbase.get<int>("upwind");
    const int & implicitUpwind           = dbase.get<int>("implicitUpwind");

    const bool addUpwinding = upwind && implicitUpwind;

    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");

    FILE *& debugFile  = dbase.get<FILE*>("debugFile");
    FILE *& pDebugFile = dbase.get<FILE*>("pDebugFile");  

  // addUpwinding=false; // *********** TURN OFF FOR NOW ***************

  // CHECK number of ghost points:
  //   Order=4 + Neumann needs one extra 
    if( orderOfAccuracy>=4  )
    {
        const int orderOfAccuracyInSpace = orderOfAccuracy;
        const int minGhostNeeded = orderOfAccuracyInSpace/2 +1;
        Range Rx=cg.numberOfDimensions();
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg=cg[grid];
            bool hasNeumannBC = false;
            ForBoundary(side,axis)
            {
        // Check for Neumann BC's
                if( mg.boundaryCondition(side,axis)==neumann )
                {
                    hasNeumannBC=true;
                    break;
                }

            }  
      // printF("\n **** hasNeumannBC = %d, minGhostNeeded=%d ***\n\n",hasNeumannBC,minGhostNeeded);

            if( hasNeumannBC )
            {
                const IntegerArray & numberOfGhostPoints = mg.numberOfGhostPoints();
                const int numGhost = min(numberOfGhostPoints(Range(0,1),Rx));
                if( numGhost < minGhostNeeded )
                {
                    printF("CgWave::formImplicitTimeSteppingMatrix:ERROR: the grid does not have enough ghost points.\n"
                    "   orderOfAccuracy=%i + Neumann Boundary Conditions needs at least %i ghost points.\n"
                    "   You could remake the grid with more ghost points to fix this error.\n",
                    orderOfAccuracyInSpace,minGhostNeeded);
                    OV_ABORT("ERROR");
                    
                } 
            }
        }
    }

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
        printf("CgWave::formImplicitTimeSteppingMatrix: Changing the implicit solver parameters... \n");
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
    if( orderOfExtrapolation>0 )
    {
        printF("\n $$$$$$$$$$ CgWave::formImplicitTimeSteppingMatrix:INFO: set orderOfExtrapolation=%d $$$$$$$$\n",orderOfExtrapolation);
        
        impSolver.set(OgesParameters::THEorderOfExtrapolation,orderOfExtrapolation);
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

    usePredefined = usePredefined && ( useMultigrid || bcApproach != useCompatibilityBoundaryConditions);

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

        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
        Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
        Index Jbv[3], &Jb1=Jbv[0], &Jb2=Jbv[1], &Jb3=Jbv[2];
        Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];
        int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];
        int jv[3], &j1=jv[0], &j2=jv[1], &j3=jv[2];
        int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];
        int m1,m2,m3; 

        const bool useCompatibility = orderOfAccuracy==4 && bcApproach==useCompatibilityBoundaryConditions;    

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
                OV_GET_SERIAL_ARRAY(int,equationNumber,equationNumberLocal);
                OV_GET_SERIAL_ARRAY(int,classify,classifyLocal);
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
        // ----- this grid is advanced with IMPLICIT time-stepping ----
  

                getIndex(mg.gridIndexRange(),I1,I2,I3);
                bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,I1,I2,I3);  
                if( ok )
                {

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

                }

        // if( orderOfAccuracy==4 )
                if( orderOfAccuracyInTime==4 ) // **BUG FIXED** April 24, 2025
                {
           // Add modified equation term 
           //     (L_2h)^2

                    if( debug & 1 && grid<3 )
                        printF("\n ***ADD Fourth-order in time coefficients to the Implicit Matrix ****\n\n");

                        const Real cLapSq = cImp(-1,1)*(c*dt)*(c*dt)*(c*dt)*(c*dt);  // CHECK ME   
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
              //              -(c^4*dt^2/12) Delta^2( cImp(1,1) *u^{n+1} + cImp(0,1) *u^n + cImp(-1,1)* u^{n-1}  )  :  fourth-order ceoff cImp(-1:1,1) 
              // For accuracy the weights depend on one parameter beta2 for second-order,
              // and a second parameter beta4 for fourth-order: (See notes in research/timeStepping/implicitTaylorSchemes.pdf)
                            ForStencil(m1,m2,m3)
                            {
                                int m  = M123(m1,m2,m3);   
                                coeffLocal(m,I1,I2,I3) += cLapSq*lapSq(m1,m2,m3);
                            }
                        }
                        else
                        {
              // OV_ABORT("implicit matrix: finish me for order=4 CURVLINEAR");
                            if( debug & 1 )
                                printF("implicit matrix: order=4 CURVLINEAR: ADD (Lap_2h)^2 correction term\n");
              // OV_GET_SERIAL_ARRAY(Real,lapSqCoeff[grid],lapSqCoeffLocal);  // for testing
                            OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rxLocal);
              // macro to make the rxLocal array look 5-dimensional 
                            #define DD(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 
              // See advWaveStencil.bf90
              // ! -- Coefficients in the Laplacian (scaled)
              // c200(i1,i2,i3) = (rx**2 + ry**2   )*dr1i**2
              // c110(i1,i2,i3) = 2.*(rx*sx + ry*sy)*dr1i*dr2i
              // c020(i1,i2,i3) = (sx**2 + sy**2   )*dr2i**2
              // c100(i1,i2,i3) = (rxx + ryy       )*dr1i
              // c010(i1,i2,i3) = (sxx + syy       )*dr2i 
                            if( numberOfDimensions==2 )
                            {
                                Index J1,J2,J3;
                                int extra=1;
                                getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
                                bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,J1,J2,J3);
                                if( ok )
                                {
                                    RealArray c200(J1,J2,J3), c110(J1,J2,J3), c020(J1,J2,J3), c100(J1,J2,J3), c010(J1,J2,J3);
                  // -- evaluate the SCALED coefficients of the Laplacian:
                  //    Delta = c200*Delta+r Delta-r + c110*Delta0r Delta0s + ...
                  // SEE ALSO
                  // research/compatibility/
                  //     lcbcDiff.mpl       : defines Qh 
                  //     writeLcbcFiles.mpl : 
                                    FOR_3D(i1,i2,i3,J1,J2,J3)
                                    {
                                        Real rx = DD(i1,i2,i3,0,0); 
                                        Real ry = DD(i1,i2,i3,0,1); 
                                        Real sx = DD(i1,i2,i3,1,0); 
                                        Real sy = DD(i1,i2,i3,1,1);         
                                        c200(i1,i2,i3) = ( SQR(rx) + SQR(ry) )/( SQR(dr[0]) );
                                        c020(i1,i2,i3) = ( SQR(sx) + SQR(sy) )/( SQR(dr[1]) );
                                        c110(i1,i2,i3) = 2.*( rx*sx + ry*sy )/( dr[0]*dr[1] ); // **CHECK FACTOR OF 1/4 
                                        Real rxr = (DD(i1+1,i2,i3,0,0)-DD(i1-1,i2,i3,0,0))/(2.*dr[0]);
                                        Real ryr = (DD(i1+1,i2,i3,0,1)-DD(i1-1,i2,i3,0,1))/(2.*dr[0]);
                                        Real sxr = (DD(i1+1,i2,i3,1,0)-DD(i1-1,i2,i3,1,0))/(2.*dr[0]);
                                        Real syr = (DD(i1+1,i2,i3,1,1)-DD(i1-1,i2,i3,1,1))/(2.*dr[0]);
                                        Real rxs = (DD(i1,i2+1,i3,0,0)-DD(i1,i2-1,i3,0,0))/(2.*dr[1]);
                                        Real rys = (DD(i1,i2+1,i3,0,1)-DD(i1,i2-1,i3,0,1))/(2.*dr[1]);
                                        Real sxs = (DD(i1,i2+1,i3,1,0)-DD(i1,i2-1,i3,1,0))/(2.*dr[1]);
                                        Real sys = (DD(i1,i2+1,i3,1,1)-DD(i1,i2-1,i3,1,1))/(2.*dr[1]);
                                        Real rxx= rx*rxr + sx*rxs;
                                        Real ryy= ry*ryr + sy*rys;
                                        Real sxx= rx*sxr + sx*sxs;
                                        Real syy= ry*syr + sy*sys;
                                        c100(i1,i2,i3) = (rxx + ryy)/(dr[0]); // rxx + ryy  // **CHECK FACTOR OF 1/2
                                        c010(i1,i2,i3) = (sxx + syy)/(dr[1]); // sxx + syy
                                    }
                                    Range R5(-2,2);
                                    RealArray cp(R5,R5);
                                    FOR_3(i1,i2,i3,I1,I2,I3)
                                    {
                    // coefficients of Lap^2 (from research/compatibility/writeLcbcFiles.mpl --> lcbcEquationsDirichlet2dOrder4.h)
                                        if( maskLocal(i1,i2,i3)>0 )
                                        {
                                            cp(-2,-2) =1./16.*c110(i1,i2,i3)*c110(i1-1,i2-1,i3);
                                            cp(-1,-2) =1./4.*c020(i1,i2,i3)*c110(i1,i2-1,i3)+c110(i1,i2,i3)*(1./4.*c020(i1-1,i2-1,i3)-1./8.*c010(i1-1,i2-1,i3))-1./8.*c010(i1,i2,i3)*c110(i1,i2-1,i3);
                                            cp( 0,-2) =c020(i1,i2,i3)*(c020(i1,i2-1,i3)-1./2.*c010(i1,i2-1,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1-1,i2-1,i3)-1./16.*c110(i1+1,i2-1,i3))+c010(i1,i2,i3)*(-1./2.*c020(i1,i2-1,i3)+1./4.*c010(i1,i2-1,i3));
                                            cp( 1,-2) =-1./4.*c020(i1,i2,i3)*c110(i1,i2-1,i3)+c110(i1,i2,i3)*(-1./4.*c020(i1+1,i2-1,i3)+1./8.*c010(i1+1,i2-1,i3))+1./8.*c010(i1,i2,i3)*c110(i1,i2-1,i3);
                                            cp( 2,-2) =1./16.*c110(i1,i2,i3)*c110(i1+1,i2-1,i3);
                                            cp(-2,-1) =1./4.*c200(i1,i2,i3)*c110(i1-1,i2,i3)+c110(i1,i2,i3)*(1./4.*c200(i1-1,i2-1,i3)-1./8.*c100(i1-1,i2-1,i3))-1./8.*c100(i1,i2,i3)*c110(i1-1,i2,i3);
                                            cp(-1,-1) =c200(i1,i2,i3)*(-1./2.*c010(i1-1,i2,i3)+c020(i1-1,i2,i3)-1./2.*c110(i1,i2,i3))+c020(i1,i2,i3)*(-1./2.*c100(i1,i2-1,i3)+c200(i1,i2-1,i3)-1./2.*c110(i1,i2,i3))+c110(i1,i2,i3)*(-1./2.*c020(i1-1,i2-1,i3)-1./2.*c200(i1-1,i2-1,i3))+c100(i1,i2,i3)*(1./4.*c010(i1-1,i2,i3)-1./2.*c020(i1-1,i2,i3))+c010(i1,i2,i3)*(1./4.*c100(i1,i2-1,i3)-1./2.*c200(i1,i2-1,i3));
                                            cp( 0,-1) =c200(i1,i2,i3)*(-1./4.*c110(i1-1,i2,i3)+1./4.*c110(i1+1,i2,i3)+c010(i1,i2,i3)-2*c020(i1,i2,i3))+c020(i1,i2,i3)*(-2*c020(i1,i2-1,i3)-2*c200(i1,i2-1,i3)+c010(i1,i2,i3)-2*c020(i1,i2,i3))+c110(i1,i2,i3)*(1./4.*c200(i1-1,i2-1,i3)-1./4.*c200(i1+1,i2-1,i3)+1./8.*c100(i1-1,i2-1,i3)+1./8.*c100(i1+1,i2-1,i3))+c100(i1,i2,i3)*(1./8.*c110(i1-1,i2,i3)+1./8.*c110(i1+1,i2,i3))+c010(i1,i2,i3)*(c020(i1,i2-1,i3)+c200(i1,i2-1,i3));
                                            cp( 1,-1) =c200(i1,i2,i3)*(-1./2.*c010(i1+1,i2,i3)+c020(i1+1,i2,i3)+1./2.*c110(i1,i2,i3))+c020(i1,i2,i3)*(1./2.*c100(i1,i2-1,i3)+c200(i1,i2-1,i3)+1./2.*c110(i1,i2,i3))+c110(i1,i2,i3)*(1./2.*c020(i1+1,i2-1,i3)+1./2.*c200(i1+1,i2-1,i3))+c100(i1,i2,i3)*(-1./4.*c010(i1+1,i2,i3)+1./2.*c020(i1+1,i2,i3))+c010(i1,i2,i3)*(-1./4.*c100(i1,i2-1,i3)-1./2.*c200(i1,i2-1,i3));
                                            cp( 2,-1)=-1./4.*c200(i1,i2,i3)*c110(i1+1,i2,i3)+c110(i1,i2,i3)*(-1./4.*c200(i1+1,i2-1,i3)-1./8.*c100(i1+1,i2-1,i3))-1./8.*c100(i1,i2,i3)*c110(i1+1,i2,i3);
                                            cp(-2, 0)=c200(i1,i2,i3)*(c200(i1-1,i2,i3)-1./2.*c100(i1-1,i2,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1-1,i2-1,i3)-1./16.*c110(i1-1,i2+1,i3))+c100(i1,i2,i3)*(-1./2.*c200(i1-1,i2,i3)+1./4.*c100(i1-1,i2,i3));
                                            cp(-1, 0)=c200(i1,i2,i3)*(-2*c020(i1-1,i2,i3)-2*c200(i1-1,i2,i3)+c100(i1,i2,i3)-2*c200(i1,i2,i3))+c020(i1,i2,i3)*(-1./4.*c110(i1,i2-1,i3)+1./4.*c110(i1,i2+1,i3)+c100(i1,i2,i3)-2*c200(i1,i2,i3))+c110(i1,i2,i3)*(1./4.*c020(i1-1,i2-1,i3)-1./4.*c020(i1-1,i2+1,i3)+1./8.*c010(i1-1,i2-1,i3)+1./8.*c010(i1-1,i2+1,i3))+c100(i1,i2,i3)*(c020(i1-1,i2,i3)+c200(i1-1,i2,i3))+c010(i1,i2,i3)*(1./8.*c110(i1,i2-1,i3)+1./8.*c110(i1,i2+1,i3));
                                            cp( 0, 0)=c200(i1,i2,i3)*(c200(i1-1,i2,i3)+c200(i1+1,i2,i3)+1./2.*c100(i1-1,i2,i3)-1./2.*c100(i1+1,i2,i3)+4*c020(i1,i2,i3)+4*c200(i1,i2,i3))+c020(i1,i2,i3)*(c020(i1,i2-1,i3)+c020(i1,i2+1,i3)+1./2.*c010(i1,i2-1,i3)-1./2.*c010(i1,i2+1,i3)+4*c020(i1,i2,i3)+4*c200(i1,i2,i3))+c110(i1,i2,i3)*(1./16.*c110(i1-1,i2-1,i3)+1./16.*c110(i1-1,i2+1,i3)+1./16.*c110(i1+1,i2-1,i3)+1./16.*c110(i1+1,i2+1,i3))+c100(i1,i2,i3)*(-1./2.*c200(i1-1,i2,i3)+1./2.*c200(i1+1,i2,i3)-1./4.*c100(i1-1,i2,i3)-1./4.*c100(i1+1,i2,i3))+c010(i1,i2,i3)*(-1./2.*c020(i1,i2-1,i3)+1./2.*c020(i1,i2+1,i3)-1./4.*c010(i1,i2-1,i3)-1./4.*c010(i1,i2+1,i3));
                                            cp( 1, 0)=c200(i1,i2,i3)*(-2*c020(i1+1,i2,i3)-2*c200(i1+1,i2,i3)-c100(i1,i2,i3)-2*c200(i1,i2,i3))+c020(i1,i2,i3)*(1./4.*c110(i1,i2-1,i3)-1./4.*c110(i1,i2+1,i3)-c100(i1,i2,i3)-2*c200(i1,i2,i3))+c110(i1,i2,i3)*(-1./4.*c020(i1+1,i2-1,i3)+1./4.*c020(i1+1,i2+1,i3)-1./8.*c010(i1+1,i2-1,i3)-1./8.*c010(i1+1,i2+1,i3))+c100(i1,i2,i3)*(-c020(i1+1,i2,i3)-c200(i1+1,i2,i3))+c010(i1,i2,i3)*(-1./8.*c110(i1,i2-1,i3)-1./8.*c110(i1,i2+1,i3));
                                            cp( 2, 0)=c200(i1,i2,i3)*(c200(i1+1,i2,i3)+1./2.*c100(i1+1,i2,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1+1,i2-1,i3)-1./16.*c110(i1+1,i2+1,i3))+c100(i1,i2,i3)*(1./2.*c200(i1+1,i2,i3)+1./4.*c100(i1+1,i2,i3));
                                            cp(-2, 1)=-1./4.*c200(i1,i2,i3)*c110(i1-1,i2,i3)+c110(i1,i2,i3)*(-1./4.*c200(i1-1,i2+1,i3)+1./8.*c100(i1-1,i2+1,i3))+1./8.*c100(i1,i2,i3)*c110(i1-1,i2,i3);
                                            cp(-1, 1)=c200(i1,i2,i3)*(1./2.*c010(i1-1,i2,i3)+c020(i1-1,i2,i3)+1./2.*c110(i1,i2,i3))+c020(i1,i2,i3)*(-1./2.*c100(i1,i2+1,i3)+c200(i1,i2+1,i3)+1./2.*c110(i1,i2,i3))+c110(i1,i2,i3)*(1./2.*c020(i1-1,i2+1,i3)+1./2.*c200(i1-1,i2+1,i3))+c100(i1,i2,i3)*(-1./4.*c010(i1-1,i2,i3)-1./2.*c020(i1-1,i2,i3))+c010(i1,i2,i3)*(-1./4.*c100(i1,i2+1,i3)+1./2.*c200(i1,i2+1,i3));
                                            cp( 0, 1)=c200(i1,i2,i3)*(1./4.*c110(i1-1,i2,i3)-1./4.*c110(i1+1,i2,i3)-c010(i1,i2,i3)-2*c020(i1,i2,i3))+c020(i1,i2,i3)*(-2*c020(i1,i2+1,i3)-2*c200(i1,i2+1,i3)-c010(i1,i2,i3)-2*c020(i1,i2,i3))+c110(i1,i2,i3)*(-1./4.*c200(i1-1,i2+1,i3)+1./4.*c200(i1+1,i2+1,i3)-1./8.*c100(i1-1,i2+1,i3)-1./8.*c100(i1+1,i2+1,i3))+c100(i1,i2,i3)*(-1./8.*c110(i1-1,i2,i3)-1./8.*c110(i1+1,i2,i3))+c010(i1,i2,i3)*(-c020(i1,i2+1,i3)-c200(i1,i2+1,i3));
                                            cp( 1, 1)=c200(i1,i2,i3)*(1./2.*c010(i1+1,i2,i3)+c020(i1+1,i2,i3)-1./2.*c110(i1,i2,i3))+c020(i1,i2,i3)*(1./2.*c100(i1,i2+1,i3)+c200(i1,i2+1,i3)-1./2.*c110(i1,i2,i3))+c110(i1,i2,i3)*(-1./2.*c020(i1+1,i2+1,i3)-1./2.*c200(i1+1,i2+1,i3))+c100(i1,i2,i3)*(1./4.*c010(i1+1,i2,i3)+1./2.*c020(i1+1,i2,i3))+c010(i1,i2,i3)*(1./4.*c100(i1,i2+1,i3)+1./2.*c200(i1,i2+1,i3));
                                            cp( 2, 1)=1./4.*c200(i1,i2,i3)*c110(i1+1,i2,i3)+c110(i1,i2,i3)*(1./4.*c200(i1+1,i2+1,i3)+1./8.*c100(i1+1,i2+1,i3))+1./8.*c100(i1,i2,i3)*c110(i1+1,i2,i3);
                                            cp(-2, 2)=1./16.*c110(i1,i2,i3)*c110(i1-1,i2+1,i3);
                                            cp(-1, 2)=-1./4.*c020(i1,i2,i3)*c110(i1,i2+1,i3)+c110(i1,i2,i3)*(-1./4.*c020(i1-1,i2+1,i3)-1./8.*c010(i1-1,i2+1,i3))-1./8.*c010(i1,i2,i3)*c110(i1,i2+1,i3);
                                            cp( 0, 2)=c020(i1,i2,i3)*(c020(i1,i2+1,i3)+1./2.*c010(i1,i2+1,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1-1,i2+1,i3)-1./16.*c110(i1+1,i2+1,i3))+c010(i1,i2,i3)*(1./2.*c020(i1,i2+1,i3)+1./4.*c010(i1,i2+1,i3));
                                            cp( 1, 2)=1./4.*c020(i1,i2,i3)*c110(i1,i2+1,i3)+c110(i1,i2,i3)*(1./4.*c020(i1+1,i2+1,i3)+1./8.*c010(i1+1,i2+1,i3))+1./8.*c010(i1,i2,i3)*c110(i1,i2+1,i3);
                                            cp( 2, 2)=1./16.*c110(i1,i2,i3)*c110(i1+1,i2+1,i3);
                                            ForStencil(m1,m2,m3)
                                            {
                                                int m  = M123(m1,m2,m3);
                                                coeffLocal(m,i1,i2,i3) += cLapSq*cp(m1,m2);
                                            }
                      // // save for testing : 
                      // ForStencil(m1,m2,m3)
                      // {
                      //   int m  = M123(m1,m2,m3);
                      //   lapSqCoeffLocal(m,i1,i2,i3) = cp(m1,m2);
                      // }          
                                        } // end if mask
                                    } // end for 3d 
                                } // end if ok
                            }
                            else
                            {
                // --------------- THREE DIMENSIONS Lap^2 ----------------
                                Index J1,J2,J3;
                                int extra=1;
                                getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
                                bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,J1,J2,J3);
                                if( ok )
                                {
                                    RealArray c200(J1,J2,J3), c110(J1,J2,J3), c020(J1,J2,J3), c100(J1,J2,J3), c010(J1,J2,J3);
                                    RealArray c002(J1,J2,J3), c101(J1,J2,J3), c011(J1,J2,J3), c001(J1,J2,J3);
                  // -- evaluate the SCALED coefficients of the Laplacian:
                  //    Delta = c200*Delta+r Delta-r + c110*Delta0r Delta0s + ...
                  // SEE ALSO
                  // research/compatibility/
                  //     lcbcDiff.mpl       : defines Qh 
                  //     writeLcbcFiles.mpl : 
                                    FOR_3D(i1,i2,i3,J1,J2,J3)
                                    {
                                        Real rx = DD(i1,i2,i3,0,0), ry = DD(i1,i2,i3,0,1), rz = DD(i1,i2,i3,0,2); 
                                        Real sx = DD(i1,i2,i3,1,0), sy = DD(i1,i2,i3,1,1), sz = DD(i1,i2,i3,1,2);         
                                        Real tx = DD(i1,i2,i3,2,0), ty = DD(i1,i2,i3,2,1), tz = DD(i1,i2,i3,2,2);         
                                        c200(i1,i2,i3) = ( SQR(rx) + SQR(ry) + SQR(rz) )/( SQR(dr[0]) );
                                        c020(i1,i2,i3) = ( SQR(sx) + SQR(sy) + SQR(sz) )/( SQR(dr[1]) );
                                        c002(i1,i2,i3) = ( SQR(tx) + SQR(ty) + SQR(tz) )/( SQR(dr[2]) );
                                        c110(i1,i2,i3) = 2.*( rx*sx + ry*sy +rz*sz )/( dr[0]*dr[1] ); 
                                        c101(i1,i2,i3) = 2.*( rx*tx + ry*ty +rz*tz )/( dr[0]*dr[2] ); 
                                        c011(i1,i2,i3) = 2.*( sx*tx + sy*ty +sz*tz )/( dr[1]*dr[2] ); 
                                        Real rxr = (DD(i1+1,i2,i3,0,0)-DD(i1-1,i2,i3,0,0))/(2.*dr[0]), ryr = (DD(i1+1,i2,i3,0,1)-DD(i1-1,i2,i3,0,1))/(2.*dr[0]), rzr = (DD(i1+1,i2,i3,0,2)-DD(i1-1,i2,i3,0,2))/(2.*dr[0]);
                                        Real sxr = (DD(i1+1,i2,i3,1,0)-DD(i1-1,i2,i3,1,0))/(2.*dr[0]), syr = (DD(i1+1,i2,i3,1,1)-DD(i1-1,i2,i3,1,1))/(2.*dr[0]), szr = (DD(i1+1,i2,i3,1,2)-DD(i1-1,i2,i3,1,2))/(2.*dr[0]);
                                        Real txr = (DD(i1+1,i2,i3,2,0)-DD(i1-1,i2,i3,2,0))/(2.*dr[0]), tyr = (DD(i1+1,i2,i3,2,1)-DD(i1-1,i2,i3,2,1))/(2.*dr[0]), tzr = (DD(i1+1,i2,i3,2,2)-DD(i1-1,i2,i3,2,2))/(2.*dr[0]);
                                        Real rxs = (DD(i1,i2+1,i3,0,0)-DD(i1,i2-1,i3,0,0))/(2.*dr[1]), rys = (DD(i1,i2+1,i3,0,1)-DD(i1,i2-1,i3,0,1))/(2.*dr[1]), rzs = (DD(i1,i2+1,i3,0,2)-DD(i1,i2-1,i3,0,2))/(2.*dr[1]);
                                        Real sxs = (DD(i1,i2+1,i3,1,0)-DD(i1,i2-1,i3,1,0))/(2.*dr[1]), sys = (DD(i1,i2+1,i3,1,1)-DD(i1,i2-1,i3,1,1))/(2.*dr[1]), szs = (DD(i1,i2+1,i3,1,2)-DD(i1,i2-1,i3,1,2))/(2.*dr[1]);
                                        Real txs = (DD(i1,i2+1,i3,2,0)-DD(i1,i2-1,i3,2,0))/(2.*dr[1]), tys = (DD(i1,i2+1,i3,2,1)-DD(i1,i2-1,i3,2,1))/(2.*dr[1]), tzs = (DD(i1,i2+1,i3,2,2)-DD(i1,i2-1,i3,2,2))/(2.*dr[1]);
                                        Real rxt = (DD(i1,i2,i3+1,0,0)-DD(i1,i2,i3-1,0,0))/(2.*dr[2]), ryt = (DD(i1,i2,i3+1,0,1)-DD(i1,i2,i3-1,0,1))/(2.*dr[2]), rzt = (DD(i1,i2,i3+1,0,2)-DD(i1,i2,i3-1,0,2))/(2.*dr[2]);
                                        Real sxt = (DD(i1,i2,i3+1,1,0)-DD(i1,i2,i3-1,1,0))/(2.*dr[2]), syt = (DD(i1,i2,i3+1,1,1)-DD(i1,i2,i3-1,1,1))/(2.*dr[2]), szt = (DD(i1,i2,i3+1,1,2)-DD(i1,i2,i3-1,1,2))/(2.*dr[2]);
                                        Real txt = (DD(i1,i2,i3+1,2,0)-DD(i1,i2,i3-1,2,0))/(2.*dr[2]), tyt = (DD(i1,i2,i3+1,2,1)-DD(i1,i2,i3-1,2,1))/(2.*dr[2]), tzt = (DD(i1,i2,i3+1,2,2)-DD(i1,i2,i3-1,2,2))/(2.*dr[2]);                
                                        Real rxx= rx*rxr + sx*rxs + tx*rxt;
                                        Real sxx= rx*sxr + sx*sxs + tx*sxt;
                                        Real txx= rx*txr + sx*txs + tx*txt;
                                        Real ryy= ry*ryr + sy*rys + ty*ryt;
                                        Real syy= ry*syr + sy*sys + ty*syt;
                                        Real tyy= ry*tyr + sy*tys + ty*tyt;
                                        Real rzz= rz*rzr + sz*rzs + tz*rzt;
                                        Real szz= rz*szr + sz*szs + tz*szt;
                                        Real tzz= rz*tzr + sz*tzs + tz*tzt;
                    // Real rxx= rx*rxr + sx*rxs + tx*rxt, ryx= rx*ryr + sx*rys + tx*ryt, rzx= rx*rzr + sx*rzs + tx*rzt;
                    // Real sxx= rx*sxr + sx*sxs + tx*sxt, syx= rx*syr + sx*sys + tx*syt, szx= rx*szr + sx*szs + tx*szt;
                    // Real txx= rx*txr + sx*txs + tx*txt, tyx= rx*tyr + sx*tys + tx*tyt, tzx= rx*tzr + sx*tzs + tx*tzt;                
                    // Real rxy= ry*rxr + sy*rxs + ty*rxt, ryy= ry*ryr + sy*rys + ty*ryt, rzy= ry*rzr + sy*rzs + ty*rzt;
                    // Real sxy= ry*sxr + sy*sxs + ty*sxt, syy= ry*syr + sy*sys + ty*syt, szy= ry*szr + sy*szs + ty*szt;
                    // Real txy= ry*txr + sy*txs + ty*txt, tyy= ry*tyr + sy*tys + ty*tyt, tzy= ry*tzr + sy*tzs + ty*tzt;  
                    // Real rxz= rz*rxr + sz*rxs + tz*rxt, ryz= rz*ryr + sz*rys + tz*ryt, rzz= rz*rzr + sz*rzs + tz*rzt;
                    // Real sxz= rz*sxr + sz*sxs + tz*sxt, syz= rz*syr + sz*sys + tz*syt, szz= rz*szr + sz*szs + tz*szt;
                    // Real txz= rz*txr + sz*txs + tz*txt, tyz= rz*tyr + sz*tys + tz*tyt, tzz= rz*tzr + sz*tzs + tz*tzt;  
                                        c100(i1,i2,i3) = (rxx + ryy + rzz)/(dr[0]); 
                                        c010(i1,i2,i3) = (sxx + syy + szz)/(dr[1]); 
                                        c001(i1,i2,i3) = (txx + tyy + tzz)/(dr[2]); 
                                    }
                                    Range R5(-2,2);
                                    RealArray cp(R5,R5,R5);
                                    FOR_3(i1,i2,i3,I1,I2,I3)
                                    {
                    // coefficients of Lap^2 (from research/compatibility/writeLcbcFiles.mpl --> lcbcEquationsDirichlet3dOrder4.h)
                                        if( maskLocal(i1,i2,i3)>0 )
                                        {
                                            cp(-2,-2,-2)=0;
                                            cp(-1,-2,-2)=0;
                                            cp( 0,-2,-2)=1./16.*c011(i1,i2,i3)*c011(i1,i2-1,i3-1);
                                            cp( 1,-2,-2)=0;
                                            cp( 2,-2,-2)=0;
                                            cp(-2,-1,-2)=0;
                                            cp(-1,-1,-2)=1./16.*c101(i1,i2,i3)*c011(i1-1,i2,i3-1)+1./16.*c011(i1,i2,i3)*c101(i1,i2-1,i3-1);
                                            cp( 0,-1,-2)=1./4.*c002(i1,i2,i3)*c011(i1,i2,i3-1)+c011(i1,i2,i3)*(-1./8.*c001(i1,i2-1,i3-1)+1./4.*c002(i1,i2-1,i3-1))-1./8.*c001(i1,i2,i3)*c011(i1,i2,i3-1);
                                            cp( 1,-1,-2)=-1./16.*c101(i1,i2,i3)*c011(i1+1,i2,i3-1)-1./16.*c011(i1,i2,i3)*c101(i1,i2-1,i3-1);
                                            cp( 2,-1,-2)=0;
                                            cp(-2, 0,-2)=1./16.*c101(i1,i2,i3)*c101(i1-1,i2,i3-1);
                                            cp(-1, 0,-2)=1./4.*c002(i1,i2,i3)*c101(i1,i2,i3-1)+c101(i1,i2,i3)*(-1./8.*c001(i1-1,i2,i3-1)+1./4.*c002(i1-1,i2,i3-1))-1./8.*c001(i1,i2,i3)*c101(i1,i2,i3-1);
                                            cp( 0, 0,-2)=c002(i1,i2,i3)*(c002(i1,i2,i3-1)-1./2.*c001(i1,i2,i3-1))+c101(i1,i2,i3)*(-1./16.*c101(i1+1,i2,i3-1)-1./16.*c101(i1-1,i2,i3-1))+c011(i1,i2,i3)*(-1./16.*c011(i1,i2+1,i3-1)-1./16.*c011(i1,i2-1,i3-1))+c001(i1,i2,i3)*(1./4.*c001(i1,i2,i3-1)-1./2.*c002(i1,i2,i3-1));
                                            cp( 1, 0,-2)=-1./4.*c002(i1,i2,i3)*c101(i1,i2,i3-1)+c101(i1,i2,i3)*(1./8.*c001(i1+1,i2,i3-1)-1./4.*c002(i1+1,i2,i3-1))+1./8.*c001(i1,i2,i3)*c101(i1,i2,i3-1);
                                            cp( 2, 0,-2)=1./16.*c101(i1,i2,i3)*c101(i1+1,i2,i3-1);
                                            cp(-2, 1,-2)=0;
                                            cp(-1, 1,-2)=-1./16.*c101(i1,i2,i3)*c011(i1-1,i2,i3-1)-1./16.*c011(i1,i2,i3)*c101(i1,i2+1,i3-1);
                                            cp( 0, 1,-2)=-1./4.*c002(i1,i2,i3)*c011(i1,i2,i3-1)+c011(i1,i2,i3)*(1./8.*c001(i1,i2+1,i3-1)-1./4.*c002(i1,i2+1,i3-1))+1./8.*c001(i1,i2,i3)*c011(i1,i2,i3-1);
                                            cp( 1, 1,-2)=1./16.*c101(i1,i2,i3)*c011(i1+1,i2,i3-1)+1./16.*c011(i1,i2,i3)*c101(i1,i2+1,i3-1);
                                            cp( 2, 1,-2)=0;
                                            cp(-2, 2,-2)=0;
                                            cp(-1, 2,-2)=0;
                                            cp( 0, 2,-2)=1./16.*c011(i1,i2,i3)*c011(i1,i2+1,i3-1);
                                            cp( 1, 2,-2)=0;
                                            cp( 2, 2,-2)=0;
                                            cp(-2,-2,-1)=0;
                                            cp(-1,-2,-1)=1./16.*c110(i1,i2,i3)*c011(i1-1,i2-1,i3)+1./16.*c011(i1,i2,i3)*c110(i1,i2-1,i3-1);
                                            cp( 0,-2,-1)=1./4.*c020(i1,i2,i3)*c011(i1,i2-1,i3)+c011(i1,i2,i3)*(-1./8.*c010(i1,i2-1,i3-1)+1./4.*c020(i1,i2-1,i3-1))-1./8.*c010(i1,i2,i3)*c011(i1,i2-1,i3);
                                            cp( 1,-2,-1)=-1./16.*c110(i1,i2,i3)*c011(i1+1,i2-1,i3)-1./16.*c011(i1,i2,i3)*c110(i1,i2-1,i3-1);
                                            cp( 2,-2,-1)=0;
                                            cp(-2,-1,-1)=1./16.*c110(i1,i2,i3)*c101(i1-1,i2-1,i3)+1./16.*c101(i1,i2,i3)*c110(i1-1,i2,i3-1);
                                            cp(-1,-1,-1)=1./4.*c200(i1,i2,i3)*c011(i1-1,i2,i3)+1./4.*c020(i1,i2,i3)*c101(i1,i2-1,i3)+1./4.*c002(i1,i2,i3)*c110(i1,i2,i3-1)+c110(i1,i2,i3)*(1./4.*c002(i1-1,i2-1,i3)-1./8.*c001(i1-1,i2-1,i3))+c101(i1,i2,i3)*(-1./8.*c010(i1-1,i2,i3-1)+1./4.*c020(i1-1,i2,i3-1))+c011(i1,i2,i3)*(1./4.*c200(i1,i2-1,i3-1)-1./8.*c100(i1,i2-1,i3-1))-1./8.*c100(i1,i2,i3)*c011(i1-1,i2,i3)-1./8.*c010(i1,i2,i3)*c101(i1,i2-1,i3)-1./8.*c001(i1,i2,i3)*c110(i1,i2,i3-1);
                                            cp( 0,-1,-1)=-1./2.*c200(i1,i2,i3)*c011(i1,i2,i3)+c020(i1,i2,i3)*(c002(i1,i2-1,i3)-1./2.*c001(i1,i2-1,i3)-1./2.*c011(i1,i2,i3))+c002(i1,i2,i3)*(-1./2.*c010(i1,i2,i3-1)+c020(i1,i2,i3-1)-1./2.*c011(i1,i2,i3))+c110(i1,i2,i3)*(-1./16.*c101(i1+1,i2-1,i3)-1./16.*c101(i1-1,i2-1,i3))+c101(i1,i2,i3)*(-1./16.*c110(i1+1,i2,i3-1)-1./16.*c110(i1-1,i2,i3-1))+c011(i1,i2,i3)*(-1./2.*c002(i1,i2-1,i3-1)-1./2.*c020(i1,i2-1,i3-1)-1./2.*c200(i1,i2-1,i3-1))+c010(i1,i2,i3)*(1./4.*c001(i1,i2-1,i3)-1./2.*c002(i1,i2-1,i3))+c001(i1,i2,i3)*(1./4.*c010(i1,i2,i3-1)-1./2.*c020(i1,i2,i3-1));
                                            cp( 1,-1,-1)=1./4.*c200(i1,i2,i3)*c011(i1+1,i2,i3)-1./4.*c020(i1,i2,i3)*c101(i1,i2-1,i3)-1./4.*c002(i1,i2,i3)*c110(i1,i2,i3-1)+c110(i1,i2,i3)*(-1./4.*c002(i1+1,i2-1,i3)+1./8.*c001(i1+1,i2-1,i3))+c101(i1,i2,i3)*(1./8.*c010(i1+1,i2,i3-1)-1./4.*c020(i1+1,i2,i3-1))+c011(i1,i2,i3)*(1./8.*c100(i1,i2-1,i3-1)+1./4.*c200(i1,i2-1,i3-1))+1./8.*c100(i1,i2,i3)*c011(i1+1,i2,i3)+1./8.*c010(i1,i2,i3)*c101(i1,i2-1,i3)+1./8.*c001(i1,i2,i3)*c110(i1,i2,i3-1);
                                            cp( 2,-1,-1)=1./16.*c110(i1,i2,i3)*c101(i1+1,i2-1,i3)+1./16.*c101(i1,i2,i3)*c110(i1+1,i2,i3-1);
                                            cp(-2, 0,-1)=1./4.*c200(i1,i2,i3)*c101(i1-1,i2,i3)+c101(i1,i2,i3)*(1./4.*c200(i1-1,i2,i3-1)-1./8.*c100(i1-1,i2,i3-1))-1./8.*c100(i1,i2,i3)*c101(i1-1,i2,i3);
                                            cp(-1, 0,-1)=c200(i1,i2,i3)*(-1./2.*c001(i1-1,i2,i3)-1./2.*c101(i1,i2,i3)+c002(i1-1,i2,i3))-1./2.*c020(i1,i2,i3)*c101(i1,i2,i3)+c002(i1,i2,i3)*(-1./2.*c100(i1,i2,i3-1)+c200(i1,i2,i3-1)-1./2.*c101(i1,i2,i3))+c110(i1,i2,i3)*(-1./16.*c011(i1-1,i2+1,i3)-1./16.*c011(i1-1,i2-1,i3))+c101(i1,i2,i3)*(-1./2.*c200(i1-1,i2,i3-1)-1./2.*c002(i1-1,i2,i3-1)-1./2.*c020(i1-1,i2,i3-1))+c011(i1,i2,i3)*(-1./16.*c110(i1,i2+1,i3-1)-1./16.*c110(i1,i2-1,i3-1))+c100(i1,i2,i3)*(1./4.*c001(i1-1,i2,i3)-1./2.*c002(i1-1,i2,i3))+c001(i1,i2,i3)*(-1./2.*c200(i1,i2,i3-1)+1./4.*c100(i1,i2,i3-1));
                                            cp( 0, 0,-1)=c200(i1,i2,i3)*(c001(i1,i2,i3)-2*c002(i1,i2,i3)-1./4.*c101(i1-1,i2,i3)+1./4.*c101(i1+1,i2,i3))+c020(i1,i2,i3)*(c001(i1,i2,i3)-2*c002(i1,i2,i3)-1./4.*c011(i1,i2-1,i3)+1./4.*c011(i1,i2+1,i3))+c002(i1,i2,i3)*(-2*c002(i1,i2,i3-1)-2*c020(i1,i2,i3-1)-2*c200(i1,i2,i3-1)+c001(i1,i2,i3)-2*c002(i1,i2,i3))+c101(i1,i2,i3)*(-1./4.*c200(i1+1,i2,i3-1)+1./4.*c200(i1-1,i2,i3-1)+1./8.*c100(i1+1,i2,i3-1)+1./8.*c100(i1-1,i2,i3-1))+c011(i1,i2,i3)*(1./8.*c010(i1,i2-1,i3-1)+1./8.*c010(i1,i2+1,i3-1)-1./4.*c020(i1,i2+1,i3-1)+1./4.*c020(i1,i2-1,i3-1))+c100(i1,i2,i3)*(1./8.*c101(i1+1,i2,i3)+1./8.*c101(i1-1,i2,i3))+c010(i1,i2,i3)*(1./8.*c011(i1,i2+1,i3)+1./8.*c011(i1,i2-1,i3))+c001(i1,i2,i3)*(c200(i1,i2,i3-1)+c020(i1,i2,i3-1)+c002(i1,i2,i3-1));
                                            cp( 1, 0,-1)=c200(i1,i2,i3)*(1./2.*c101(i1,i2,i3)-1./2.*c001(i1+1,i2,i3)+c002(i1+1,i2,i3))+1./2.*c020(i1,i2,i3)*c101(i1,i2,i3)+c002(i1,i2,i3)*(1./2.*c100(i1,i2,i3-1)+c200(i1,i2,i3-1)+1./2.*c101(i1,i2,i3))+c110(i1,i2,i3)*(1./16.*c011(i1+1,i2-1,i3)+1./16.*c011(i1+1,i2+1,i3))+c101(i1,i2,i3)*(1./2.*c200(i1+1,i2,i3-1)+1./2.*c020(i1+1,i2,i3-1)+1./2.*c002(i1+1,i2,i3-1))+c011(i1,i2,i3)*(1./16.*c110(i1,i2+1,i3-1)+1./16.*c110(i1,i2-1,i3-1))+c100(i1,i2,i3)*(-1./4.*c001(i1+1,i2,i3)+1./2.*c002(i1+1,i2,i3))+c001(i1,i2,i3)*(-1./4.*c100(i1,i2,i3-1)-1./2.*c200(i1,i2,i3-1));
                                            cp( 2, 0,-1)=-1./4.*c200(i1,i2,i3)*c101(i1+1,i2,i3)+c101(i1,i2,i3)*(-1./4.*c200(i1+1,i2,i3-1)-1./8.*c100(i1+1,i2,i3-1))-1./8.*c100(i1,i2,i3)*c101(i1+1,i2,i3);
                                            cp(-2, 1,-1)=-1./16.*c110(i1,i2,i3)*c101(i1-1,i2+1,i3)-1./16.*c101(i1,i2,i3)*c110(i1-1,i2,i3-1);
                                            cp(-1, 1,-1)=-1./4.*c200(i1,i2,i3)*c011(i1-1,i2,i3)+1./4.*c020(i1,i2,i3)*c101(i1,i2+1,i3)-1./4.*c002(i1,i2,i3)*c110(i1,i2,i3-1)+c110(i1,i2,i3)*(-1./4.*c002(i1-1,i2+1,i3)+1./8.*c001(i1-1,i2+1,i3))+c101(i1,i2,i3)*(1./4.*c020(i1-1,i2,i3-1)+1./8.*c010(i1-1,i2,i3-1))+c011(i1,i2,i3)*(-1./4.*c200(i1,i2+1,i3-1)+1./8.*c100(i1,i2+1,i3-1))+1./8.*c100(i1,i2,i3)*c011(i1-1,i2,i3)+1./8.*c010(i1,i2,i3)*c101(i1,i2+1,i3)+1./8.*c001(i1,i2,i3)*c110(i1,i2,i3-1);
                                            cp( 0, 1,-1)=1./2.*c200(i1,i2,i3)*c011(i1,i2,i3)+c020(i1,i2,i3)*(1./2.*c011(i1,i2,i3)-1./2.*c001(i1,i2+1,i3)+c002(i1,i2+1,i3))+c002(i1,i2,i3)*(1./2.*c010(i1,i2,i3-1)+c020(i1,i2,i3-1)+1./2.*c011(i1,i2,i3))+c110(i1,i2,i3)*(1./16.*c101(i1+1,i2+1,i3)+1./16.*c101(i1-1,i2+1,i3))+c101(i1,i2,i3)*(1./16.*c110(i1+1,i2,i3-1)+1./16.*c110(i1-1,i2,i3-1))+c011(i1,i2,i3)*(1./2.*c200(i1,i2+1,i3-1)+1./2.*c002(i1,i2+1,i3-1)+1./2.*c020(i1,i2+1,i3-1))+c010(i1,i2,i3)*(-1./4.*c001(i1,i2+1,i3)+1./2.*c002(i1,i2+1,i3))+c001(i1,i2,i3)*(-1./2.*c020(i1,i2,i3-1)-1./4.*c010(i1,i2,i3-1));
                                            cp( 1, 1,-1)=-1./4.*c200(i1,i2,i3)*c011(i1+1,i2,i3)-1./4.*c020(i1,i2,i3)*c101(i1,i2+1,i3)+1./4.*c002(i1,i2,i3)*c110(i1,i2,i3-1)+c110(i1,i2,i3)*(1./4.*c002(i1+1,i2+1,i3)-1./8.*c001(i1+1,i2+1,i3))+c101(i1,i2,i3)*(-1./4.*c020(i1+1,i2,i3-1)-1./8.*c010(i1+1,i2,i3-1))+c011(i1,i2,i3)*(-1./8.*c100(i1,i2+1,i3-1)-1./4.*c200(i1,i2+1,i3-1))-1./8.*c100(i1,i2,i3)*c011(i1+1,i2,i3)-1./8.*c010(i1,i2,i3)*c101(i1,i2+1,i3)-1./8.*c001(i1,i2,i3)*c110(i1,i2,i3-1);
                                            cp( 2, 1,-1)=-1./16.*c110(i1,i2,i3)*c101(i1+1,i2+1,i3)-1./16.*c101(i1,i2,i3)*c110(i1+1,i2,i3-1);
                                            cp(-2, 2,-1)=0;
                                            cp(-1, 2,-1)=1./16.*c110(i1,i2,i3)*c011(i1-1,i2+1,i3)+1./16.*c011(i1,i2,i3)*c110(i1,i2+1,i3-1);
                                            cp( 0, 2,-1)=-1./4.*c020(i1,i2,i3)*c011(i1,i2+1,i3)+c011(i1,i2,i3)*(-1./8.*c010(i1,i2+1,i3-1)-1./4.*c020(i1,i2+1,i3-1))-1./8.*c010(i1,i2,i3)*c011(i1,i2+1,i3);
                                            cp( 1, 2,-1)=-1./16.*c110(i1,i2,i3)*c011(i1+1,i2+1,i3)-1./16.*c011(i1,i2,i3)*c110(i1,i2+1,i3-1);
                                            cp( 2, 2,-1)=0;
                                            cp(-2,-2, 0)=1./16.*c110(i1,i2,i3)*c110(i1-1,i2-1,i3);
                                            cp(-1,-2, 0)=1./4.*c020(i1,i2,i3)*c110(i1,i2-1,i3)+c110(i1,i2,i3)*(1./4.*c020(i1-1,i2-1,i3)-1./8.*c010(i1-1,i2-1,i3))-1./8.*c010(i1,i2,i3)*c110(i1,i2-1,i3);
                                            cp( 0,-2, 0)=c020(i1,i2,i3)*(-1./2.*c010(i1,i2-1,i3)+c020(i1,i2-1,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1+1,i2-1,i3)-1./16.*c110(i1-1,i2-1,i3))+c011(i1,i2,i3)*(-1./16.*c011(i1,i2-1,i3+1)-1./16.*c011(i1,i2-1,i3-1))+c010(i1,i2,i3)*(1./4.*c010(i1,i2-1,i3)-1./2.*c020(i1,i2-1,i3));
                                            cp( 1,-2, 0)=-1./4.*c020(i1,i2,i3)*c110(i1,i2-1,i3)+c110(i1,i2,i3)*(1./8.*c010(i1+1,i2-1,i3)-1./4.*c020(i1+1,i2-1,i3))+1./8.*c010(i1,i2,i3)*c110(i1,i2-1,i3);
                                            cp( 2,-2, 0)=1./16.*c110(i1,i2,i3)*c110(i1+1,i2-1,i3);
                                            cp(-2,-1, 0)=1./4.*c200(i1,i2,i3)*c110(i1-1,i2,i3)+c110(i1,i2,i3)*(-1./8.*c100(i1-1,i2-1,i3)+1./4.*c200(i1-1,i2-1,i3))-1./8.*c100(i1,i2,i3)*c110(i1-1,i2,i3);
                                            cp(-1,-1, 0)=c200(i1,i2,i3)*(c020(i1-1,i2,i3)-1./2.*c010(i1-1,i2,i3)-1./2.*c110(i1,i2,i3))+c020(i1,i2,i3)*(c200(i1,i2-1,i3)-1./2.*c100(i1,i2-1,i3)-1./2.*c110(i1,i2,i3))-1./2.*c002(i1,i2,i3)*c110(i1,i2,i3)+c110(i1,i2,i3)*(-1./2.*c200(i1-1,i2-1,i3)-1./2.*c020(i1-1,i2-1,i3)-1./2.*c002(i1-1,i2-1,i3))+c101(i1,i2,i3)*(-1./16.*c011(i1-1,i2,i3-1)-1./16.*c011(i1-1,i2,i3+1))+c011(i1,i2,i3)*(-1./16.*c101(i1,i2-1,i3+1)-1./16.*c101(i1,i2-1,i3-1))+c100(i1,i2,i3)*(-1./2.*c020(i1-1,i2,i3)+1./4.*c010(i1-1,i2,i3))+c010(i1,i2,i3)*(-1./2.*c200(i1,i2-1,i3)+1./4.*c100(i1,i2-1,i3));
                                            cp( 0,-1, 0)=c200(i1,i2,i3)*(c010(i1,i2,i3)-2*c020(i1,i2,i3)-1./4.*c110(i1-1,i2,i3)+1./4.*c110(i1+1,i2,i3))+c020(i1,i2,i3)*(-2*c002(i1,i2-1,i3)-2*c020(i1,i2-1,i3)-2*c200(i1,i2-1,i3)+c010(i1,i2,i3)-2*c020(i1,i2,i3))+c002(i1,i2,i3)*(c010(i1,i2,i3)-2*c020(i1,i2,i3)-1./4.*c011(i1,i2,i3-1)+1./4.*c011(i1,i2,i3+1))+c110(i1,i2,i3)*(1./8.*c100(i1-1,i2-1,i3)-1./4.*c200(i1+1,i2-1,i3)+1./4.*c200(i1-1,i2-1,i3)+1./8.*c100(i1+1,i2-1,i3))+c011(i1,i2,i3)*(1./8.*c001(i1,i2-1,i3-1)+1./8.*c001(i1,i2-1,i3+1)-1./4.*c002(i1,i2-1,i3+1)+1./4.*c002(i1,i2-1,i3-1))+c100(i1,i2,i3)*(1./8.*c110(i1+1,i2,i3)+1./8.*c110(i1-1,i2,i3))+c010(i1,i2,i3)*(c200(i1,i2-1,i3)+c002(i1,i2-1,i3)+c020(i1,i2-1,i3))+c001(i1,i2,i3)*(1./8.*c011(i1,i2,i3+1)+1./8.*c011(i1,i2,i3-1));
                                            cp( 1,-1, 0)=c200(i1,i2,i3)*(1./2.*c110(i1,i2,i3)-1./2.*c010(i1+1,i2,i3)+c020(i1+1,i2,i3))+c020(i1,i2,i3)*(1./2.*c110(i1,i2,i3)+c200(i1,i2-1,i3)+1./2.*c100(i1,i2-1,i3))+1./2.*c002(i1,i2,i3)*c110(i1,i2,i3)+c110(i1,i2,i3)*(1./2.*c200(i1+1,i2-1,i3)+1./2.*c020(i1+1,i2-1,i3)+1./2.*c002(i1+1,i2-1,i3))+c101(i1,i2,i3)*(1./16.*c011(i1+1,i2,i3-1)+1./16.*c011(i1+1,i2,i3+1))+c011(i1,i2,i3)*(1./16.*c101(i1,i2-1,i3+1)+1./16.*c101(i1,i2-1,i3-1))+c100(i1,i2,i3)*(-1./4.*c010(i1+1,i2,i3)+1./2.*c020(i1+1,i2,i3))+c010(i1,i2,i3)*(-1./4.*c100(i1,i2-1,i3)-1./2.*c200(i1,i2-1,i3));
                                            cp( 2,-1, 0)=-1./4.*c200(i1,i2,i3)*c110(i1+1,i2,i3)+c110(i1,i2,i3)*(-1./4.*c200(i1+1,i2-1,i3)-1./8.*c100(i1+1,i2-1,i3))-1./8.*c100(i1,i2,i3)*c110(i1+1,i2,i3);
                                            cp(-2, 0, 0)=c200(i1,i2,i3)*(c200(i1-1,i2,i3)-1./2.*c100(i1-1,i2,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1-1,i2+1,i3)-1./16.*c110(i1-1,i2-1,i3))+c101(i1,i2,i3)*(-1./16.*c101(i1-1,i2,i3-1)-1./16.*c101(i1-1,i2,i3+1))+c100(i1,i2,i3)*(-1./2.*c200(i1-1,i2,i3)+1./4.*c100(i1-1,i2,i3));
                                            cp(-1, 0, 0)=c200(i1,i2,i3)*(c100(i1,i2,i3)-2*c200(i1,i2,i3)-2*c002(i1-1,i2,i3)-2*c020(i1-1,i2,i3)-2*c200(i1-1,i2,i3))+c020(i1,i2,i3)*(c100(i1,i2,i3)-2*c200(i1,i2,i3)+1./4.*c110(i1,i2+1,i3)-1./4.*c110(i1,i2-1,i3))+c002(i1,i2,i3)*(c100(i1,i2,i3)-2*c200(i1,i2,i3)+1./4.*c101(i1,i2,i3+1)-1./4.*c101(i1,i2,i3-1))+c110(i1,i2,i3)*(1./8.*c010(i1-1,i2+1,i3)+1./4.*c020(i1-1,i2-1,i3)-1./4.*c020(i1-1,i2+1,i3)+1./8.*c010(i1-1,i2-1,i3))+c101(i1,i2,i3)*(1./8.*c001(i1-1,i2,i3-1)+1./8.*c001(i1-1,i2,i3+1)+1./4.*c002(i1-1,i2,i3-1)-1./4.*c002(i1-1,i2,i3+1))+c100(i1,i2,i3)*(c200(i1-1,i2,i3)+c020(i1-1,i2,i3)+c002(i1-1,i2,i3))+c010(i1,i2,i3)*(1./8.*c110(i1,i2+1,i3)+1./8.*c110(i1,i2-1,i3))+c001(i1,i2,i3)*(1./8.*c101(i1,i2,i3+1)+1./8.*c101(i1,i2,i3-1));
                                            cp( 0, 0, 0)=c200(i1,i2,i3)*(-1./2.*c100(i1+1,i2,i3)+c200(i1-1,i2,i3)+c200(i1+1,i2,i3)+1./2.*c100(i1-1,i2,i3)+4*c200(i1,i2,i3)+4*c002(i1,i2,i3)+4*c020(i1,i2,i3))+c020(i1,i2,i3)*(4*c200(i1,i2,i3)+4*c002(i1,i2,i3)+4*c020(i1,i2,i3)+c020(i1,i2+1,i3)+1./2.*c010(i1,i2-1,i3)-1./2.*c010(i1,i2+1,i3)+c020(i1,i2-1,i3))+c002(i1,i2,i3)*(4*c200(i1,i2,i3)+4*c002(i1,i2,i3)+4*c020(i1,i2,i3)-1./2.*c001(i1,i2,i3+1)+c002(i1,i2,i3-1)+c002(i1,i2,i3+1)+1./2.*c001(i1,i2,i3-1))+c110(i1,i2,i3)*(1./16.*c110(i1+1,i2+1,i3)+1./16.*c110(i1-1,i2+1,i3)+1./16.*c110(i1+1,i2-1,i3)+1./16.*c110(i1-1,i2-1,i3))+c101(i1,i2,i3)*(1./16.*c101(i1+1,i2,i3-1)+1./16.*c101(i1+1,i2,i3+1)+1./16.*c101(i1-1,i2,i3-1)+1./16.*c101(i1-1,i2,i3+1))+c011(i1,i2,i3)*(1./16.*c011(i1,i2-1,i3+1)+1./16.*c011(i1,i2+1,i3-1)+1./16.*c011(i1,i2+1,i3+1)+1./16.*c011(i1,i2-1,i3-1))+c100(i1,i2,i3)*(-1./2.*c200(i1-1,i2,i3)-1./4.*c100(i1-1,i2,i3)-1./4.*c100(i1+1,i2,i3)+1./2.*c200(i1+1,i2,i3))+c010(i1,i2,i3)*(-1./4.*c010(i1,i2+1,i3)-1./4.*c010(i1,i2-1,i3)-1./2.*c020(i1,i2-1,i3)+1./2.*c020(i1,i2+1,i3))+c001(i1,i2,i3)*(-1./4.*c001(i1,i2,i3-1)-1./4.*c001(i1,i2,i3+1)+1./2.*c002(i1,i2,i3+1)-1./2.*c002(i1,i2,i3-1));
                                            cp( 1, 0, 0)=c200(i1,i2,i3)*(-c100(i1,i2,i3)-2*c200(i1,i2,i3)-2*c002(i1+1,i2,i3)-2*c020(i1+1,i2,i3)-2*c200(i1+1,i2,i3))+c020(i1,i2,i3)*(-c100(i1,i2,i3)-2*c200(i1,i2,i3)-1./4.*c110(i1,i2+1,i3)+1./4.*c110(i1,i2-1,i3))+c002(i1,i2,i3)*(-c100(i1,i2,i3)-2*c200(i1,i2,i3)-1./4.*c101(i1,i2,i3+1)+1./4.*c101(i1,i2,i3-1))+c110(i1,i2,i3)*(-1./8.*c010(i1+1,i2+1,i3)-1./8.*c010(i1+1,i2-1,i3)+1./4.*c020(i1+1,i2+1,i3)-1./4.*c020(i1+1,i2-1,i3))+c101(i1,i2,i3)*(-1./8.*c001(i1+1,i2,i3+1)-1./8.*c001(i1+1,i2,i3-1)+1./4.*c002(i1+1,i2,i3+1)-1./4.*c002(i1+1,i2,i3-1))+c100(i1,i2,i3)*(-c020(i1+1,i2,i3)-c002(i1+1,i2,i3)-c200(i1+1,i2,i3))+c010(i1,i2,i3)*(-1./8.*c110(i1,i2+1,i3)-1./8.*c110(i1,i2-1,i3))+c001(i1,i2,i3)*(-1./8.*c101(i1,i2,i3+1)-1./8.*c101(i1,i2,i3-1));
                                            cp( 2, 0, 0)=c200(i1,i2,i3)*(1./2.*c100(i1+1,i2,i3)+c200(i1+1,i2,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1+1,i2+1,i3)-1./16.*c110(i1+1,i2-1,i3))+c101(i1,i2,i3)*(-1./16.*c101(i1+1,i2,i3-1)-1./16.*c101(i1+1,i2,i3+1))+c100(i1,i2,i3)*(1./4.*c100(i1+1,i2,i3)+1./2.*c200(i1+1,i2,i3));
                                            cp(-2, 1, 0)=-1./4.*c200(i1,i2,i3)*c110(i1-1,i2,i3)+c110(i1,i2,i3)*(-1./4.*c200(i1-1,i2+1,i3)+1./8.*c100(i1-1,i2+1,i3))+1./8.*c100(i1,i2,i3)*c110(i1-1,i2,i3);
                                            cp(-1, 1, 0)=c200(i1,i2,i3)*(c020(i1-1,i2,i3)+1./2.*c010(i1-1,i2,i3)+1./2.*c110(i1,i2,i3))+c020(i1,i2,i3)*(-1./2.*c100(i1,i2+1,i3)+c200(i1,i2+1,i3)+1./2.*c110(i1,i2,i3))+1./2.*c002(i1,i2,i3)*c110(i1,i2,i3)+c110(i1,i2,i3)*(1./2.*c200(i1-1,i2+1,i3)+1./2.*c020(i1-1,i2+1,i3)+1./2.*c002(i1-1,i2+1,i3))+c101(i1,i2,i3)*(1./16.*c011(i1-1,i2,i3+1)+1./16.*c011(i1-1,i2,i3-1))+c011(i1,i2,i3)*(1./16.*c101(i1,i2+1,i3-1)+1./16.*c101(i1,i2+1,i3+1))+c100(i1,i2,i3)*(-1./4.*c010(i1-1,i2,i3)-1./2.*c020(i1-1,i2,i3))+c010(i1,i2,i3)*(-1./4.*c100(i1,i2+1,i3)+1./2.*c200(i1,i2+1,i3));
                                            cp( 0, 1, 0)=c200(i1,i2,i3)*(-c010(i1,i2,i3)-2*c020(i1,i2,i3)+1./4.*c110(i1-1,i2,i3)-1./4.*c110(i1+1,i2,i3))+c020(i1,i2,i3)*(-2*c002(i1,i2+1,i3)-2*c020(i1,i2+1,i3)-2*c200(i1,i2+1,i3)-c010(i1,i2,i3)-2*c020(i1,i2,i3))+c002(i1,i2,i3)*(-c010(i1,i2,i3)-2*c020(i1,i2,i3)+1./4.*c011(i1,i2,i3-1)-1./4.*c011(i1,i2,i3+1))+c110(i1,i2,i3)*(-1./4.*c200(i1-1,i2+1,i3)-1./8.*c100(i1-1,i2+1,i3)-1./8.*c100(i1+1,i2+1,i3)+1./4.*c200(i1+1,i2+1,i3))+c011(i1,i2,i3)*(-1./8.*c001(i1,i2+1,i3-1)+1./4.*c002(i1,i2+1,i3+1)-1./4.*c002(i1,i2+1,i3-1)-1./8.*c001(i1,i2+1,i3+1))+c100(i1,i2,i3)*(-1./8.*c110(i1+1,i2,i3)-1./8.*c110(i1-1,i2,i3))+c010(i1,i2,i3)*(-c200(i1,i2+1,i3)-c002(i1,i2+1,i3)-c020(i1,i2+1,i3))+c001(i1,i2,i3)*(-1./8.*c011(i1,i2,i3+1)-1./8.*c011(i1,i2,i3-1));
                                            cp( 1, 1, 0)=c200(i1,i2,i3)*(-1./2.*c110(i1,i2,i3)+1./2.*c010(i1+1,i2,i3)+c020(i1+1,i2,i3))+c020(i1,i2,i3)*(1./2.*c100(i1,i2+1,i3)+c200(i1,i2+1,i3)-1./2.*c110(i1,i2,i3))-1./2.*c002(i1,i2,i3)*c110(i1,i2,i3)+c110(i1,i2,i3)*(-1./2.*c020(i1+1,i2+1,i3)-1./2.*c200(i1+1,i2+1,i3)-1./2.*c002(i1+1,i2+1,i3))+c101(i1,i2,i3)*(-1./16.*c011(i1+1,i2,i3-1)-1./16.*c011(i1+1,i2,i3+1))+c011(i1,i2,i3)*(-1./16.*c101(i1,i2+1,i3-1)-1./16.*c101(i1,i2+1,i3+1))+c100(i1,i2,i3)*(1./4.*c010(i1+1,i2,i3)+1./2.*c020(i1+1,i2,i3))+c010(i1,i2,i3)*(1./4.*c100(i1,i2+1,i3)+1./2.*c200(i1,i2+1,i3));
                                            cp( 2, 1, 0)=1./4.*c200(i1,i2,i3)*c110(i1+1,i2,i3)+c110(i1,i2,i3)*(1./8.*c100(i1+1,i2+1,i3)+1./4.*c200(i1+1,i2+1,i3))+1./8.*c100(i1,i2,i3)*c110(i1+1,i2,i3);
                                            cp(-2, 2, 0)=1./16.*c110(i1,i2,i3)*c110(i1-1,i2+1,i3);
                                            cp(-1, 2, 0)=-1./4.*c020(i1,i2,i3)*c110(i1,i2+1,i3)+c110(i1,i2,i3)*(-1./8.*c010(i1-1,i2+1,i3)-1./4.*c020(i1-1,i2+1,i3))-1./8.*c010(i1,i2,i3)*c110(i1,i2+1,i3);
                                            cp( 0, 2, 0)=c020(i1,i2,i3)*(c020(i1,i2+1,i3)+1./2.*c010(i1,i2+1,i3))+c110(i1,i2,i3)*(-1./16.*c110(i1+1,i2+1,i3)-1./16.*c110(i1-1,i2+1,i3))+c011(i1,i2,i3)*(-1./16.*c011(i1,i2+1,i3-1)-1./16.*c011(i1,i2+1,i3+1))+c010(i1,i2,i3)*(1./4.*c010(i1,i2+1,i3)+1./2.*c020(i1,i2+1,i3));
                                            cp( 1, 2, 0)=1./4.*c020(i1,i2,i3)*c110(i1,i2+1,i3)+c110(i1,i2,i3)*(1./8.*c010(i1+1,i2+1,i3)+1./4.*c020(i1+1,i2+1,i3))+1./8.*c010(i1,i2,i3)*c110(i1,i2+1,i3);
                                            cp( 2, 2, 0)=1./16.*c110(i1,i2,i3)*c110(i1+1,i2+1,i3);
                                            cp(-2,-2, 1)=0;
                                            cp(-1,-2, 1)=-1./16.*c110(i1,i2,i3)*c011(i1-1,i2-1,i3)-1./16.*c011(i1,i2,i3)*c110(i1,i2-1,i3+1);
                                            cp( 0,-2, 1)=-1./4.*c020(i1,i2,i3)*c011(i1,i2-1,i3)+c011(i1,i2,i3)*(1./8.*c010(i1,i2-1,i3+1)-1./4.*c020(i1,i2-1,i3+1))+1./8.*c010(i1,i2,i3)*c011(i1,i2-1,i3);
                                            cp( 1,-2, 1)=1./16.*c110(i1,i2,i3)*c011(i1+1,i2-1,i3)+1./16.*c011(i1,i2,i3)*c110(i1,i2-1,i3+1);
                                            cp( 2,-2, 1)=0;
                                            cp(-2,-1, 1)=-1./16.*c110(i1,i2,i3)*c101(i1-1,i2-1,i3)-1./16.*c101(i1,i2,i3)*c110(i1-1,i2,i3+1);
                                            cp(-1,-1, 1)=-1./4.*c200(i1,i2,i3)*c011(i1-1,i2,i3)-1./4.*c020(i1,i2,i3)*c101(i1,i2-1,i3)+1./4.*c002(i1,i2,i3)*c110(i1,i2,i3+1)+c110(i1,i2,i3)*(1./8.*c001(i1-1,i2-1,i3)+1./4.*c002(i1-1,i2-1,i3))+c101(i1,i2,i3)*(1./8.*c010(i1-1,i2,i3+1)-1./4.*c020(i1-1,i2,i3+1))+c011(i1,i2,i3)*(1./8.*c100(i1,i2-1,i3+1)-1./4.*c200(i1,i2-1,i3+1))+1./8.*c100(i1,i2,i3)*c011(i1-1,i2,i3)+1./8.*c010(i1,i2,i3)*c101(i1,i2-1,i3)+1./8.*c001(i1,i2,i3)*c110(i1,i2,i3+1);
                                            cp( 0,-1, 1)=1./2.*c200(i1,i2,i3)*c011(i1,i2,i3)+c020(i1,i2,i3)*(c002(i1,i2-1,i3)+1./2.*c001(i1,i2-1,i3)+1./2.*c011(i1,i2,i3))+c002(i1,i2,i3)*(-1./2.*c010(i1,i2,i3+1)+c020(i1,i2,i3+1)+1./2.*c011(i1,i2,i3))+c110(i1,i2,i3)*(1./16.*c101(i1+1,i2-1,i3)+1./16.*c101(i1-1,i2-1,i3))+c101(i1,i2,i3)*(1./16.*c110(i1+1,i2,i3+1)+1./16.*c110(i1-1,i2,i3+1))+c011(i1,i2,i3)*(1./2.*c200(i1,i2-1,i3+1)+1./2.*c002(i1,i2-1,i3+1)+1./2.*c020(i1,i2-1,i3+1))+c010(i1,i2,i3)*(-1./2.*c002(i1,i2-1,i3)-1./4.*c001(i1,i2-1,i3))+c001(i1,i2,i3)*(-1./4.*c010(i1,i2,i3+1)+1./2.*c020(i1,i2,i3+1));
                                            cp( 1,-1, 1)=-1./4.*c200(i1,i2,i3)*c011(i1+1,i2,i3)+1./4.*c020(i1,i2,i3)*c101(i1,i2-1,i3)-1./4.*c002(i1,i2,i3)*c110(i1,i2,i3+1)+c110(i1,i2,i3)*(-1./8.*c001(i1+1,i2-1,i3)-1./4.*c002(i1+1,i2-1,i3))+c101(i1,i2,i3)*(-1./8.*c010(i1+1,i2,i3+1)+1./4.*c020(i1+1,i2,i3+1))+c011(i1,i2,i3)*(-1./8.*c100(i1,i2-1,i3+1)-1./4.*c200(i1,i2-1,i3+1))-1./8.*c100(i1,i2,i3)*c011(i1+1,i2,i3)-1./8.*c010(i1,i2,i3)*c101(i1,i2-1,i3)-1./8.*c001(i1,i2,i3)*c110(i1,i2,i3+1);
                                            cp( 2,-1, 1)=-1./16.*c110(i1,i2,i3)*c101(i1+1,i2-1,i3)-1./16.*c101(i1,i2,i3)*c110(i1+1,i2,i3+1);
                                            cp(-2, 0, 1)=-1./4.*c200(i1,i2,i3)*c101(i1-1,i2,i3)+c101(i1,i2,i3)*(-1./4.*c200(i1-1,i2,i3+1)+1./8.*c100(i1-1,i2,i3+1))+1./8.*c100(i1,i2,i3)*c101(i1-1,i2,i3);
                                            cp(-1, 0, 1)=c200(i1,i2,i3)*(1./2.*c001(i1-1,i2,i3)+1./2.*c101(i1,i2,i3)+c002(i1-1,i2,i3))+1./2.*c020(i1,i2,i3)*c101(i1,i2,i3)+c002(i1,i2,i3)*(-1./2.*c100(i1,i2,i3+1)+1./2.*c101(i1,i2,i3)+c200(i1,i2,i3+1))+c110(i1,i2,i3)*(1./16.*c011(i1-1,i2+1,i3)+1./16.*c011(i1-1,i2-1,i3))+c101(i1,i2,i3)*(1./2.*c200(i1-1,i2,i3+1)+1./2.*c002(i1-1,i2,i3+1)+1./2.*c020(i1-1,i2,i3+1))+c011(i1,i2,i3)*(1./16.*c110(i1,i2+1,i3+1)+1./16.*c110(i1,i2-1,i3+1))+c100(i1,i2,i3)*(-1./2.*c002(i1-1,i2,i3)-1./4.*c001(i1-1,i2,i3))+c001(i1,i2,i3)*(-1./4.*c100(i1,i2,i3+1)+1./2.*c200(i1,i2,i3+1));
                                            cp( 0, 0, 1)=c200(i1,i2,i3)*(-c001(i1,i2,i3)-2*c002(i1,i2,i3)+1./4.*c101(i1-1,i2,i3)-1./4.*c101(i1+1,i2,i3))+c020(i1,i2,i3)*(-c001(i1,i2,i3)-2*c002(i1,i2,i3)+1./4.*c011(i1,i2-1,i3)-1./4.*c011(i1,i2+1,i3))+c002(i1,i2,i3)*(-2*c002(i1,i2,i3+1)-2*c020(i1,i2,i3+1)-2*c200(i1,i2,i3+1)-c001(i1,i2,i3)-2*c002(i1,i2,i3))+c101(i1,i2,i3)*(-1./4.*c200(i1-1,i2,i3+1)+1./4.*c200(i1+1,i2,i3+1)-1./8.*c100(i1-1,i2,i3+1)-1./8.*c100(i1+1,i2,i3+1))+c011(i1,i2,i3)*(1./4.*c020(i1,i2+1,i3+1)-1./8.*c010(i1,i2-1,i3+1)-1./8.*c010(i1,i2+1,i3+1)-1./4.*c020(i1,i2-1,i3+1))+c100(i1,i2,i3)*(-1./8.*c101(i1+1,i2,i3)-1./8.*c101(i1-1,i2,i3))+c010(i1,i2,i3)*(-1./8.*c011(i1,i2+1,i3)-1./8.*c011(i1,i2-1,i3))+c001(i1,i2,i3)*(-c200(i1,i2,i3+1)-c020(i1,i2,i3+1)-c002(i1,i2,i3+1));
                                            cp( 1, 0, 1)=c200(i1,i2,i3)*(-1./2.*c101(i1,i2,i3)+1./2.*c001(i1+1,i2,i3)+c002(i1+1,i2,i3))-1./2.*c020(i1,i2,i3)*c101(i1,i2,i3)+c002(i1,i2,i3)*(-1./2.*c101(i1,i2,i3)+1./2.*c100(i1,i2,i3+1)+c200(i1,i2,i3+1))+c110(i1,i2,i3)*(-1./16.*c011(i1+1,i2-1,i3)-1./16.*c011(i1+1,i2+1,i3))+c101(i1,i2,i3)*(-1./2.*c200(i1+1,i2,i3+1)-1./2.*c020(i1+1,i2,i3+1)-1./2.*c002(i1+1,i2,i3+1))+c011(i1,i2,i3)*(-1./16.*c110(i1,i2+1,i3+1)-1./16.*c110(i1,i2-1,i3+1))+c100(i1,i2,i3)*(1./2.*c002(i1+1,i2,i3)+1./4.*c001(i1+1,i2,i3))+c001(i1,i2,i3)*(1./4.*c100(i1,i2,i3+1)+1./2.*c200(i1,i2,i3+1));
                                            cp( 2, 0, 1)=1./4.*c200(i1,i2,i3)*c101(i1+1,i2,i3)+c101(i1,i2,i3)*(1./4.*c200(i1+1,i2,i3+1)+1./8.*c100(i1+1,i2,i3+1))+1./8.*c100(i1,i2,i3)*c101(i1+1,i2,i3);
                                            cp(-2, 1, 1)=1./16.*c110(i1,i2,i3)*c101(i1-1,i2+1,i3)+1./16.*c101(i1,i2,i3)*c110(i1-1,i2,i3+1);
                                            cp(-1, 1, 1)=1./4.*c200(i1,i2,i3)*c011(i1-1,i2,i3)-1./4.*c020(i1,i2,i3)*c101(i1,i2+1,i3)-1./4.*c002(i1,i2,i3)*c110(i1,i2,i3+1)+c110(i1,i2,i3)*(-1./8.*c001(i1-1,i2+1,i3)-1./4.*c002(i1-1,i2+1,i3))+c101(i1,i2,i3)*(-1./8.*c010(i1-1,i2,i3+1)-1./4.*c020(i1-1,i2,i3+1))+c011(i1,i2,i3)*(-1./8.*c100(i1,i2+1,i3+1)+1./4.*c200(i1,i2+1,i3+1))-1./8.*c100(i1,i2,i3)*c011(i1-1,i2,i3)-1./8.*c010(i1,i2,i3)*c101(i1,i2+1,i3)-1./8.*c001(i1,i2,i3)*c110(i1,i2,i3+1);
                                            cp( 0, 1, 1)=-1./2.*c200(i1,i2,i3)*c011(i1,i2,i3)+c020(i1,i2,i3)*(-1./2.*c011(i1,i2,i3)+1./2.*c001(i1,i2+1,i3)+c002(i1,i2+1,i3))+c002(i1,i2,i3)*(-1./2.*c011(i1,i2,i3)+c020(i1,i2,i3+1)+1./2.*c010(i1,i2,i3+1))+c110(i1,i2,i3)*(-1./16.*c101(i1+1,i2+1,i3)-1./16.*c101(i1-1,i2+1,i3))+c101(i1,i2,i3)*(-1./16.*c110(i1+1,i2,i3+1)-1./16.*c110(i1-1,i2,i3+1))+c011(i1,i2,i3)*(-1./2.*c020(i1,i2+1,i3+1)-1./2.*c200(i1,i2+1,i3+1)-1./2.*c002(i1,i2+1,i3+1))+c010(i1,i2,i3)*(1./2.*c002(i1,i2+1,i3)+1./4.*c001(i1,i2+1,i3))+c001(i1,i2,i3)*(1./2.*c020(i1,i2,i3+1)+1./4.*c010(i1,i2,i3+1));
                                            cp( 1, 1, 1)=1./4.*c200(i1,i2,i3)*c011(i1+1,i2,i3)+1./4.*c020(i1,i2,i3)*c101(i1,i2+1,i3)+1./4.*c002(i1,i2,i3)*c110(i1,i2,i3+1)+c110(i1,i2,i3)*(1./8.*c001(i1+1,i2+1,i3)+1./4.*c002(i1+1,i2+1,i3))+c101(i1,i2,i3)*(1./4.*c020(i1+1,i2,i3+1)+1./8.*c010(i1+1,i2,i3+1))+c011(i1,i2,i3)*(1./8.*c100(i1,i2+1,i3+1)+1./4.*c200(i1,i2+1,i3+1))+1./8.*c100(i1,i2,i3)*c011(i1+1,i2,i3)+1./8.*c010(i1,i2,i3)*c101(i1,i2+1,i3)+1./8.*c001(i1,i2,i3)*c110(i1,i2,i3+1);
                                            cp( 2, 1, 1)=1./16.*c110(i1,i2,i3)*c101(i1+1,i2+1,i3)+1./16.*c101(i1,i2,i3)*c110(i1+1,i2,i3+1);
                                            cp(-2, 2, 1)=0;
                                            cp(-1, 2, 1)=-1./16.*c110(i1,i2,i3)*c011(i1-1,i2+1,i3)-1./16.*c011(i1,i2,i3)*c110(i1,i2+1,i3+1);
                                            cp( 0, 2, 1)=1./4.*c020(i1,i2,i3)*c011(i1,i2+1,i3)+c011(i1,i2,i3)*(1./4.*c020(i1,i2+1,i3+1)+1./8.*c010(i1,i2+1,i3+1))+1./8.*c010(i1,i2,i3)*c011(i1,i2+1,i3);
                                            cp( 1, 2, 1)=1./16.*c110(i1,i2,i3)*c011(i1+1,i2+1,i3)+1./16.*c011(i1,i2,i3)*c110(i1,i2+1,i3+1);
                                            cp( 2, 2, 1)=0;
                                            cp(-2,-2, 2)=0;
                                            cp(-1,-2, 2)=0;
                                            cp( 0,-2, 2)=1./16.*c011(i1,i2,i3)*c011(i1,i2-1,i3+1);
                                            cp( 1,-2, 2)=0;
                                            cp( 2,-2, 2)=0;
                                            cp(-2,-1, 2)=0;
                                            cp(-1,-1, 2)=1./16.*c101(i1,i2,i3)*c011(i1-1,i2,i3+1)+1./16.*c011(i1,i2,i3)*c101(i1,i2-1,i3+1);
                                            cp( 0,-1, 2)=-1./4.*c002(i1,i2,i3)*c011(i1,i2,i3+1)+c011(i1,i2,i3)*(-1./8.*c001(i1,i2-1,i3+1)-1./4.*c002(i1,i2-1,i3+1))-1./8.*c001(i1,i2,i3)*c011(i1,i2,i3+1);
                                            cp( 1,-1, 2)=-1./16.*c101(i1,i2,i3)*c011(i1+1,i2,i3+1)-1./16.*c011(i1,i2,i3)*c101(i1,i2-1,i3+1);
                                            cp( 2,-1, 2)=0;
                                            cp(-2, 0, 2)=1./16.*c101(i1,i2,i3)*c101(i1-1,i2,i3+1);
                                            cp(-1, 0, 2)=-1./4.*c002(i1,i2,i3)*c101(i1,i2,i3+1)+c101(i1,i2,i3)*(-1./8.*c001(i1-1,i2,i3+1)-1./4.*c002(i1-1,i2,i3+1))-1./8.*c001(i1,i2,i3)*c101(i1,i2,i3+1);
                                            cp( 0, 0, 2)=c002(i1,i2,i3)*(1./2.*c001(i1,i2,i3+1)+c002(i1,i2,i3+1))+c101(i1,i2,i3)*(-1./16.*c101(i1+1,i2,i3+1)-1./16.*c101(i1-1,i2,i3+1))+c011(i1,i2,i3)*(-1./16.*c011(i1,i2-1,i3+1)-1./16.*c011(i1,i2+1,i3+1))+c001(i1,i2,i3)*(1./4.*c001(i1,i2,i3+1)+1./2.*c002(i1,i2,i3+1));
                                            cp( 1, 0, 2)=1./4.*c002(i1,i2,i3)*c101(i1,i2,i3+1)+c101(i1,i2,i3)*(1./8.*c001(i1+1,i2,i3+1)+1./4.*c002(i1+1,i2,i3+1))+1./8.*c001(i1,i2,i3)*c101(i1,i2,i3+1);
                                            cp( 2, 0, 2)=1./16.*c101(i1,i2,i3)*c101(i1+1,i2,i3+1);
                                            cp(-2, 1, 2)=0;
                                            cp(-1, 1, 2)=-1./16.*c101(i1,i2,i3)*c011(i1-1,i2,i3+1)-1./16.*c011(i1,i2,i3)*c101(i1,i2+1,i3+1);
                                            cp( 0, 1, 2)=1./4.*c002(i1,i2,i3)*c011(i1,i2,i3+1)+c011(i1,i2,i3)*(1./4.*c002(i1,i2+1,i3+1)+1./8.*c001(i1,i2+1,i3+1))+1./8.*c001(i1,i2,i3)*c011(i1,i2,i3+1);
                                            cp( 1, 1, 2)=1./16.*c101(i1,i2,i3)*c011(i1+1,i2,i3+1)+1./16.*c011(i1,i2,i3)*c101(i1,i2+1,i3+1);
                                            cp( 2, 1, 2)=0;
                                            cp(-2, 2, 2)=0;
                                            cp(-1, 2, 2)=0;
                                            cp( 0, 2, 2)=1./16.*c011(i1,i2,i3)*c011(i1,i2+1,i3+1);
                                            cp( 1, 2, 2)=0;
                                            cp( 2, 2, 2)=0;
                                            ForStencil(m1,m2,m3)
                                            {
                                                int m  = M123(m1,m2,m3);
                                                coeffLocal(m,i1,i2,i3) += cLapSq*cp(m1,m2,m3);
                                            }
                      // // save for testing : 
                      // ForStencil(m1,m2,m3)
                      // {
                      //   int m  = M123(m1,m2,m3);
                      //   lapSqCoeffLocal(m,i1,i2,i3) = cp(m1,m2,m3);
                      // }          
                                        } // end if mask
                                    } // end for 3d 
                                } // end if ok 
                            } // end 3d 
                        } // end curvilinear

                }

            }
            else
            {
        // ----- this grid is advanced with EXPLICIT time-stepping ----
        // set the matrix the IDENTITY
                if( debug & 1 )
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


      // ------------------------------
      // ----- BOUNDARY CONDITIONS ----
      // ------------------------------
            IntegerArray bcLocal(2,3);
            ParallelGridUtility::getLocalBoundaryConditions( mg, bcLocal );

            const int extrapOrder = orderOfAccuracy+1;
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
      // --- FILL BOUNDARY CONDITIONS ----
            ForBoundary(side,axis)
            {

                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Ib1,Ib2,Ib3);

                if( !ok ) continue;

        // // limit extrapolation width if there are not enough grid points *wdh* Dec 4, 2023
        // const int maxExtrapWidth = mg.gridIndexRange(1,axis)-mg.gridIndexRange(0,axis)+1;
        // const int extrapWidth = min(orderOfAccuracy+1,maxExtrapWidth);
        // setExtrapolationWidth(extrapWidth);
        // if( extrapWidth<orderOfAccuracy+1 )
        // {
        //   printF("impCoeff: Limiting extrapolationWidth to %d for (grid,side,axis)=(%d,%d,%d) since the grid is so coarse.\n",extrapWidth,grid,side,axis);
        // }


        // Set the index-shift for this side
                is1=is2=is3=0;
                isv[axis]=1-2*side;   // +1 on left and -1 on right  

                const int axisp1 = (axis+1) % numberOfDimensions;    

                if( mg.boundaryCondition(side,axis)==dirichlet         ||
            // mg.boundaryCondition(side,axis)==CgWave::absorbing ||  // ** DO THIS FOR NOW : absorbing terminated with Dirichlet
                        mg.boundaryCondition(side,axis)==exactBC )
                {
          // ------------ FILL DIRICHLET BC ------------
                        if( bcLocal(side,axis)==exactBC )
                        {
                            if( debug & 1 )
                                printF("+++++ IMPLICIT BC: FILL MATRIX BC=%d FOR (grid,side,axis)=(%d,%d,%d) EXACT BC\n",bcLocal(side,axis),grid,side,axis);
              // ---- exact BC ----
              // Set a Dirichlet condition on the boundary and ghost points 
                            getBoundaryIndex(mg.gridIndexRange(),side,axis,Jb1,Jb2,Jb3,numberOfGhostLines); // extended boundary index
                            if( side==0 )
                                Jbv[axis] = Range(Jbv[axis].getBase()-numberOfGhostLines,Jbv[axis].getBound());  // assign ghost on left
                            else
                                Jbv[axis] = Range(Jbv[axis].getBase(),Jbv[axis].getBound()+numberOfGhostLines);  // asign ghost on right
                            bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Jb1,Jb2,Jb3);
                            if( ok )
                            {
                                const IntegerArray & gid = mg.gridIndexRange();
                                FOR_3D(i1,i2,i3,Jb1,Jb2,Jb3)
                                {
                                    if( maskLocal(i1,i2,i3) > 0 ) // *wdh* avoid changing an "active" point such as interpolation
                                    {
                                        bool isGhost = i1<gid(0,0) || i1>gid(1,0) ||
                                                                      i2<gid(0,1) || i2>gid(1,1) ||
                                                                      i3<gid(0,2) || i3>gid(1,2);
                                        if( isGhost )
                                        {
                                            setClassify(e,i1,i2,i3, SparseRepForMGF::ghost1);  // this is a "real" equation with a RHS (e.g. not extrapolation)
                                        }
                                        coeffLocal(    M,i1,i2,i3) = 0.0;  // zero out any existing equations
                                        coeffLocal(mDiag,i1,i2,i3) = 1.0;
                                    }
                                }  
                            } // end if ok   
                        }
                        else
                        {
                            if( debug & 1 )
                                printF("+++++ IMPLICIT BC: FILL MATRIX BC=%d FOR (grid,side,axis)=(%d,%d,%d) DIRICHLET/ABSORBING\n",bcLocal(side,axis),grid,side,axis);
                            bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Ib1,Ib2,Ib3);
                            if( ok )
                                {
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3)
                                {
                                    if( maskLocal(i1,i2,i3) > 0 ) // *wdh* avoid changing an "active" point such as interpolation
                                    {
                                        coeffLocal(    M,i1,i2,i3) = 0.0;  // zero out any existing equations
                                        coeffLocal(mDiag,i1,i2,i3) = 1.0;
                                    }
                                }
                            }
                        }
            // coeffLocal(    M,Ib1,Ib2,Ib3) = 0.0;  // zero out any existing equations
            // coeffLocal(mDiag,Ib1,Ib2,Ib3) = 1.0;          
                        const bool useCompatibility = orderOfAccuracy==4 && bcApproach==useCompatibilityBoundaryConditions;
                        if( useCompatibility && bcLocal(side,axis)==dirichlet)
                        {
                            if( debug & 1 )
                                printF("+++++ IMPLICIT Dirichlet BC: FILL MATRIX COMPATIBILITY BC (grid,side,axis)=(%d,%d,%d)\n",grid,side,axis);
              // ----- Set a Dirichlet BC on the extended boundary points ----
              //
              //           G--G--B--I--I--I-- ... --I--B--G--G
              //           D  D                           D  D
              // D = set Dirichlet BC here
              //
                            if( true ) // *wdh* Dec 4, 2023
                            {
                                getBoundaryIndex(mg.gridIndexRange(),side,axis,Jb1,Jb2,Jb3,numberOfGhostLines); // extended boundary index
                                if( false )
                                {
                                    printf("IMP:CBC: Set Dirichlet BC on extended boundary (grid,side,axis)=(%d,%d,%d)\n",grid,side,axis); 
                                    printf("fill matrix Dirichlet: Jb1=[%4d,%4d], Jb2=[%4d,%4d] Jb3=[%4d,%4d]\n",
                                            Jb1.getBase(),Jb1.getBound(),
                                            Jb2.getBase(),Jb2.getBound(),
                                            Jb3.getBase(),Jb3.getBound()
                                            );
                                }
                                bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Jb1,Jb2,Jb3);
                                if( ok )
                                {
                                    const IntegerArray & gid = mg.gridIndexRange();
                                    FOR_3D(i1,i2,i3,Jb1,Jb2,Jb3) 
                                    {
                                        if( maskLocal(i1,i2,i3) > 0 ) // *wdh* avoid changing an "active" point such as interpolation 
                                        {       
                                            bool isGhost = i1<gid(0,0) || i1>gid(1,0) ||
                                                                          i2<gid(0,1) || i2>gid(1,1) ||
                                                                          i3<gid(0,2) || i3>gid(1,2);
                                            if( isGhost )
                                            {
                        // printf("IMP:CBC: Set Dirichlet BC on extended boundary (grid,side,axis)=(%d,%d,%d) (i1,i2,i3)=(%d,%d,%d)\n",
                        //          grid,side,axis,i1,i2,i3);
                                                setClassify(e,i1,i2,i3, SparseRepForMGF::ghost1);  
                                                coeffLocal(M,i1,i2,i3) = 0.0;  // zero out any existing equations
                                                ForStencil(m1,m2,m3)
                                                {
                                                    int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                                    if( m1==0 && m2==0 && m3==0 )
                                                        coeffLocal(m,i1,i2,i3) = 1.0;  // Dirichlet BC
                                                    else
                                                        coeffLocal(m,i1,i2,i3) = 0.;
                                                    int j1=i1+m1, j2=i2+m2, j3=i3+m3;                   
                                                    setEquationNumber(m, e,i1,i2,i3,  cc,j1,j2,j3 );    
                                                }
                                            }
                                        }
                                    }
                                } // end if ok
                            }
              // At a Dirichlet-Dirichlet Corner we cannot use a CBCs on both sides since this is a duplicate
              //         |
              //      C--+
              //         |
              //      D--X--+--+--+--
              //         |  |  |
              //         D  C  C
              // D = use dirichlet BC on extended ghost line(s)
              // C = use CBC
              // 
              // IntegerArray skipCorner(2);
                            for( int dira=1; dira<numberOfDimensions; dira++ ) // loop over adjacent directions
                            {
                                int axisa = (axis+dira) % numberOfDimensions;  // adjacent axis
                                for( int sidea=0; sidea<=1; sidea++ )          // adjacent side
                                {
                                    if( bcLocal(sidea,axisa)==dirichlet ||
                                            bcLocal(sidea,axisa)==exactBC )
                                    {
                    // skipCorner(sidea)=1;
                                        Ibv[axisa] = sidea== 0 ? Range(Ibv[axisa].getBase()+1,Ibv[axisa].getBound()  ) 
                                                                                      : Range(Ibv[axisa].getBase()  ,Ibv[axisa].getBound()-1);
                    // Now done above: 
                    // getBoundaryIndex(mg.gridIndexRange(),side,axis,Jb1,Jb2,Jb3);
                    // Jbv[axisa] = mg.gridIndexRange(sidea,axisa);
                    // FOR_3D(i1,i2,i3,Jb1,Jb2,Jb3) 
                    // {
                    //   // for( int ghost=1; ghost<=numberOfGhostLines; ghost++ )
                    //   for( int ghost=1; ghost<=numberOfGhostLines; ghost++ )
                    //   {
                    //     int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                    //     // printf("IMP:CBC: Set Dirichlet BC at ghost=%d adjacent to a corner: (i1m,i2m,i3m)=(%d,%d,%d)\n",
                    //     //         ghost,i1m,i2m,i3m);
                    //     setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost1);  
                    //     ForStencil(m1,m2,m3)
                    //     {
                    //       int m  = M123(m1,m2,m3);        // the single-component coeff-index
                    //       if( m1==0 && m2==0 && m3==0 )
                    //         coeffLocal(m,i1m,i2m,i3m) = 1.0;  // Dirichlet BC
                    //       else
                    //         coeffLocal(m,i1m,i2m,i3m) = 0.;
                    //       int j1=i1m+m1, j2=i2m+m2, j3=i3m+m3;                   
                    //       setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );    
                    //     }
                    //   }
                    // }
                                    }
                                }
                            }
                            bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Ib1,Ib2,Ib3);
                            if( ok )
                            {
                                realSerialArray lapCoeff(M0,Ib1,Ib2,Ib3);
                // Fourth-order laplacian: 
                                mgop.assignCoefficients(MappedGridOperators::laplacianOperator ,lapCoeff, Ib1,Ib2,Ib3,0,0);
                                if( false )
                                {
                                    printf("Fill matrix Dirichlet CBC1: Ib1=[%d,%d] Ib2=[%d,%d] Ib3=[%d,%d]\n",
                                              Ib1.getBase(),Ib1.getBound(),
                                              Ib2.getBase(),Ib2.getBound(),
                                              Ib3.getBase(),Ib3.getBound()
                                              );
                                }
                                const Real cSq = c*c;
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    if( maskLocal(i1,i2,i3) > 0 ) 
                                    {
                                        int i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
                    // Specify that this a "real" equation on the first ghost line: 
                    // (A "real" equation has a possible non-zero right-hand-side)
                                        setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost1);              
                                        ForStencil(m1,m2,m3)
                                        {
                                            int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                            coeffLocal(m,i1m,i2m,i3m) = cSq*lapCoeff(m,i1,i2,i3);
                      // Specify that the above coeff value is the coefficient of component cc at the grid point (j1,j2,j3).
                                            int j1=i1+m1, j2=i2+m2, j3=i3+m3;                       // the stencil is centred on the boundary pt (i1,i2,i3)
                                            setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                        }
                                    }
                                } // end FOR_3D
                // fill ghost 2 with extrapolation
                                if( true )
                                {
                                    for( int ghost=2; ghost<=numberOfGhostLines; ghost++ )
                                    {
                      // *wdh* Nov 22, 2023: 
                      // getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                      // coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0;
                                            bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Ib1,Ib2,Ib3);
                                            if( ok )
                                            {
                                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                                {
                                                    if( maskLocal(i1,i2,i3) != 0 )
                                                    { 
                                                        int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                            // --- fill in the coefficients of the extrapolation formula ---
                                                        coeffLocal(M0,i1m,i2m,i3m)=0; // zero out all *wdh* Nov 22, 2023: 
                                                        for( int m=0; m<=extrapOrder; m++ )
                                                        {
                                                            coeffLocal(m,i1m,i2m,i3m) = extrapCoeff[m];
                                                            int j1=i1m + m*is1, j2=i2m + m*is2, j3=i3m + m*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                            setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                                        }
                                                    }
                                                } // end FOR_3D
                                            }
                                    } // end for ghost
                                }
                            } // end if ok
              // OV_ABORT("finish me");
                        }
                        else if( bcLocal(side,axis)!=exactBC )
                        {
              // --- EXTRAPOLATE GHOST LINES ---
                            for( int ghost=1; ghost<=numberOfGhostLines; ghost++ )
                            {
                  // *wdh* Nov 22, 2023: 
                  // getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                  // coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0;
                                    bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Ib1,Ib2,Ib3);
                                    if( ok )
                                    {
                                        FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                        {
                                            if( maskLocal(i1,i2,i3) != 0 )
                                            { 
                                                int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                        // --- fill in the coefficients of the extrapolation formula ---
                                                coeffLocal(M0,i1m,i2m,i3m)=0; // zero out all *wdh* Nov 22, 2023: 
                                                for( int m=0; m<=extrapOrder; m++ )
                                                {
                                                    coeffLocal(m,i1m,i2m,i3m) = extrapCoeff[m];
                                                    int j1=i1m + m*is1, j2=i2m + m*is2, j3=i3m + m*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                    setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                                }
                                            }
                                        } // end FOR_3D
                                    }
                            } // end for ghost
                        }

                }
                else if( mg.boundaryCondition(side,axis)==neumann )
                {
          // ------------ FILL NEUMANN BC ------------
                        if( debug & 1 )
                            printF("+++++ IMPLICIT BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) NEUMANN\n",grid,side,axis);
                        mg.update(MappedGrid::THEvertexBoundaryNormal);
                        OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal); 
            // TEST
            // getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                        realSerialArray xCoeff(M0,Ib1,Ib2,Ib3), yCoeff(M0,Ib1,Ib2,Ib3), zCoeff; 
                        mgop.assignCoefficients(MappedGridOperators::xDerivative ,xCoeff, Ib1,Ib2,Ib3,0,0);
                        mgop.assignCoefficients(MappedGridOperators::yDerivative ,yCoeff, Ib1,Ib2,Ib3,0,0);
                        if( numberOfDimensions==3 )
                        {
                            zCoeff.redim(M0,Ib1,Ib2,Ib3);
                            mgop.assignCoefficients(MappedGridOperators::zDerivative ,zCoeff, Ib1,Ib2,Ib3,0,0);
                        }
                        const bool useCompatibility = orderOfAccuracy==4 && bcApproach==useCompatibilityBoundaryConditions;
            // [Jb1,Jb2,Jb3] : boundary points, but exclude ends with Dirichlet BC's
                        Jb1=Ib1; Jb2=Ib2; Jb3=Ib3;
                        if( useCompatibility )
                        {
              // exclude end points with Dirichlet BC's
                            for( int dira=1; dira<numberOfDimensions; dira++ )
                            {
                                int axisa = ( axis + dira ) % numberOfDimensions; // adjacent axis
                                int ia = mg.gridIndexRange(0,axisa);
                                int ib = mg.gridIndexRange(1,axisa);
                                if( bcLocal(0,axisa)==dirichlet ) 
                                    ia++;
                                if( bcLocal(1,axisa)==dirichlet ) 
                                    ib--;
                                Jbv[axisa] = Range(ia,ib);
                            }
                        }
                        printf("fill matrix Neumann: Jb1=[%4d,%4d], Jb2=[%4d,%4d] Jb3=[%4d,%4d]\n",
                                        Jb1.getBase(),Jb1.getBound(),
                                        Jb2.getBase(),Jb2.getBound(),
                                        Jb3.getBase(),Jb3.getBound()
                                        );
                        FOR_3D(i1,i2,i3,Jb1,Jb2,Jb3) // loop over points on the boundary
                        {
                            if( maskLocal(i1,i2,i3)>0 )
                            {      
                                int i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
                // Specify that this a "real" equation on the first ghost line: 
                // (A "real" equation has a possible non-zero right-hand-side)
                                setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost1);              
                                coeffLocal(M,i1m,i2m,i3m) = 0.0;  // zero out any existing equations
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
                            }
                        } // end FOR_3D
                        if( useCompatibility )
                        {
                            if( numberOfDimensions==3 )
                            {
                                printF("+++++ IMPLICIT Neumann BC in 3D --- FINISH ME FOR 3D +++++++\n");
                                printF("There  is a potential problem in 3D with a zero pivot from the CBC2 for a box domain\n");
                // OV_ABORT("error");
                            }
                            Range R5 = Range(-2,2);
                            const Real cSq = c*c;
              // if( false && isRectangular ) // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMP %%%%%%%%%%%%%%%%%%%
                            if( isRectangular ) 
                            {
                                if( debug & 1 )
                                    printF("+++++ IMPLICIT Neumann BC: FILL MATRIX COMPATIBILITY BC CARTESIAN (grid,side,axis)=(%d,%d,%d)\n",grid,side,axis);
                                RealArray xxCoeff(R5), yyCoeff(R5), zzCoeff(R5);
                                xxCoeff=0.; yyCoeff=0.; zzCoeff=0.;
                // D+D-
                                xxCoeff(-1) = 1./(dx[0]*dx[0]); yyCoeff(-1) = 1./(dx[1]*dx[1]); zzCoeff(-1) = 1./(dx[2]*dx[2]); 
                                xxCoeff( 0) =-2./(dx[0]*dx[0]); yyCoeff( 0) =-2./(dx[1]*dx[1]); zzCoeff( 0) =-2./(dx[2]*dx[2]); 
                                xxCoeff( 1) = 1./(dx[0]*dx[0]); yyCoeff( 1) = 1./(dx[1]*dx[1]); zzCoeff( 1) = 1./(dx[2]*dx[2]); 
                                RealArray xCoeff(R5), yCoeff(R5), zCoeff(R5);
                                xCoeff=0.; yCoeff=0.; zCoeff=0.;
                // Dz stencil 
                                xCoeff(-1) = -1./(2.*dx[0]); yCoeff(-1) = -1./(2.*dx[1]); zCoeff(-1) = -1./(2.*dx[2]); 
                                xCoeff( 0) =  0.;            yCoeff( 0) =  0.;            zCoeff( 0) =  0.;
                                xCoeff( 1) =  1./(2.*dx[0]); yCoeff( 1) =  1./(2.*dx[1]); zCoeff( 1) =  1./(2.*dx[2]); 
                // Identity stencil 
                                RealArray iCoeff(R5);
                                iCoeff=0.;
                                iCoeff(0)=1.; 
                                RealArray xxxCoeff(R5), yyyCoeff(R5), zzzCoeff(R5);
                                xxxCoeff=0.; yyyCoeff=0.; zzzCoeff=0.;
                // D0(D+D-) stencil 
                                xxxCoeff(-2) = -1./(2.*dx[0]*dx[0]*dx[0]); 
                                xxxCoeff(-1) = +2./(2.*dx[0]*dx[0]*dx[0]); 
                                xxxCoeff( 0) =  0.;                        
                                xxxCoeff( 1) = -2./(2.*dx[0]*dx[0]*dx[0]); 
                                xxxCoeff( 2) =  1./(2.*dx[0]*dx[0]*dx[0]); 
                                yyyCoeff(-2) = -1./(2.*dx[1]*dx[1]*dx[1]); 
                                yyyCoeff(-1) = +2./(2.*dx[1]*dx[1]*dx[1]); 
                                yyyCoeff( 0) =  0.;                        
                                yyyCoeff( 1) = -2./(2.*dx[1]*dx[1]*dx[1]); 
                                yyyCoeff( 2) =  1./(2.*dx[1]*dx[1]*dx[1]); 
                                zzzCoeff(-2) = -1./(2.*dx[2]*dx[2]*dx[2]); 
                                zzzCoeff(-1) = +2./(2.*dx[2]*dx[2]*dx[2]); 
                                zzzCoeff( 0) =  0.;                        
                                zzzCoeff( 1) = -2./(2.*dx[2]*dx[2]*dx[2]); 
                                zzzCoeff( 2) =  1./(2.*dx[2]*dx[2]*dx[2]); 
                // Outward normal is (an1,an2,an3)
                                const Real an1 = axis==0 ? -1.*(1-2*side) : 0.; 
                                const Real an2 = axis==1 ? -1.*(1-2*side) : 0.; 
                                const Real an3 = axis==2 ? -1.*(1-2*side) : 0.; 
                                printF("Neuman CBC Delta(u.n) : (an1,an2,an3)=(%12.4e,%12.4e,%12.4e), is1=%d, is2=%d, is3=%d\n",an1,an2,an3,is1,is2,is3);
                                FOR_3D(i1,i2,i3,Jb1,Jb2,Jb3) // loop over points on the boundary
                                {
                                    if( maskLocal(i1,i2,i3)>0 )
                                    {
                                        int i1m=i1-2*is1, i2m=i2-2*is2, i3m=i3-2*is3; //  2nd ghost point is (i1m,i2m,i3m)
                    // Specify that this a "real" equation on the first ghost line: 
                    // (A "real" equation has a possible non-zero right-hand-side)
                                        setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost1);              
                                        coeffLocal(M,i1m,i2m,i3m) = 0.0;  // zero out any existing equations
                                        ForStencil(m1,m2,m3)
                                        {
                                            int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                            if( false ) // try this **TEMP***
                                            {
                        // This mimics replacing cross terms where we can:
                        // If
                        //     u.x = g
                        // Then
                        //     u.xyy = g.yy, u.xzz = g.zz
                        //  
                        // This seems to work
                                                if( axis==0 )
                                                    coeffLocal(m,i1m,i2m,i3m) = an1*cSq*(  xxxCoeff(m1)*  iCoeff(m2)*  iCoeff(m3) );
                                                else if( axis==1 )
                                                    coeffLocal(m,i1m,i2m,i3m) = an2*cSq*(  iCoeff(m1)  *yyyCoeff(m2)*  iCoeff(m3) );
                                                else 
                                                    coeffLocal(m,i1m,i2m,i3m) = an3*cSq*(  iCoeff(m1)  *  iCoeff(m2)*zzzCoeff(m3) );              
                                            }
                                            else
                                            {
                                                if( axis==0 )
                                                { // left or right side: u_xxx + u_xyy + u_xzz
                                                    coeffLocal(m,i1m,i2m,i3m) = 
                                                                                an1*cSq*(  xxxCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                                                                                                    +  xCoeff(m1)*  yyCoeff(m2)*   iCoeff(m3) 
                                                                                              );
                                                    if( numberOfDimensions==3 )
                                                    {
                                                        coeffLocal(m,i1m,i2m,i3m) += an1*cSq*(  xCoeff(m1)*   iCoeff(m2)*   zzCoeff(m3) );
                                                    }
                                                }
                                                else if( axis==1 )
                                                { // top or bottom: u_xxy + u_yyy + u_yzz
                                                    coeffLocal(m,i1m,i2m,i3m) = 
                                                                                  an2*cSq*(  xxCoeff(m1)*  yCoeff(m2)*   iCoeff(m3) 
                                                                                                    +  iCoeff(m1)*yyyCoeff(m2)*   iCoeff(m3) 
                                                                                                  ); 
                                                    if( numberOfDimensions==3 )
                                                    {
                                                        coeffLocal(m,i1m,i2m,i3m) += an2*cSq*(  iCoeff(m1)*   yCoeff(m2)*   zzCoeff(m3) );
                                                    }                                      
                                                } 
                                                else
                                                { // front or back: u_xxz + u_yyz + u_zzz
                                                    coeffLocal(m,i1m,i2m,i3m) = 
                                                                                an3*cSq*(  xxCoeff(m1)*   iCoeff(m2)*   zCoeff(m3) 
                                                                                                    + iCoeff(m1)*  yyCoeff(m2)*   zCoeff(m3) 
                                                                                                    + iCoeff(m1)*   iCoeff(m2)* zzzCoeff(m3) 
                                                                                                  );                    
                                                }               
                        // if( i2==Jb2.getBase() && i3==Jb3.getBase() && coeffLocal(m,i1m,i2m,i3m) != 0. )
                        // {
                        //   printF("(i1m,i2m,i3m)=(%4d,%4d,%4d) m=%d coeff=%10.2e\n",i1m,i2m,i3m,m,coeffLocal(m,i1m,i2m,i3m));
                        // }
                                            }
                      // Specify that the above coeff value is the coefficient of component cc at the grid point (j1,j2,j3).
                                            int j1=i1+m1, j2=i2+m2, j3=i3+m3;                       // the stencil is centred on the boundary pt (i1,i2,i3)
                                            setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                        }
                                    }
                                } // end FOR_3D
                            }
                            else
                            {
                                printF("FILL MATRIX WITH Neumann CBC2 order=4 curvilinear (grid,side,axis)=(%d,%d,%d)\n",grid,side,axis);
                                if( isRectangular ) // %%%%%%%%%%%%%%%%%%%%%%%% TEMP %%%%%%%%%%%%%%%%%%%%
                                { // rectangular grid grid-spacings: 
                                    mg.update(MappedGrid::THEinverseVertexDerivative );
                  // unit square grid spacings: 
                                    for( int dir=0; dir<3; dir++ )
                                    {
                                        dr[dir]=mg.gridSpacing(dir);   
                                    }
                                }   
                                OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rxLocal);
                // macro to make the rxLocal array look 5-dimensional 
                                #define DD(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2))     
                // ----- COMPUTE DERIVATIVES OF METRICS -----
                                Range Rd2=numberOfDimensions*numberOfDimensions;
                                RealArray ddx(Ib1,Ib2,Ib3,Rd2), ddy(Ib1,Ib2,Ib3,Rd2);
                                mgop.derivative( MappedGridOperators::xDerivative,rxLocal,ddx,Ib1,Ib2,Ib3,Rd2);            
                                mgop.derivative( MappedGridOperators::yDerivative,rxLocal,ddy,Ib1,Ib2,Ib3,Rd2);
                                RealArray ddxx(Ib1,Ib2,Ib3,Rd2), ddxy(Ib1,Ib2,Ib3,Rd2), ddyy(Ib1,Ib2,Ib3,Rd2);    
                                mgop.derivative( MappedGridOperators::xxDerivative,rxLocal,ddxx,Ib1,Ib2,Ib3,Rd2);            
                                mgop.derivative( MappedGridOperators::xyDerivative,rxLocal,ddxy,Ib1,Ib2,Ib3,Rd2);            
                                mgop.derivative( MappedGridOperators::yyDerivative,rxLocal,ddyy,Ib1,Ib2,Ib3,Rd2);             
                                #define DDX(i1,i2,i3,m1,m2) ddx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                                #define DDY(i1,i2,i3,m1,m2) ddy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 
                                #define DDXX(i1,i2,i3,m1,m2) ddxx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                                #define DDXY(i1,i2,i3,m1,m2) ddxy(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                                #define DDYY(i1,i2,i3,m1,m2) ddyy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 
                // Define stencils for parametric derivatives 
                                Range R5 = Range(-2,2);
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
                                if( numberOfDimensions==2 )
                                {      
                                    FOR_3D(i1,i2,i3,Jb1,Jb2,Jb3) // loop over points on the boundary
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
                      // Real rxxxx = DDXXX(i1,i2,i3,0,0);  // rx.xxx
                      // Real rxxyy = DDXYY(i1,i2,i3,0,1);  // ry.xyy
                      // Real ryyyy = DDYYY(i1,i2,i3,0,1);  // ry.yyy   
                                            Real sx    =    DD(i1,i2,i3,1,0);
                                            Real sy    =    DD(i1,i2,i3,1,1);
                                            Real sxx   =   DDX(i1,i2,i3,1,0);
                                            Real sxy   =   DDY(i1,i2,i3,1,0);
                                            Real syy   =   DDY(i1,i2,i3,1,1);  
                                            Real sxxx  =  DDXX(i1,i2,i3,1,0);
                                            Real sxxy  =  DDXY(i1,i2,i3,1,0);
                                            Real sxyy  =  DDYY(i1,i2,i3,1,0);  // sx.yy
                                            Real syyy  =  DDYY(i1,i2,i3,1,1);  
                      // Real sxxxx = DDXXX(i1,i2,i3,1,0);  // sx.xxx
                      // Real sxxyy = DDXYY(i1,i2,i3,1,1);  // sy.xyy
                      // Real syyyy = DDYYY(i1,i2,i3,1,1);  // sy.yyy             
                                            const Real an1 =  normal(i1,i2,i3,0);
                                            const Real an2 =  normal(i1,i2,i3,1);
                      // ---- COEFFICIENTS OF 2D n.Grad( LAPLACIAN ) : from laplacianCoefficients.mpl ----
                                            Real urrr = an1 * (pow(rx, 0.3e1) + ry * ry * rx) + an2 * (ry * rx * rx + pow(ry, 0.3e1));
                                            Real urrs = an1 * (3 * rx * rx * sx + ry * (sy * rx + ry * sx) + sy * ry * rx) + an2 * (sy * rx * rx + 2 * ry * sx * rx + 3 * ry * ry * sy);
                                            Real urss = an1 * (3 * rx * sx * sx + ry * sx * sy + sy * (sy * rx + ry * sx)) + an2 * (2 * sy * sx * rx + ry * sx * sx + 3 * ry * sy * sy);
                                            Real usss = an1 * (pow(sx, 0.3e1) + sx * sy * sy) + an2 * (sy * sx * sx + pow(sy, 0.3e1));
                                            Real urr = an1 * (3 * rx * rxx + ryy * rx + 2 * ry * rxy) + an2 * (2 * rxy * rx + ry * rxx + 3 * ry * ryy);
                                            Real urs = an1 * (3 * sxx * rx + syy * rx + 3 * sx * rxx + 2 * sy * rxy + 2 * ry * sxy + ryy * sx) + an2 * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx + 3 * syy * ry + 3 * sy * ryy);
                                            Real uss = an1 * (3 * sx * sxx + sx * syy + 2 * sxy * sy) + an2 * (2 * sx * sxy + sxx * sy + 3 * sy * syy);
                                            Real ur = an1 * (rxxx + rxyy) + an2 * (rxxy + ryyy);
                                            Real us = an1 * (sxxx + sxyy) + an2 * (sxxy + syyy);
                                            int i1m=i1-2*is1, i2m=i2-2*is2, i3m=i3-2*is3; //  2nd ghost point is (i1m,i2m,i3m)
                      // Specify that this a "real" equation on the first ghost line: 
                      // (A "real" equation has a possible non-zero right-hand-side)
                                            setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost1);              
                      // printF(" (i1,i2)=(%d,%d) coeff: urrr=%10.2e urrs=%10.2e urss=%10.2e usss=%10.2e urr=%10.2e urs=%10.2e uss=%10.2e ur=%10.2e us=%10.2e\n",
                      //      i1,i2,urrr,urrs,urss,usss,urr,urs,uss,ur,us);
                                            coeffLocal(M,i1m,i2m,i3m) = 0.0;  // zero out any existing equations
                                            ForStencil(m1,m2,m3)
                                            {
                                                int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                                coeffLocal(m,i1m,i2m,i3m) = 
                                                                                      cSq*(  urrr * rrrCoeff(m1)*   iCoeff(m2)
                                                                                                + urrs *  rrCoeff(m1)*   sCoeff(m2)
                                                                                                + urss *   rCoeff(m1)*  ssCoeff(m2)
                                                                                                + usss *   iCoeff(m1)* sssCoeff(m2)
                                                                                                + urr  *  rrCoeff(m1)*   iCoeff(m2)
                                                                                                + urs  *   rCoeff(m1)*   sCoeff(m2)
                                                                                                + uss  *   iCoeff(m1)*  ssCoeff(m2)
                                                                                                + ur   *   rCoeff(m1)*   iCoeff(m2)
                                                                                                + us   *   iCoeff(m1)*   sCoeff(m2)
                                                                                                );          
                        // Specify that the above coeff value is the coefficient of component cc at the grid point (j1,j2,j3).
                                                int j1=i1+m1, j2=i2+m2, j3=i3+m3;                       // the stencil is centred on the boundary pt (i1,i2,i3)
                                                setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                            }
                                        } // end if mask
                                    } // end FOR_3D
                                }
                                else
                                {
                  // --- THREE DIMENSIONS -----
                  // OV_ABORT("Neumann CBC order=4 curvilinear: finish me for 3D");
                                    RealArray rxza(I1,I2,I3,Rd2);
                                    mgop.derivative( MappedGridOperators::zDerivative,rxLocal,rxza,I1,I2,I3,Rd2);            
                                    #define DDZ(i1,i2,i3,m1,m2) rxza(i1,i2,i3,(m1)+numberOfDimensions*(m2))
                                    const int extra=1; 
                                    RealArray ddxz(I1,I2,I3,Rd2), ddyz(I1,I2,I3,Rd2), ddzz(I1,I2,I3,Rd2);         
                                    mgop.derivative( MappedGridOperators::xzDerivative,rxLocal,ddxz,I1,I2,I3,Rd2);            
                                    mgop.derivative( MappedGridOperators::yzDerivative,rxLocal,ddyz,I1,I2,I3,Rd2);            
                                    mgop.derivative( MappedGridOperators::zzDerivative,rxLocal,ddzz,I1,I2,I3,Rd2);             
                                    #define DDXZ(i1,i2,i3,m1,m2) ddxz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                                    #define DDYZ(i1,i2,i3,m1,m2) ddyz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                                    #define DDZZ(i1,i2,i3,m1,m2) ddzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                  // RealArray ddxxx(I1,I2,I3,Rd2), ddxyy(I1,I2,I3,Rd2), ddyyy(I1,I2,I3,Rd2);
                  // RealArray ddxzz(I1,I2,I3,Rd2), ddzzz(I1,I2,I3,Rd2), ddyzz(I1,I2,I3,Rd2);
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
                                    FOR_3D(i1,i2,i3,Jb1,Jb2,Jb3) // loop over points on the boundary
                                    {
                                        if( maskLocal(i1,i2,i3)>0 )
                                        {
                      // declareMetricDerivatives3d(r,0)
                      // declareMetricDerivatives3d(s,1)
                      // declareMetricDerivatives3d(t,2)
                                            const Real an1 = normal(i1,i2,i3,0);
                                            const Real an2 = normal(i1,i2,i3,1);
                                            const Real an3 = normal(i1,i2,i3,2);
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
                      // Real rxxxx  =  DDXXX(i1,i2,i3,0,0);  // .xxxx
                      // Real rxxyy  =  DDXYY(i1,i2,i3,0,0);  // .xxyy
                      // Real ryyyy  =  DDYYY(i1,i2,i3,0,1);  // .yyy
                      // Real rxxzz  =  DDXZZ(i1,i2,i3,0,0);  // .xxzz
                      // Real rzzzz  =  DDZZZ(i1,i2,i3,0,2);  // .zzzz
                      // Real ryyzz  =  DDYZZ(i1,i2,i3,0,1);  // .yyzz                        
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
                      // Real sxxxx  =  DDXXX(i1,i2,i3,1,0);  // .xxxx
                      // Real sxxyy  =  DDXYY(i1,i2,i3,1,0);  // .xxyy
                      // Real syyyy  =  DDYYY(i1,i2,i3,1,1);  // .yyy
                      // Real sxxzz  =  DDXZZ(i1,i2,i3,1,0);  // .xxzz
                      // Real szzzz  =  DDZZZ(i1,i2,i3,1,2);  // .zzzz
                      // Real syyzz  =  DDYZZ(i1,i2,i3,1,1);  // .yyzz                        
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
                      // Real txxxx  =  DDXXX(i1,i2,i3,2,0);  // .xxxx
                      // Real txxyy  =  DDXYY(i1,i2,i3,2,0);  // .xxyy
                      // Real tyyyy  =  DDYYY(i1,i2,i3,2,1);  // .yyy
                      // Real txxzz  =  DDXZZ(i1,i2,i3,2,0);  // .xxzz
                      // Real tzzzz  =  DDZZZ(i1,i2,i3,2,2);  // .zzzz
                      // Real tyyzz  =  DDYZZ(i1,i2,i3,2,1);  // .yyzz                
                      // ---- COEFFICIENTS OF 3D n.Grad( LAPLACIAN ): from laplacianCoefficients.mpl ----
                                            Real urrr = an1 * (pow(rx, 0.3e1) + ry * ry * rx) + an2 * (ry * rx * rx + pow(ry, 0.3e1)) + an1 * rz * rz * rx + an2 * rz * rz * ry + an3 * (rz * rx * rx + rz * ry * ry + pow(rz, 0.3e1));
                                            Real urrs = an1 * (3 * rx * rx * sx + ry * (sy * rx + ry * sx) + sy * ry * rx) + an2 * (sy * rx * rx + 2 * ry * sx * rx + 3 * ry * ry * sy) + an1 * (rz * (sz * rx + rz * sx) + sz * rz * rx) + an2 * (rz * (sz * ry + rz * sy) + sz * rz * ry) + an3 * (sz * rx * rx + 2 * rz * sx * rx + sz * ry * ry + 2 * rz * sy * ry + 3 * rz * rz * sz);
                                            Real urss = an1 * (3 * rx * sx * sx + ry * sx * sy + sy * (sy * rx + ry * sx)) + an2 * (2 * sy * sx * rx + ry * sx * sx + 3 * ry * sy * sy) + an1 * (rz * sz * sx + sz * (sz * rx + rz * sx)) + an2 * (rz * sz * sy + sz * (sz * ry + rz * sy)) + an3 * (2 * sz * sx * rx + 2 * sz * sy * ry + rz * sx * sx + rz * sy * sy + 3 * rz * sz * sz);
                                            Real usss = an1 * (pow(sx, 0.3e1) + sx * sy * sy) + an2 * (sy * sx * sx + pow(sy, 0.3e1)) + an1 * sz * sz * sx + an2 * sz * sz * sy + an3 * (sz * sx * sx + sz * sy * sy + pow(sz, 0.3e1));
                                            Real urrt = an1 * (ry * (ty * rx + ry * tx) + ty * ry * rx + 3 * rx * rx * tx) + an2 * (ty * rx * rx + 2 * ry * tx * rx + 3 * ry * ry * ty) + an1 * (rz * (tz * rx + rz * tx) + tz * rz * rx) + an2 * (rz * (tz * ry + rz * ty) + tz * rz * ry) + an3 * (tz * rx * rx + 2 * rz * tx * rx + tz * ry * ry + 2 * rz * ty * ry + 3 * rz * rz * tz);
                                            Real urtt = an1 * (ry * tx * ty + ty * (ty * rx + ry * tx) + 3 * rx * tx * tx) + an2 * (2 * ty * tx * rx + ry * tx * tx + 3 * ry * ty * ty) + an1 * (rz * tx * tz + tz * (tz * rx + rz * tx)) + an2 * (rz * ty * tz + tz * (tz * ry + rz * ty)) + an3 * (2 * tz * tx * rx + 2 * tz * ty * ry + rz * tx * tx + rz * ty * ty + 3 * rz * tz * tz);
                                            Real uttt = an1 * (pow(tx, 0.3e1) + tx * ty * ty) + an2 * (ty * tx * tx + pow(ty, 0.3e1)) + an1 * tx * tz * tz + an2 * ty * tz * tz + an3 * (tz * tx * tx + tz * ty * ty + pow(tz, 0.3e1));
                                            Real usst = an1 * (sy * (ty * sx + sy * tx) + ty * sx * sy + 3 * sx * sx * tx) + an2 * (ty * sx * sx + 2 * sy * tx * sx + 3 * sy * sy * ty) + an1 * (sz * (tz * sx + sz * tx) + tz * sz * sx) + an2 * (sz * (tz * sy + sz * ty) + tz * sz * sy) + an3 * (tz * sx * sx + 2 * sz * tx * sx + tz * sy * sy + 2 * sz * ty * sy + 3 * sz * sz * tz);
                                            Real ustt = an1 * (sy * tx * ty + ty * (ty * sx + sy * tx) + 3 * sx * tx * tx) + an2 * (2 * ty * tx * sx + sy * tx * tx + 3 * sy * ty * ty) + an1 * (sz * tx * tz + tz * (tz * sx + sz * tx)) + an2 * (sz * ty * tz + tz * (tz * sy + sz * ty)) + an3 * (2 * tz * tx * sx + 2 * tz * ty * sy + sz * tx * tx + sz * ty * ty + 3 * sz * tz * tz);
                                            Real urr = an1 * (3 * rx * rxx + ryy * rx + 2 * ry * rxy) + an2 * (2 * rxy * rx + ry * rxx + 3 * ry * ryy) + an1 * (rzz * rx + 2 * rz * rxz) + an2 * (rzz * ry + 2 * rz * ryz) + an3 * (2 * rxz * rx + rz * rxx + 2 * ryz * ry + rz * ryy + 3 * rz * rzz);
                                            Real urs = an1 * (3 * sxx * rx + syy * rx + 3 * sx * rxx + 2 * sy * rxy + 2 * ry * sxy + ryy * sx) + an2 * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx + 3 * syy * ry + 3 * sy * ryy) + an1 * (szz * rx + 2 * sz * rxz + 2 * rz * sxz + rzz * sx) + an2 * (szz * ry + 2 * sz * ryz + 2 * rz * syz + rzz * sy) + an3 * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + 2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * sxx + rz * syy + 3 * rz * szz + 3 * sz * rzz);
                                            Real uss = an1 * (3 * sx * sxx + sx * syy + 2 * sxy * sy) + an2 * (2 * sx * sxy + sxx * sy + 3 * sy * syy) + an1 * (szz * sx + 2 * sz * sxz) + an2 * (szz * sy + 2 * sz * syz) + an3 * (2 * sxz * sx + sz * sxx + 2 * syz * sy + sz * syy + 3 * sz * szz);
                                            Real urt = an1 * (3 * txx * rx + tyy * rx + 3 * tx * rxx + 2 * ty * rxy + 2 * ry * txy + ryy * tx) + an2 * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx + 3 * tyy * ry + 3 * ty * ryy) + an1 * (tzz * rx + 2 * tz * rxz + 2 * rz * txz + rzz * tx) + an2 * (tzz * ry + 2 * tz * ryz + 2 * rz * tyz + rzz * ty) + an3 * (2 * txz * rx + tz * rxx + 2 * tx * rxz + 2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * txx + rz * tyy + 3 * tzz * rz + 3 * tz * rzz);
                                            Real utt = an1 * (3 * tx * txx + tx * tyy + 2 * txy * ty) + an2 * (2 * tx * txy + txx * ty + 3 * ty * tyy) + an1 * (tx * tzz + 2 * txz * tz) + an2 * (ty * tzz + 2 * tyz * tz) + an3 * (2 * tx * txz + txx * tz + 2 * ty * tyz + tyy * tz + 3 * tz * tzz);
                                            Real ust = an1 * (3 * txx * sx + tyy * sx + 3 * tx * sxx + 2 * ty * sxy + 2 * sy * txy + syy * tx) + an2 * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx + 3 * tyy * sy + 3 * ty * syy) + an1 * (tzz * sx + 2 * tz * sxz + 2 * sz * txz + szz * tx) + an2 * (tzz * sy + 2 * tz * syz + 2 * sz * tyz + szz * ty) + an3 * (2 * txz * sx + tz * sxx + 2 * tx * sxz + 2 * tyz * sy + tz * syy + 2 * ty * syz + sz * txx + sz * tyy + 3 * tzz * sz + 3 * tz * szz);
                                            Real ur = an1 * (rxxx + rxyy) + an2 * (rxxy + ryyy) + an1 * rxzz + an2 * ryzz + an3 * (rzzz + ryyz + rxxz);
                                            Real us = an1 * (sxxx + sxyy) + an2 * (sxxy + syyy) + an1 * sxzz + an2 * syzz + an3 * (szzz + syyz + sxxz);
                                            Real ut = an1 * (txyy + txxx) + an2 * (tyyy + txxy) + an1 * txzz + an2 * tyzz + an3 * (tzzz + tyyz + txxz);
                                            int i1m=i1-2*is1, i2m=i2-2*is2, i3m=i3-2*is3; //  2nd ghost point is (i1m,i2m,i3m)
                      // Specify that this a "real" equation on the first ghost line: 
                      // (A "real" equation has a possible non-zero right-hand-side)
                                            setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost1);  
                      // if( i2==Jb2.getBase() && i3==Jb3.getBase()  )
                      // {
                      //   printF(" (i1,i2,i3)=(%3d,%3d,%3d) coeff: urrr=%10.2e urrs=%10.2e urss=%10.2e usss=%10.2e uttt=%10.2e urr=%10.2e urs=%10.2e uss=%10.2e ur=%10.2e us=%10.2e\n",
                      //       i1,i2,i3,urrr,urrs,urss,usss,uttt,urr,urs,uss,ur,us);
                      //   printF("   urrr,urrs,urss,usss,urrt,urtt,uttt,usst,ustt=%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e \n",
                      //   urrr,urrs,urss,usss,urrt,urtt,uttt,usst,ustt);
                      //   printF("   urr,urs,urt,uss,ust,utt,ur,us,ut=%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e \n",urr,urs,urt,uss,ust,utt,ur,us,ut);
                      //   printF(" rxx,rxy,rxz,ryy,ryz,rzz=%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n",rxx,rxy,rxz,ryy,ryz,rzz);
                      //   printF(" sxx,sxy,sxz,syy,syz,szz=%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n",sxx,sxy,sxz,syy,syz,szz);
                      //   printF(" txx,txy,txz,tyy,tyz,tzz=%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n",txx,txy,txz,tyy,tyz,tzz);
                      // }
                                            coeffLocal(M,i1m,i2m,i3m) = 0.0;  // zero out any existing equations
                                            ForStencil(m1,m2,m3)
                                            {
                                                int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                                coeffLocal(m,i1m,i2m,i3m) = 
                                                                                      cSq*( 
                                                                                                    urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
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
                        // Specify that the above coeff value is the coefficient of component cc at the grid point (j1,j2,j3).
                                                int j1=i1+m1, j2=i2+m2, j3=i3+m3;                       // the stencil is centred on the boundary pt (i1,i2,i3)
                                                setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                        // if( i2==Jb2.getBase() && i3==Jb3.getBase() && coeffLocal(m,i1m,i2m,i3m) != 0. )
                        // {
                        //   printF("(i1m,i2m,i3m)=(%4d,%4d,%4d) m=%d coeff=%10.2e\n",i1m,i2m,i3m,m,coeffLocal(m,i1m,i2m,i3m));
                        // }
                                            } // end for stencil
                                        } // end if mask
                                    } // end FOR_3D
                                }
                            }
                        }
                        else
                        {
              // fill ghost 2 with extrapolation
                            for( int ghost=2; ghost<=numberOfGhostLines; ghost++ )
                            {
                  // *wdh* Nov 22, 2023: 
                  // getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                  // coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0;
                                    bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Ib1,Ib2,Ib3);
                                    if( ok )
                                    {
                                        FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                        {
                                            if( maskLocal(i1,i2,i3) != 0 )
                                            { 
                                                int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                        // --- fill in the coefficients of the extrapolation formula ---
                                                coeffLocal(M0,i1m,i2m,i3m)=0; // zero out all *wdh* Nov 22, 2023: 
                                                for( int m=0; m<=extrapOrder; m++ )
                                                {
                                                    coeffLocal(m,i1m,i2m,i3m) = extrapCoeff[m];
                                                    int j1=i1m + m*is1, j2=i2m + m*is2, j3=i3m + m*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                    setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                                }
                                            }
                                        } // end FOR_3D
                                    }
                            } // end for ghost
                        }

                }
                else if( mg.boundaryCondition(side,axis)==CgWave::absorbing || 
                                  mg.boundaryCondition(side,axis)==CgWave::abcEM2 )
                {
                    Real ca = c;
                    if( solveHelmholtz )
                    {
            // Adjust c for the EM2 absorbing BC to account for time-discretization errors
            //   D+t (Dx ) w + A+( ... )            
                        ca = c*tan(frequencyArray(0)*dt/2)/(frequencyArraySave(0)*dt/2); 
                    }

                    if( debug & 1 )
                    {
                        if( mg.boundaryCondition(side,axis)==CgWave::absorbing )
                            printF("+++++ IMPLICIT BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) ABSORBING/EM2, c=%12.4e, cEM2=%12.4e\n",grid,side,axis,c,ca);
                        else if( mg.boundaryCondition(side,axis)==CgWave::abcEM2 )
                            printF("+++++ IMPLICIT BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) EM2, c=%12.4e, cEM2=%12.4e\n",grid,side,axis,c,ca);
                    }


                    if( numberOfDimensions==3 )
                    {
                        OV_ABORT(" IMPLICIT BC: FILL MATRIX -- FINISH ME in 3D"); 
                    }

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

          // res = -is*(unx-ucx)/dt + (.5*c)*( unxx + ucxx) + (.25*c)*( unyy + ucyy );

          // printF("fill Implicit with EM2 BC -- c=%12.4e, ca=%12.4e dt=%12.4e\n",c,ca,dt );

                    Range Rw(-halfWidth1,halfWidth1);
                    RealArray abcCoeff(Rw,Rw,Rw), abcCoeff2(Rw,Rw,Rw);
                    abcCoeff=0.; abcCoeff2=0.; 

                    bool useCompatiblityRBC= orderOfAccuracy>=4; // if true apply a CBC for the 2nd ghost line on an EM2 boundary



          // Engquist-Majda order2 scheme: 
          //  -is*D+t ( D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = f 
          //  -is*D+t ( D0x W^n) + L*( .5 W^{n+1} + .5*W^n ) = f
          //  -is ( D0x W^{n+1} - W^n )/dt + L*( .5 W^{n+1} + .5*W^n ) = f 
          //  -is*D0x W^{n+1}/dt + .5*L W^{n+1} = -is*D0xW^n/dt - .5*L W^n + f 
                    if( isRectangular )
                    {

                        if( orderOfAccuracy==2 ) 
                        {
                            const Real cxx = axis==0 ? ca/(dx[0]*dx[0]) : .5*ca/(dx[0]*dx[0]);
                            const Real cyy = axis==1 ? ca/(dx[1]*dx[1]) : .5*ca/(dx[1]*dx[1]);
                            const Real cx  = axis==0 ? -is1*1./(2.*dx[0]*dt) : 0.;                
                            const Real cy  = axis==1 ? -is2*1./(2.*dx[1]*dt) : 0.;

                            abcCoeff( 0,-1,0) = .5*(         1.*cyy )    - cy;
                            abcCoeff(-1, 0,0) = .5*( 1.*cxx         )    - cx;
                            abcCoeff( 0, 0,0) = .5*(-2.*cxx -2.*cyy );
                            abcCoeff(+1, 0,0) = .5*( 1.*cxx         )    + cx;
                            abcCoeff( 0,+1,0) = .5*(         1.*cyy )    + cy;

                        }
                        else if( orderOfAccuracy==4 )
                        {

              // -is * D+t D4x W^n + c*(Dxx4+.5*Dyy4) (W^{n+1} + W^n)/2

                            if( 1==1 )
                            {
                // exactly match matlab code
                                int is; 
                                const Real dx1=dx[0], dy1=dx[1];
                                const Real dx2=dx1*dx1, dx3=dx2*dx1, dx4=dx3*dx1;
                                const Real dy2=dy1*dy1, dy3=dy2*dy1, dy4=dy3*dy1;  
                                Real cx=0., cxx=0., cy=0., cyy=0.;              
                                if( axis==0 ) 
                                {
                                    is = is1;
                                    cx=1./(12.*dx1); cy=0;            cxx=1./(12.*dx2); cyy=.5/(12.*dy2); 
                                }
                                else
                                {
                                    is = is2; 
                                    cx=0;            cy=1./(12.*dy1); cxx=.5/(12.*dx2); cyy=1./(12.*dy2); 
                                }


                                abcCoeff( 0,-2,0) = (1./(dt))*(      + 1.*cy) -is*.5*ca*(          - 1.*cyy ); 
                                abcCoeff( 0,-1,0) = (1./(dt))*(      - 8.*cy) -is*.5*ca*(          +16.*cyy ); 
                                abcCoeff(-2, 0,0) = (1./(dt))*( 1.*cx       ) -is*.5*ca*(  -1.*cxx          ); 
                                abcCoeff(-1, 0,0) = (1./(dt))*(-8.*cx       ) -is*.5*ca*(  16.*cxx          );
                                abcCoeff( 0, 0,0) = (1./(dt))*( 0.*cx       ) -is*.5*ca*( -30.*cxx -30.*cyy );
                                abcCoeff(+1, 0,0) = (1./(dt))*( 8.*cx       ) -is*.5*ca*(  16.*cxx          );
                                abcCoeff(+2, 0,0) = (1./(dt))*(-1.*cx       ) -is*.5*ca*(  -1.*cxx          ); 
                                abcCoeff( 0,+1,0) = (1./(dt))*(      + 8.*cy) -is*.5*ca*(          +16.*cyy ); 
                                abcCoeff( 0,+2,0) = (1./(dt))*(      - 1.*cy) -is*.5*ca*(          - 1.*cyy ); 

                                if( useCompatiblityRBC )
                                {
                                    Real cxxx=0., cxxxx=0., cyyy=0., cyyyy=0., cxxyy=0.;
                                    if( axis==0 ) 
                                    {
                                        is=is1;
                                        cxxx=1./(2.*dx3 ); cyyy=0;          cxxxx=1./(dx4 ); cyyyy=0;         cxxyy=.5/(dx2 *dy2 ); 
                                    }
                                    else 
                                    {
                                        is=is2;
                                        cxxx=0;            cyyy=1/(2.*dy3 ); cxxxx=0;        cyyyy=1./(dy4 ); cxxyy=.5/(dx2 *dy2 );
                                    }

                                    abcCoeff2(-1,-1,0) =                               -is*.5*ca*(                           cxxyy); 
                                    abcCoeff2(+1,-1,0) =                               -is*.5*ca*(                           cxxyy); 
                                    abcCoeff2( 0,-2,0) = (1./(dt))*(        - 1.*cyyy) -is*.5*ca*(             1.*cyyyy           ); 
                                    abcCoeff2( 0,-1,0) = (1./(dt))*(        + 2.*cyyy) -is*.5*ca*(            -4.*cyyyy  -2.*cxxyy); 
                                    abcCoeff2(-2, 0,0) = (1./(dt))*(-1.*cxxx         ) -is*.5*ca*(   1.*cxxxx                     ); 
                                    abcCoeff2(-1, 0,0) = (1./(dt))*( 2.*cxxx         ) -is*.5*ca*(  -4.*cxxxx            -2.*cxxyy);
                                    abcCoeff2( 0, 0,0) = (1./(dt))*( 0.*cxxx         ) -is*.5*ca*(   6.*cxxxx  +6.*cyyyy +4.*cxxyy);
                                    abcCoeff2(+1, 0,0) = (1./(dt))*(-2.*cxxx         ) -is*.5*ca*(  -4.*cxxxx            -2.*cxxyy);
                                    abcCoeff2(+2, 0,0) = (1./(dt))*( 1.*cxxx         ) -is*.5*ca*(   1.*cxxxx                     ); 
                                    abcCoeff2( 0,+1,0) = (1./(dt))*(        - 2.*cyyy) -is*.5*ca*(            -4.*cyyyy  -2.*cxxyy); 
                                    abcCoeff2( 0,+2,0) = (1./(dt))*(          1.*cyyy) -is*.5*ca*(             1.*cyyyy           );
                                    abcCoeff2(-1,+1,0) =                               -is*.5*ca*(                           cxxyy); 
                                    abcCoeff2(+1,+1,0) =                               -is*.5*ca*(                           cxxyy);  


                                }

                            }
                            else
                            {
                                const Real cxx = axis==0 ? ca/(12.*dx[0]*dx[0]) : .5*ca/(12.*dx[0]*dx[0]);
                                const Real cyy = axis==1 ? ca/(12.*dx[1]*dx[1]) : .5*ca/(12.*dx[1]*dx[1]);
                                const Real cx  = axis==0 ? -is1*1./(12.*dx[0]*dt) : 0.;                
                                const Real cy  = axis==1 ? -is2*1./(12.*dx[1]*dt) : 0.;

                                abcCoeff( 0,-2,0) = .5*(              -cyy )         +   cy;
                                abcCoeff( 0,-1,0) = .5*(           16.*cyy )         -8.*cy; 
                                abcCoeff(-2, 0,0) = .5*(     -cxx          )  +   cx       ;         
                                abcCoeff(-1, 0,0) = .5*(  16.*cxx          )  -8.*cx       ;         
                                abcCoeff( 0, 0,0) = .5*( -30.*cxx -30.*cyy )               ;
                                abcCoeff(+1, 0,0) = .5*(  16.*cxx          )  +8.*cx       ; 
                                abcCoeff(+2, 0,0) = .5*(     -cxx          )  -   cx       ;
                                abcCoeff( 0,+1,0) = .5*(           16.*cyy )         +8.*cy;
                                abcCoeff( 0,+2,0) = .5*(              -cyy )         -   cy;

                                if( useCompatiblityRBC )
                                {
                  // ghost line 2:
                  // Use:
                  // D_x^2( EM2 BC ) (left/right)
                  // D_y^2( EM2 BC )  (bottom/top) 

                  // -is * D+t D0x D+dD-x W^n + c*Dxx4 (W^{n+1} + W^n)/2

                                    Real cxxxx = axis==0 ?    ca/(dx[0]*dx[0]*dx[0]*dx[0]) : 0.;
                                    Real cxxyy =           .5*ca/(dx[0]*dx[0]*dx[1]*dx[1]);
                                    Real cyyyy = axis==1 ?    ca/(dx[1]*dx[1]*dx[1]*dx[1]) : 0.;
                                    Real cxxx  = axis==0 ? -is1*1./(2.*dx[0]*dx[0]*dx[0]*dt) : 0.;                
                                    Real cyyy  = axis==1 ? -is2*1./(2.*dx[1]*dx[1]*dx[1]*dt) : 0.;

                  // ***check me**
                                    abcCoeff2(-1,-1,0) = .5*(                          cxxyy )          ;
                                    abcCoeff2(+1,-1,0) = .5*(                          cxxyy )          ;
                                    abcCoeff2( 0,-2,0) = .5*(                cyyyy           )  -   cyyy;
                                    abcCoeff2( 0,-1,0) = .5*(            -4.*cyyyy -2.*cxxyy )  +2.*cyyy; 
                                    abcCoeff2(-2, 0,0) = .5*(      cxxxx                     )  -   cxxx;         
                                    abcCoeff2(-1, 0,0) = .5*(  -4.*cxxxx           -2.*cxxyy )  +2.*cxxx;         
                                    abcCoeff2( 0, 0,0) = .5*(  +6.*cxxxx +6.*cyyyy +4.*cxxyy )          ;
                                    abcCoeff2(+1, 0,0) = .5*(  -4.*cxxxx           -2.*cxxyy )  -2.*cxxx; 
                                    abcCoeff2(+2, 0,0) = .5*(      cxxxx                     )  +   cxxx;
                                    abcCoeff2( 0,+1,0) = .5*(            -4.*cyyyy -2.*cxxyy )  -2.*cyyy;
                                    abcCoeff2( 0,+2,0) = .5*(                cyyyy           )  +   cyyy;
                                    abcCoeff2(-1,+1,0) = .5*(                          cxxyy )          ;
                                    abcCoeff2(+1,+1,0) = .5*(                          cxxyy )          ;

                                }
                            }

                        }
                        else
                        {
                            OV_ABORT("FINISH ME - order>4 EM implicit BCs");
                        }
                        if( false )
                        {
              // OLD WAY 
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
                    }
                    else
                    {
                          OV_ABORT(" IMP BC ABORBING -- FINISH ME: curvilinear"); 
                    }

                    FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                    {
                        if( maskLocal(i1,i2,i3) > 0 ) 
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

                            if( useCompatiblityRBC )
                            {
                                int i1m=i1-2*is1, i2m=i2-2*is2, i3m=i3-2*is3; //  2nd ghost point is (i1m,i2m,i3m)
                                setClassify(e,i1m,i2m,i3m, SparseRepForMGF::ghost2);  
                                ForStencil(m1,m2,m3)
                                {
                                    int m  = M123(m1,m2,m3);        // the single-component coeff-index

                                    coeffLocal(m,i1m,i2m,i3m) = abcCoeff2(m1,m2,m3);

                  // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
                                    int j1=i1+m1, j2=i2+m2, j3=i3+m3;                       // the stencil is centred on the boundary pt (i1,i2,i3)
                                    setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                }                 

                            }
                        }

                    } // end FOR_3D

          // fill additional ghost with extrapolation
                    int ghostStart= useCompatiblityRBC ? orderOfAccuracy/2+1 : 2;
                    printF("  ...useCompatiblityRBC=%d, ghostStart=%d\n'",(int)useCompatiblityRBC,ghostStart);
                    for( int ghost=ghostStart; ghost<=numberOfGhostLines; ghost++ )
                    {
              // *wdh* Nov 22, 2023: 
              // getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
              // coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0;
                            bool ok=ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,Ib1,Ib2,Ib3);
                            if( ok )
                            {
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    if( maskLocal(i1,i2,i3) != 0 )
                                    { 
                                        int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                    // --- fill in the coefficients of the extrapolation formula ---
                                        coeffLocal(M0,i1m,i2m,i3m)=0; // zero out all *wdh* Nov 22, 2023: 
                                        for( int m=0; m<=extrapOrder; m++ )
                                        {
                                            coeffLocal(m,i1m,i2m,i3m) = extrapCoeff[m];
                                            int j1=i1m + m*is1, j2=i2m + m*is2, j3=i3m + m*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                            setEquationNumber(m, e,i1m,i2m,i3m,  cc,j1,j2,j3 );      // macro to set equationNumber
                                        }
                                    }
                                } // end FOR_3D
                            }
                    } // end for ghost

                }
                else if(  mg.boundaryCondition(side,axis)> 0 )
                {
                    printF("fill implicit matrix:ERROR: unknown boundaryCondition=%d \n",mg.boundaryCondition(side,axis));
                    OV_ABORT("error");

                }
            } // end for boundary


            if( (addUpwinding && debug>3)  )
            {
        // ::display(impCoeff[grid],sPrintF("implicit time-stepping matrix on grid=%d",grid));
                displayCoeff(impCoeff[grid],sPrintF("AFTER FILL BCS: implicit time-stepping matrix on grid=%d",grid));
                OV_ABORT("stop here for now");
            }




      // --- FILL CORNERS AND EDGES ----

            if( useCompatibility && numberOfDimensions==2 )
            {
        // --- corner points in two dimensions ---
                for( int side2=0; side2<=1; side2++ )
                for( int side1=0; side1<=1; side1++ )
                {

                    const int bc1 = mg.boundaryCondition(side1,0);
                    const int bc2 = mg.boundaryCondition(side2,1);
                    if( bc1>0 && bc2>0 && bc1!=exactBC && bc2!=exactBC )
                    {
                        if( bc1!=dirichlet && bc1!=exactBC && bc1!=neumann )
                        {
                            printF("Un-supported corner bc1 = %d\n",bc1);
                            OV_ABORT("error");
                        }
                        if( bc2!=dirichlet && bc2!=exactBC && bc2!=neumann )
                        {
                            printF("Un-supported corner bc2 = %d\n",bc2);
                            OV_ABORT("error");
                        }

            // symSign = 1 : even symmetry
            //          -1 : odd symmetry
                        Real symSign = +1.; 
                        if( bc1==dirichlet || bc1==exactBC )
                            symSign = -symSign;
                        if( bc2==dirichlet || bc2==exactBC )
                            symSign = -symSign;

                        const int is1 = 1-2*side1;
                        const int is2 = 1-2*side2;
        
            // corner point:
                        const int i1 = mg.gridIndexRange(side1,0);
                        const int i2 = mg.gridIndexRange(side2,1);
                        const int i3 = mg.gridIndexRange(    0,2);
                        if( maskLocal(i1,i2,i3) > 0 )   
                        {
                            for( int m2=1; m2<=numberOfGhostLines; m2++ )
                            for( int m1=1; m1<=numberOfGhostLines; m1++ )
                            {              

                                int j1 = i1-is1*m1, j2=i2-is2*m2, j3=i3; // ghost 
                                int k1 = i1+is1*m1, k2=i2+is2*m2, k3=i3; // interior point 

                // printF("IMPLICIT-MATRIX: Set symmetry condition at corner point (j1,j2,j3)=(%4d,%4d,%4d)\n",j1,j2,j3);

                // Specify that this a "real" equation on the first ghost line: 
                // (A "real" equation has a possible non-zero right-hand-side)
                                setClassify(e,j1,j2,j3, SparseRepForMGF::ghost1);  

                                coeffLocal(M,j1,j2,j3) = 0.0;  // zero out any existing equations

                // The even or odd symmetry equation is 
                //  u(j1,j2,j3,0) - symSign*u(k1,k2,k3,0) = RHS
                                int m=0; 
                                coeffLocal(m,j1,j2,j3) = 1.;
                                setEquationNumber(m, e,j1,j2,j3,  cc,j1,j2,j3 ); 

                                m=1; 
                                coeffLocal(m,j1,j2,j3) = -symSign;
                                setEquationNumber(m, e,j1,j2,j3,  cc,k1,k2,k3 ); 

                            }              

                        }       


                    }

                }
            }
            else if( useCompatibility && numberOfDimensions==3 )  
            {
                int side1,side2,side3;
                int n1a,n1b, n2a,n2b, n3a,n3b;
                for( int edgeDirection=0; edgeDirection<=2; edgeDirection++ ) // direction parallel to the edge
                {
          // There are a total of 4 edges along a given edge direction
                    for( int sidea=0; sidea<=1; sidea++ )
                    for( int sideb=0; sideb<=1; sideb++ )
                    {
                        if( edgeDirection==0 )
                        {
                            side1=0;     side2=sidea; side3=sideb;
                        }
                        else if( edgeDirection==1 )
                        {
                            side1=sideb; side2=0;     side3=sidea;
                        }
                        else
                        {
                            side1=sidea; side2=sideb; side3=0;
                        }
                        int is1=1-2*(side1);
                        int is2=1-2*(side2);
                        int is3=1-2*(side3);
                        int bc1,bc2; 

                        if( edgeDirection==2 )
                        {
                            is3=0;
                            n1a=mg.gridIndexRange(side1,0);
                            n1b=mg.gridIndexRange(side1,0);
                            n2a=mg.gridIndexRange(side2,1);
                            n2b=mg.gridIndexRange(side2,1);
                            n3a=mg.gridIndexRange(    0,2);
                            n3b=mg.gridIndexRange(    1,2);
                            bc1=mg.boundaryCondition(side1,0);
                            bc2=mg.boundaryCondition(side2,1);
                        }
                        else if( edgeDirection==1 )
                        {
                            is2=0;
                            n1a=mg.gridIndexRange(side1,0);
                            n1b=mg.gridIndexRange(side1,0);
                            n2a=mg.gridIndexRange(    0,1);
                            n2b=mg.gridIndexRange(    1,1);
                            n3a=mg.gridIndexRange(side3,2);
                            n3b=mg.gridIndexRange(side3,2);
                            bc1=mg.boundaryCondition(side1,0);
                            bc2=mg.boundaryCondition(side3,2);
                        }
                        else 
                        {
                            is1=0; 
                            n1a=mg.gridIndexRange(    0,0);
                            n1b=mg.gridIndexRange(    1,0);
                            n2a=mg.gridIndexRange(side2,1);
                            n2b=mg.gridIndexRange(side2,1);
                            n3a=mg.gridIndexRange(side3,2);
                            n3b=mg.gridIndexRange(side3,2);
                            bc1=mg.boundaryCondition(side2,1);
                            bc2=mg.boundaryCondition(side3,2);
                        }
                  
                        if( bc1>0 && bc2>0 && bc1!=exactBC && bc2!=exactBC )
                        {
              // -- this is an edge between two physical boundaries --
                            if( bc1!=dirichlet && bc1!=exactBC && bc1!=neumann )
                            {
                                printF("Implicit: Un-supported edge bc1 = %d\n",bc1);
                                OV_ABORT("error");
                            }
                            if( bc2!=dirichlet && bc2!=exactBC && bc2!=neumann )
                            {
                                printF("Implicit: Un-supported edge bc2 = %d\n",bc2);
                                OV_ABORT("error");
                            }

                            Real symSign = +1.; 
                            if( bc1==dirichlet || bc1==exactBC )
                                symSign = -symSign;
                            if( bc2==dirichlet || bc2==exactBC )
                                symSign = -symSign;

              // --- loop over points on the edge ---
                            for( int i3=n3a; i3<=n3b; i3++ )
                            for( int i2=n2a; i2<=n2b; i2++ )
                            for( int i1=n1a; i1<=n1b; i1++ )            
                            { 
                                if( maskLocal(i1,i2,i3)>0 )
                                {
                                    for( int m2=1; m2<=numberOfGhostLines; m2++ )
                                    for( int m1=1; m1<=numberOfGhostLines; m1++ )
                                    { 
                                        int j1,j2,j3, k1,k2,k3;
                                        if( edgeDirection==0 )
                                        {
                      // edge lies along i1=const
                                            j1 = i1; j2=i2-is2*m1; j3=i3-is3*m2; // ghost 
                                            k1 = i1; k2=i2+is2*m1; k3=i3+is3*m2; // interior point  
                                        }            
                                        else if( edgeDirection==1 )
                                        {
                      // edge lies along i2=const
                                            j1 = i1-is1*m1; j2=i2; j3=i3-is3*m2; // ghost 
                                            k1 = i1+is1*m1; k2=i2; k3=i3+is3*m2; // interior point 
                                        } 
                                        else
                                        {
                      // edge lies along i3=constant
                                            j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; // ghost 
                                            k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; // interior point 
                                        } 

                    // printF("IMPLICIT-MATRIX: Set symmetry condition at EDGE point (j1,j2,j3)=(%d,%d,%d)\n",j1,j2,j3);

                    // Specify that this a "real" equation on the first ghost line: 
                    // (A "real" equation has a possible non-zero right-hand-side)
                                        setClassify(e,j1,j2,j3, SparseRepForMGF::ghost1);  

                                        coeffLocal(M,j1,j2,j3) = 0.0;  // zero out any existing equations

                    // The even or odd symmetry equation is 
                    //  u(j1,j2,j3,0) - symSign*u(k1,k2,k3,0) = RHS
                                        int m=0; 
                                        coeffLocal(m,j1,j2,j3) = 1.;
                                        setEquationNumber(m, e,j1,j2,j3,  cc,j1,j2,j3 ); 

                                        m=1; 
                                        coeffLocal(m,j1,j2,j3) = -symSign;
                                        setEquationNumber(m, e,j1,j2,j3,  cc,k1,k2,k3 );                     
                                    }
                                }
                            } // end for 3d
                        } // end if bc1>0 and bc2>0 
                    }
                } // end for edgeDirection

        // --- Vertices in 3D dimensions ---
                for( int side3=0; side3<=1; side3++ )
                for( int side2=0; side2<=1; side2++ )
                for( int side1=0; side1<=1; side1++ )
                {
                    const int bc1 = mg.boundaryCondition(side1,0);
                    const int bc2 = mg.boundaryCondition(side2,1);
                    const int bc3 = mg.boundaryCondition(side3,2);
                    if( bc1>0 && bc2>0 && bc3>0 && bc1!=exactBC && bc2!=exactBC && bc3!=exactBC )
                    {
                        if( bc1!=dirichlet && bc1!=exactBC && bc1!=neumann )
                        {
                            printF("Implicit: Un-supported vertex bc1 = %d\n",bc1);
                            OV_ABORT("error");
                        }
                        if( bc2!=dirichlet && bc2!=exactBC && bc2!=neumann )
                        {
                            printF("Implicit: Un-supported vertex bc2 = %d\n",bc2);
                            OV_ABORT("error");
                        }
                        if( bc3!=dirichlet && bc3!=exactBC && bc3!=neumann )
                        {
                            printF("Implicit: Un-supported vertex bc3 = %d\n",bc3);
                            OV_ABORT("error");
                        }

             // symSign = 1 : even symmetry
            //          -1 : odd symmetry
                        Real symSign = +1.; 
                        if( bc1==dirichlet || bc1==exactBC )
                            symSign = -symSign;
                        if( bc2==dirichlet || bc2==exactBC )
                            symSign = -symSign;
                        if( bc3==dirichlet || bc3==exactBC )
                            symSign = -symSign;            

                        const int is1 = 1-2*side1;
                        const int is2 = 1-2*side2;
                        const int is3 = 1-2*side3;
        
            // Vertex coordinates: 
                        const int i1 = mg.gridIndexRange(side1,0);
                        const int i2 = mg.gridIndexRange(side2,1);
                        const int i3 = mg.gridIndexRange(side3,2);
                        if( maskLocal(i1,i2,i3) > 0 )   
                        {
                            for( int m3=1; m3<=numberOfGhostLines; m3++ )
                            for( int m2=1; m2<=numberOfGhostLines; m2++ )
                            for( int m1=1; m1<=numberOfGhostLines; m1++ )
                            {              

                                int j1=i1-is1*m1, j2=i2-is2*m2, j3=i3-is3*m3; // ghost 
                                int k1=i1+is1*m1, k2=i2+is2*m2, k3=i3+is3*m3; // interior point 

                // printF("IMPLICIT-MATRIX: Set symmetry condition at a 3D VERTEX (j1,j2,j3)=(%4d,%4d,%4d)\n",j1,j2,j3);

                // Specify that this a "real" equation on the first ghost line: 
                // (A "real" equation has a possible non-zero right-hand-side)
                                setClassify(e,j1,j2,j3, SparseRepForMGF::ghost1);  

                                coeffLocal(M,j1,j2,j3) = 0.0;  // zero out any existing equations

                // The even or odd symmetry equation is 
                //  u(j1,j2,j3,0) - symSign*u(k1,k2,k3,0) = RHS
                                int m=0; 
                                coeffLocal(m,j1,j2,j3) = 1.;
                                setEquationNumber(m, e,j1,j2,j3,  cc,j1,j2,j3 ); 

                                m=1; 
                                coeffLocal(m,j1,j2,j3) = -symSign;
                                setEquationNumber(m, e,j1,j2,j3,  cc,k1,k2,k3 ); 

                            }  // end if for m1,m2,m3            

                        }   // end if maskLocal > 0     


                    } // end if bc1>0 and bc2>0 and bc3>0 
                } // end if for side1, sid2, side 3

            } // end if useCompat and nd=3 

        } // end for grid


    // OV_ABORT("stop here for now"); // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if( debug & 16 )
        {
            int grid=0;
            displayCoeff(impCoeff[grid],sPrintF("BEFORE FINISH BC: implicit time-stepping matrix on grid=%d",grid),debugFile);
        } 

        
        if( !useCompatibility )
        {
      // *wdh* add this back Oct 1, 2024

            BoundaryConditionParameters extrapParams;
            extrapParams.orderOfExtrapolation=orderOfAccuracy+1;      // ** TEST

            impCoeff.finishBoundaryConditions(extrapParams); 
            printF("FILL IMP CORNER BCS -- EXTRAP TO ORDER %d\n",extrapParams.orderOfExtrapolation);
        }
        else
        {
            impCoeff.finishBoundaryConditions(); 
        }

        if( debug & 8 )
        {
            int grid=0;
            displayCoeff(impCoeff[grid],sPrintF("AFTER FINISH BC: implicit time-stepping matrix on grid=%d",grid),debugFile);
            fflush(debugFile);
            OV_ABORT("stop here for now");
        }    

        impSolver.setCoefficientArray( impCoeff );   // supply coefficients to Oges

    }
    timing(timeForInitialize) += getCPU()-cpu0;

  
    return 0;
}
