// This file automatically generated from solveHelmholtzDirect.bC with bpp.
// ----Form and solve the Helmholtz equation using a direct or iterative solver ----

#include "CgWaveHoltz.h"
#include "CgWave.h"
// #include "Overture.h" 
// #include "MappedGridOperators.h"
#include "Oges.h"
#include "ParallelUtility.h"
#include "CompositeGridOperators.h"
#include "SparseRep.h"

#include "PlotStuff.h"
#include "GL_GraphicsInterface.h"


#include "gridFunctionNorms.h"


#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )


// Use this for indexing into coefficient matrices representing systems of equations
#define CE(c,e) (stencilSize*((c)+numberOfComponentsForCoefficients*(e)))
#define M123(m1,m2,m3) (m1+halfWidth1+width*(m2+halfWidth2+width*(m3+halfWidth3)))
#define M123CE(m1,m2,m3,c,e) (M123(m1,m2,m3)+CE(c,e))


#define ForStencil(m1,m2,m3)   for( m3=-halfWidth3; m3<=halfWidth3; m3++) for( m2=-halfWidth2; m2<=halfWidth2; m2++) for( m1=-halfWidth1; m1<=halfWidth1; m1++) 

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

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

// =======================================================================
//  Macro to zero out the matrix coefficients for equations e1,e1+1,..,e2
// =======================================================================
#define zeroMatrixCoefficients( coeff,e1,e2, i1,i2,i3 )for( int m=CE(0,e1); m<=CE(0,e2+1)-1; m++ ) coeff(m,i1,i2,i3)=0.


// ==========================================================================
// Macro: fill the matrix with extrapolation for a given ghost=1,2,3,...
// ==========================================================================

// ==========================================================================
// Macro: fill the matrix with extrapolation formula
//
// Input: eq,uc : equation and component number 
// Input: i1g,i2g,i3g : point to extrapolate (unused point)
// Input: is1,is2,is3 : direction to extrpolate 
// 
// ==========================================================================


// ==========================================================================
// Macro: setup the variables needed to fill a sparse matrix on a mappedGrid
// ==========================================================================

// From advWave.bf90

// ! ===========================================================================================
// ! Macro: compute the coefficients in the sosup dissipation for curvilinear grids
// ! ===========================================================================================
// #beginMacro getSosupDissipationCoeff2d(upwindCoefficient)
//  do dir=0,1
//    ! diss-coeff ~= 1/(change in x along direction r(dir) )
//    ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
//    upwindCoefficient(dir) = betaUpwind*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
//  end do
// #endMacro

// #beginMacro getSosupDissipationCoeff3d(upwindCoefficient)
//  do dir=0,2
//    ! diss-coeff ~= 1/(change in x along direction r(dir) )
//    ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
//    upwindCoefficient(dir) = betaUpwind*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2  + rsxy(i1,i2,i3,dir,2)**2 )/dr(dir) 
//  end do
//  ! write(*,'(" upwindCoefficient =",3(1pe10.2))') (upwindCoefficient(dir),dir=0,2)
// #endMacro



// ===============================================================
// Macro: setup OgesParmeters                                  
// ===============================================================


// ===============================================================
// Macro: fill the coefficients for upwind dissipation
//
// Upwind dissipation takes the form
//      w.tt = L w + beta * M w_t 
//          M =  -(h^2 D+D-)^2 
//  Setting
//      w = u exp( I*sigma*omega*t ) 
//      u = ur + I*ui 
// Gives
//     -omega^2 u = L u + I beta*sigma*omega M u 
// or 
//   -omega^2 ur = L ur  - beta*sigma*omega M ui
//   -omega^2 ui = L ui  + beta*sigma*omega M ur
//  
// ===============================================================

// =============================================================
// Macro: initialize the upwind coefficients
// =============================================================

// ============================================================================================
// ============================================================================================

// ============================================================================================
/// With superGrid there are two sets of metrics  
// ============================================================================================

// ============================================================================================
/// \brief Form and solve the Helmholtz equation using a direct or iterative solver.
// ============================================================================================
int CgWaveHoltz::solveHelmholtzDirect( realCompositeGridFunction & u, realCompositeGridFunction & f  )
{
    int debug=1; 

    const real & omega     = dbase.get<real>("omega");

    CgWave & cgWave        = *dbase.get<CgWave*>("cgWave");
    const real & c         = cgWave.dbase.get<real>("c");

    const int & orderOfAccuracy          = cgWave.dbase.get<int>("orderOfAccuracy");
    const int & filterTimeDerivative     = cgWave.dbase.get<int>("filterTimeDerivative");
    const int & useSuperGrid             = cgWave.dbase.get<int>("useSuperGrid");
    const IntegerArray & superGrid       = cgWave.dbase.get<IntegerArray>("superGrid");
    const int & adjustPlotsForSuperGrid  = cgWave.dbase.get<int>("adjustPlotsForSuperGrid");
    const int & adjustErrorsForSuperGrid = cgWave.dbase.get<int>("adjustErrorsForSuperGrid");
    const int & solveForScatteredField   = cgWave.dbase.get<int>("solveForScatteredField");
    
    const Real & damp                    = cgWave.dbase.get<Real>("damp");
    const Real & dampSave                = cgWave.dbase.get<Real>("dampSave");
    const int & upwind                   = cgWave.dbase.get<int>("upwind");  

    const int & numberOfFrequencies      = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray     = dbase.get<RealArray>("frequencyArray");

    const bool solveComplex = filterTimeDerivative; 


    printF("\n ============== DIRECT SOLVE OF THE HELMHOLTZ EQUATION =============\n");
    printF("    c=%.4g, numberOfFrequencies=%d, orderOfAccuracy=%d, solveComplex=%d\n", c,numberOfFrequencies,orderOfAccuracy,(int)solveComplex);
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
        printF(" freq(%d) = %.6g\n",freq,frequencyArray(freq));
    }
    if( damp!=0. )
        printF(" Damping is on: damp=%14.6e, dampSave=%14.6e\n",damp,dampSave);
    printF(" useSuperGrid=%d, superGrid[grid] = [",useSuperGrid);
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        printF("%i,",superGrid(grid));
    printF("]\n");
  // ::display(superGrid,"superGrid","%2i");
    printF(" ===================================================================\n");

  // ::display(frequencyArray,"frequencyArray");

    if( omega != frequencyArray(0) )
    {
        printF("solveHelmholtzDirect: ERROR omega=%g, != frequencyArray(0)=%g (TEMP SANITY CHECK -- FIX ME)\n",omega,frequencyArray(0));
        OV_ABORT("error"); 
    }


  // --- we save the number of iterations ---
    if( !dbase.has_key("helmholtzSolverIterations") )
        dbase.put<int>("helmholtzSolverIterations")=0;
    int & helmholtzSolverIterations =  dbase.get<int>("helmholtzSolverIterations");
  // --- we save the name of the sparse solver we use ---
    if( !dbase.has_key("helmholtzSolverName") )
        dbase.put<aString>("helmholtzSolverName");
    aString & helmholtzSolverName =  dbase.get<aString>("helmholtzSolverName");
    

    CompositeGrid & cg = *u.getCompositeGrid();
    const int & numberOfComponentGrids = cg.numberOfComponentGrids(); 
    const int & numberOfDimensions     = cg.numberOfDimensions();



  // --------------------------------------------------------
  // ------ Fill in the forcing and boundary conditions -----
  // --------------------------------------------------------
    cgWave.getHelmholtzForcing( f );



  // --- DO THIS FOR NOW : **FIX ME** -----

    const int planeWave=1, gaussianPlaneWave=2, boxHelmholtz=3;

    Real kxBoxHelmholtz=1., kyBoxHelmholtz=1., kzBoxHelmholtz=1.;
    int knownSolutionOption=0; // no known solution
    if( cgWave.dbase.has_key("userDefinedKnownSolutionData") )
    {
    // printF("++++ solveHelmholtzDirect: userDefinedKnownSolutionData is found\n");
        DataBase & db =  cgWave.dbase.get<DataBase>("userDefinedKnownSolutionData");
        const aString & userKnownSolution = db.get<aString>("userKnownSolution");
        if( userKnownSolution=="planeWave"  )
        {
            knownSolutionOption=planeWave;                   // this number must match in bcOptWave.bf90
        }
        else if( userKnownSolution=="gaussianPlaneWave"  ) 
        {
            knownSolutionOption=gaussianPlaneWave;           // this number must match in bcOptWave.bf90
        }    
        else if( userKnownSolution=="boxHelmholtz"  ) 
        {
            knownSolutionOption=boxHelmholtz;                // this number must match in bcOptWave.bf90
            kxBoxHelmholtz = cgWave.dbase.get<Real>("kxBoxHelmholtz");
            kyBoxHelmholtz = cgWave.dbase.get<Real>("kyBoxHelmholtz");
            kzBoxHelmholtz = cgWave.dbase.get<Real>("kzBoxHelmholtz"); 
            printF("\n @@@@ solveHelmholtzDirect: boxHelmholtz: kxBoxHelmholtz=%g\n\n",kxBoxHelmholtz);     
        }
        else if( userKnownSolution=="polyPeriodic"  ) 
        {
            knownSolutionOption=4;                   // this number must match in bcOptWave.bf90
        } 

    } 
    else
    {
     // printF("+++ solveHelmholtzDirect: userDefinedKnownSolutionData is NOT found\n");
    } 

    IntegerArray boundaryConditions(2,3,numberOfComponentGrids);
    boundaryConditions = OgesParameters::dirichlet; // default
    RealArray bcData(2,2,3,numberOfComponentGrids);
    Range all; 

  // Solve constCoeff(0,grid)*I +constCoeff(1,grid)*Laplacian 
  // We solve:  omega^2 I + c^2 Delta = f 
    RealArray constantCoeff(2,numberOfComponentGrids);
  // constantCoeff(0,all) = -SQR(omega);    // FIX SIGN IN CGWAVE 
  // constantCoeff(1,all) = -SQR(c); 

  // Assign boundary conditions for Oges
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        ForBoundary(side,axis)
        {
              if( mg.boundaryCondition(side,axis)==CgWave::dirichlet ||
                      mg.boundaryCondition(side,axis)==CgWave::absorbing )  // *********** DO THIS FOR NOW use dirichlet BCs for supergrid
              {
                  boundaryConditions(side,axis,grid) = OgesParameters::dirichlet;
              }
              else if( mg.boundaryCondition(side,axis)==CgWave::neumann )
              { 
                  boundaryConditions(side,axis,grid) = OgesParameters::neumann;
              }
              else if( mg.boundaryCondition(side,axis) > 0 )
              {
                  printF("CgWave::solveHelmholtzDirect:ERROR: unknown boundaryCondition=%d for (side,axis,grid)=(%d,%d,%d)\n", mg.boundaryCondition(side,axis),side,axis,grid);
                  OV_ABORT("ERROR");
              }

        }

    }  

    CompositeGridOperators op(cg);
    op.setOrderOfAccuracy(orderOfAccuracy);

  // -- set up extrapolation formula --
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
      printF("solveHelmholtzDirect:: unexpected extrapOrder=%d\n",extrapOrder);
      OV_ABORT("ERROR");
    }

    if( !solveComplex )
    {
    //  --- solve the Real Helmholtz problem ----

        for( int freq=0; freq<numberOfFrequencies; freq++ )
        {
            const Real omega = frequencyArray(freq);

            Oges solver( cg );                     // create a solver

                if( dbase.has_key("helmholtzOgesParameters") )
                {  
                    printf("solveHelmholtzDirect: OgesParameters found. Changing the direct Helmholtz solver parameters. \n");
                    OgesParameters & par = dbase.get<OgesParameters>("helmholtzOgesParameters");
                    solver.setOgesParameters(par);
                }
                else
                {
                    int solverType=OgesParameters::yale; 
          // solverType=OgesParameters::PETSc;
          // solverType=OgesParameters::PETScNew; // parallel
                    solver.set(OgesParameters::THEsolverType,solverType); 
                    if( solverType==OgesParameters::PETSc )
                      solver.set(OgesParameters::THEsolverMethod,OgesParameters::biConjugateGradientStabilized);
          // solver.set(OgesParameters::THEparallelSolverMethod,OgesParameters::gmres);
          // solver.set(OgesParameters::THErelativeTolerance,max(tol,REAL_EPSILON*10.));
          // solver.set(OgesParameters::THEmaximumNumberOfIterations,10000);
          // if( iluLevels>=0 )
          //   solver.set(OgesParameters::THEnumberOfIncompleteLULevels,iluLevels);
                }
                if( upwind )
                {
                    Real fillinRatio=0;
          //solver.get(OgesParameters::THEfillinRatio,fillinRatio );
          // printF("Current fillinRatio = %g\n",fillinRatio);
                    int stencilWidth = orderOfAccuracy+1;
                    int stencilSize=int( pow(stencilWidth,cg.numberOfDimensions())+1 );  // add 1 for interpolation equations
                    fillinRatio = (stencilSize+2)*2 + 20.; // what should this be ?
                    printF("Set new fillinRatio = %g\n",fillinRatio);
                    solver.set(OgesParameters::THEfillinRatio,fillinRatio );
                }

            if( freq==0 )
            {
                helmholtzSolverName = solver.parameters.getSolverName();
                printF("\n === Direct Helmholtz Solver:\n %s\n =====\n",(const char*)helmholtzSolverName);
            }

  

      // -- use predefined equations : 
            constantCoeff(0,all) = SQR(omega);    
            constantCoeff(1,all) = SQR(c); 
            bcData=0.;
            solver.setEquationAndBoundaryConditions( OgesParameters::heatEquationOperator ,op,boundaryConditions,bcData,constantCoeff );

      // realCompositeGridFunction u(cg,all,all,all);
      // realCompositeGridFunction f(cg,all,all,all);

      // --- fill in boundary conditions --- 
      //  Boundary conditions will be filled in by getHelmholtzForcing, but these are not correct for numberOfFrequencies>1 
            Index Ib1,Ib2,Ib3; 
            Index Ig1,Ig2,Ig3; 
            for( int grid=0; grid<numberOfComponentGrids; grid++ )
            {
                MappedGrid & mg = cg[grid];
                for( int axis=0; axis<cg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        if( mg.boundaryCondition(side,axis)==CgWave::dirichlet )
                        {
                              getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                              OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);
                              OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);

                              if( knownSolutionOption==boxHelmholtz )
                              {
                  // const Real kx = kxBoxHelmholtz + twoPi*freq;
                  // const Real ky = kyBoxHelmholtz + twoPi*freq;
                  // const Real kz = kzBoxHelmholtz + twoPi*freq;
                                    const Real rpar[] = {0.,kxBoxHelmholtz/twoPi,kyBoxHelmholtz/twoPi,kzBoxHelmholtz/twoPi}; // 
                                    Real omega,kx,ky,kz; 
                      // This macro is used in: 
                      //    userDefinedKnownSolution.bC
                      //    userDefinedForcing.bC
                      //    solveHelmholtz.bC 
                                            omega = frequencyArray(freq);
                                            kx  = rpar[1]*twoPi*(freq*.5+1.);
                                            ky  = rpar[2]*twoPi*(freq*.5+1.);
                                            kz  = rpar[3]*twoPi*(freq*.5+1.);  
                      // kx  = (rpar[1]+freq)*twoPi;
                      // ky  = (rpar[2]+freq)*twoPi;
                      // kz  = (rpar[3]+freq)*twoPi; 
                                                                  
                                    if( mg.numberOfDimensions()==2 )
                                        fLocal(Ib1,Ib2,Ib3,freq)= sin( kx*xLocal(Ib1,Ib2,Ib3,0) )*sin( ky*xLocal(Ib1,Ib2,Ib3,1) );
                                    else
                                        fLocal(Ib1,Ib2,Ib3,freq)= sin( kx*xLocal(Ib1,Ib2,Ib3,0) )*sin( ky*xLocal(Ib1,Ib2,Ib3,1) )*sin( kz*xLocal(Ib1,Ib2,Ib3,2) );
                              }

                        }
                        else if( mg.boundaryCondition(side,axis)==CgWave::neumann )
                        {
                              getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
                              OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);
               // fLocal(Ig1,Ig2,Ig3)=0.;              // Neumann BC -- fill ghost values
                        }
                        else if( mg.boundaryCondition(side,axis)>0 )
                        {
                            printF("CgWave::solveHelmholtzDirect:ERROR: unknown boundaryCondition=%d for (side,axis,grid)=(%d,%d,%d)\n", mg.boundaryCondition(side,axis),side,axis,grid);
                            OV_ABORT("ERROR");
                        }   
                    }
                }
            }

            real time0=getCPU();

      // ------- SOLVE THE HELMHOLTZ EQUATIONS -----
            if( numberOfFrequencies==1 )
            {
                u=0.;  // initial guess for iterative solvers
                solver.solve( u,f ); 
            }
            else
            {
                realCompositeGridFunction uTemp(cg,all,all,all), fTemp(cg,all,all,all);
                uTemp=0.;  // initial guess for iterative solvers
                Index I1,I2,I3;
        // copy f(.,.,.,freq) to fTemp
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);
                    OV_GET_SERIAL_ARRAY(Real,fTemp[grid],fTempLocal);
                    getIndex(cg[grid].dimension(),I1,I2,I3);
                    fTempLocal(I1,I2,I3)=fLocal(I1,I2,I3,freq);
                }

                solver.solve( uTemp,fTemp ); 

                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
                    OV_GET_SERIAL_ARRAY(Real,uTemp[grid],uTempLocal);
                    getIndex(cg[grid].dimension(),I1,I2,I3);
                    uLocal(I1,I2,I3,freq)=uTempLocal(I1,I2,I3);
                }
            }  


            real time= ParallelUtility::getMaxValue(getCPU()-time0);
            helmholtzSolverIterations += solver.getNumberOfIterations();
            printF("\n*** freq=%d: omega=%.5g, max residual=%8.2e, time for direct Helmholtz solve = %8.2e (s) (iterations=%i) ***\n",
                          freq,omega,solver.getMaximumResidual(),time,helmholtzSolverIterations);
      

        } // end for freq 
    }
    else
    {
    // ----------------------------------------------------------
    // ----------------------------------------------------------
    // ------------- Solve the complex Helmholtz ----------------
    // ----------------------------------------------------------
    // ----------------------------------------------------------

        const int numberOfComponents = 2;

        IntegerArray & superGrid = cgWave.dbase.get<IntegerArray>("superGrid");

        Range all;
        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
        Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
        Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];
        int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];    
        int i1,i2,i3, j1,j2,j3, m1,m2,m3;

    // make a grid function to hold the coefficients
        assert( orderOfAccuracy==2 );

        int stencilWidth = orderOfAccuracy+1;
    // if( upwind )
    //   stencilWidth +=2;
    // printF(" **** SETTING stencilWidth=%d *****\n",stencilWidth);


        int stencilSize=int( pow(stencilWidth,cg.numberOfDimensions())+1 );  // add 1 for interpolation equations

      
        bool isAllRectangular=true;
        bool isAllCurvilinear=true;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            bool isRectangular = mg.isRectangular();
            isAllRectangular = isAllRectangular && isRectangular;
            isAllCurvilinear  = isAllCurvilinear &&  !isRectangular;
        }   
        if( upwind  )
        {
      // Upwind: there are 4 extra coefficients in 2D for EACH equation
      //                    X 4
      //                    |                Cartesian: 
      //                    O                 3  O  4
      //                    |                    |
      //            X---O---O---O---X         O--O--O
      //            1       |       2            |
      //                    O                 1  O  2
      //                    |
      //                    X 3   
      // Cartesian case: Store extra upwind in unused points
      // Curvilinear case: we need extra storage in the coeff array


            if( !isAllRectangular )
            {
                int extraStencil = 3;   // we already have 1 extra for interpolation
                stencilSize += extraStencil; 
                printF(" **** SETTING stencilSize=%d, extraStencil=%d (for wider upwind formula) *****\n",stencilSize,extraStencil);
            }
        }

        int stencilDimension=stencilSize*SQR(numberOfComponents);
        realCompositeGridFunction coeff(cg,stencilDimension,all,all,all); 
    // make this grid function a coefficient matrix:
    // int numberOfGhostLines= 1;
        int numberOfGhostLines=(stencilWidth-1)/2;
        if( upwind )
        {
            numberOfGhostLines += 1;   // for extrapolate an extra ghost line when upwinding
        }    

        coeff.setIsACoefficientMatrix(TRUE,stencilSize,numberOfGhostLines,numberOfComponents);
        coeff=0.;

        op.setStencilSize(stencilSize);
            
    // create grid functions to hold two components (real and imag parts): 
    // realCompositeGridFunction uc(cg,all,all,all,numberOfComponents),
        realCompositeGridFunction fc(cg,all,all,all,numberOfComponents);

        op.setNumberOfComponentsForCoefficients(numberOfComponents);
    // uc.setOperators(op);                              // associate differential operators with u
        coeff.setOperators(op);

        const int eq1=0, eq2=1;   // equation numbers
        const int urc=0, uic=1;      // component numbers
        const int ic1=0, ic2=1;      // component numbers

    // Here are the interior equations: 
    //    c^2 * Delta u + omega^2 u + sigma*omega*damp*v = f(x)
    //    c^2 * Delta v + omega^2 v - sigma*omega*damp*u = 0 

    // Real damp  =  1.;    // ** FIX ME ***  
        Real sigma     = -1.;       // SIGN FOR exp( sigma*I*omega*t)    ** FIX ME ***  
        Real omegaSign = -1.; // ** FIX ME ****************************************

        const Real a1 = omega*omega;
        const Real a2 = omega*omega;
        const Real b1 =  sigma*omega*damp;
        const Real b2 = -sigma*omega*damp; 
        const Real cSq = c*c;


        if( !useSuperGrid || isAllCurvilinear )
        {
      // we can use the standard operators
            coeff=( cSq*op.laplacianCoefficients(eq1,urc) + a1*op.identityCoefficients(eq1,urc) + b1*op.identityCoefficients(eq1,uic)+
                            cSq*op.laplacianCoefficients(eq2,uic) + a2*op.identityCoefficients(eq2,uic) + b2*op.identityCoefficients(eq2,urc)   );
        }
        else
        {
      // -- Some Cartesian Grids + SuperGrid
            const int numberOfComponentsForCoefficients =2;
          
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const bool isRectangular = mg.isRectangular(); 

                realMappedGridFunction & mgCoeff = coeff[grid];
                MappedGridOperators & mgop       = op[grid];

                if( !isRectangular || !superGrid(grid) )
                {
          // Curvilinear grid, or Cartesian and no superGrid
                    if( 1==0 )
                    { // trouble here : why ?
                        mgCoeff =( cSq*mgop.laplacianCoefficients(eq1,urc) + a1*mgop.identityCoefficients(eq1,urc) + b1*mgop.identityCoefficients(eq1,uic)+
                                              cSq*mgop.laplacianCoefficients(eq2,uic) + a2*mgop.identityCoefficients(eq2,uic) + b2*mgop.identityCoefficients(eq2,urc)   );
                    }
                    else
                    {
                      
                        Range M0 = stencilSize;  // for first equation
                        Range M1 = M0 + CE(1,1); // for second equation

                        OV_GET_SERIAL_ARRAY(Real,mgCoeff,coeffLocal);

                        getIndex(mg.dimension(),I1,I2,I3);
                        RealArray dCoeff(M0,I1,I2,I3);

            // --- Coefficients for Laplacian --
                        mgop.assignCoefficients( MappedGridOperators::laplacianOperator,dCoeff,I1,I2,I3,0,0);

                        coeffLocal(M0,I1,I2,I3) = cSq*( dCoeff(M0,I1,I2,I3) );  
                        coeffLocal(M1,I1,I2,I3) = cSq*( dCoeff(M0,I1,I2,I3) ); 

            // ---- terms involving the identity ------
                        mgop.assignCoefficients( MappedGridOperators::identityOperator,dCoeff,I1,I2,I3,0,0);

                        coeffLocal(M0,I1,I2,I3) += a1*dCoeff(M0,I1,I2,I3);   // omega^2 I 
                        coeffLocal(M1,I1,I2,I3) += a2*dCoeff(M0,I1,I2,I3);   // omega^2 I 

                        if( damp!=0 )
                        { // ---- damping terms ---
                            coeffLocal(M0+CE(1,0),I1,I2,I3) += b1*dCoeff(M0,I1,I2,I3);   // b1*ui in ur eqn
                            coeffLocal(M0+CE(0,1),I1,I2,I3) += b2*dCoeff(M0,I1,I2,I3);   // b2*ur in ui eqn
                        }


                    }

                }
                else
                {
          // Cartesian grid + super grid 
                    printF("&&&&&& solveHelmholtzDirect: Adjust Cartesian grid=%d for superGrid &&&&&&&\n",grid);

          // useAbsorbingLayer(axis,grid) = 1 if this axis has a superGridLayer 
                    IntegerArray & useAbsorbingLayer = cgWave.dbase.get<IntegerArray>("useAbsorbingLayer");

                    Range M0 = stencilSize;  // for first equation
                    Range M1 = M0 + CE(1,1); // for second equation

          // Range M = mgCoeff.dimension(0);

                    OV_GET_SERIAL_ARRAY(Real,mgCoeff,coeffLocal);

                    getIndex(mg.dimension(),I1,I2,I3);
                    RealArray ddCoeff(M0,I1,I2,I3), dCoeff(M0,I1,I2,I3);

          // --- Coefficients for xx and x derivatives ---
                    mgop.assignCoefficients( MappedGridOperators::xxDerivative,ddCoeff,I1,I2,I3,0,0);
                    mgop.assignCoefficients( MappedGridOperators::xDerivative,  dCoeff,I1,I2,I3,0,0);

                    if( useAbsorbingLayer(0,grid) )
                    {
            // -- scale coefficients using superGrid functions --
                        RealArray *& etaxSuperGrid = cgWave.dbase.get<RealArray*>("etaxSuperGrid" );
                        RealArray & etax = etaxSuperGrid[grid];
                        for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )
                        {
                            for( int m=0; m<stencilSize; m++ )
                            {
                                ddCoeff(m,i1,I2,I3) *= etax(i1,0);  // scale by "(r.x)^2"
                                dCoeff(m,i1,I2,I3)  *= etax(i1,1);  // scale by "r.xx"
                            }
                        }
                    }

                    coeffLocal(M0,I1,I2,I3) = cSq*( ddCoeff(M0,I1,I2,I3)+ dCoeff(M0,I1,I2,I3) );    // transformed .xx derivative
                    coeffLocal(M1,I1,I2,I3) = cSq*( ddCoeff(M0,I1,I2,I3)+ dCoeff(M0,I1,I2,I3) ); 

          // --- Coefficients for yy and y derivatives ---
                    mgop.assignCoefficients( MappedGridOperators::yyDerivative,ddCoeff,I1,I2,I3,0,0);
                    mgop.assignCoefficients( MappedGridOperators::yDerivative,  dCoeff,I1,I2,I3,0,0);

                    if( useAbsorbingLayer(1,grid) )
                    {
            // -- scale coefficients using superGrid functions --  
                        RealArray *& etaySuperGrid = cgWave.dbase.get<RealArray*>("etaySuperGrid" );
                        RealArray & etay = etaySuperGrid[grid];
                        for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )
                        {
                            for( int m=0; m<stencilSize; m++ )
                            {
                                ddCoeff(m,I1,i2,I3) *= etay(i2,0);
                                dCoeff(m,I1,i2,I3)  *= etay(i2,1);
                            }
                        }
                    }

                    coeffLocal(M0,I1,I2,I3) += cSq*( ddCoeff(M0,I1,I2,I3)+ dCoeff(M0,I1,I2,I3) );    // transformed .yy derivative
                    coeffLocal(M1,I1,I2,I3) += cSq*( ddCoeff(M0,I1,I2,I3)+ dCoeff(M0,I1,I2,I3) ); 

          // ---- terms involving the identity ------
                    mgop.assignCoefficients( MappedGridOperators::identityOperator,dCoeff,I1,I2,I3,0,0);

                    coeffLocal(M0,I1,I2,I3) += a1*dCoeff(M0,I1,I2,I3);   // omega^2 I 
                    coeffLocal(M1,I1,I2,I3) += a2*dCoeff(M0,I1,I2,I3);   // omega^2 I 

                    if( damp!=0 )
                    { // ---- damping terms ---
                        coeffLocal(M0+CE(1,0),I1,I2,I3) += b1*dCoeff(M0,I1,I2,I3);   // b1*ui in ur eqn
                        coeffLocal(M0+CE(0,1),I1,I2,I3) += b2*dCoeff(M0,I1,I2,I3);   // b2*ur in ui eqn
                    }

                }


            }  
      // OV_ABORT("solveHelmholtzDirect: Finish Cartesian grid with Supergrid"); 

        }




        Real upwindCoefficient[3]; // holds upwind-diss coeff
        int upwindWidth=5;
    // upwindWidth=3;    // ************USE LOWER ORDER UPWIND FOR NOW ********
        int upwindHalfWidth=(upwindWidth-1)/2;

        Range Ruw(-upwindHalfWidth,upwindHalfWidth);
        RealArray upwindWeights(Ruw);    

        if( upwindWidth==5 )
        {
      //  here are weights in -( -h^2 D+D-)^2 
            upwindWeights(-2) = -1.;
            upwindWeights(-1) =  4.;
            upwindWeights( 0) = -6.;
            upwindWeights( 1) =  4.;
            upwindWeights( 2) = -1.; 
        }
        else if( upwindWidth==3 )
        {
      //  here are weights in -( -h^2 D+D-)
            upwindWeights(-1) =  1.;
            upwindWeights( 0) = -2.;
            upwindWeights( 1) =  1.;
        }   

        if( upwind )
        {
      // --------- add upwind dissipation --------

          

            printF("\n ########## ADD UPWIND TO DIRECT HELMHOLTZ SOLVER Upwind-stencil=%d ####### \n\n", upwindWidth);

      // Upwind dissipation takes the form
      //      w.tt = L w + beta * M w_t 
      //          M =  -(h^2 D+D-)^2 
      //  Setting
      //      w = u exp( I*sigma*omega*t ) 
      //      u = ur + I*ui 
      // Gives
      //     -omega^2 u = L u + I beta*sigma*omega M u 
      // or 
      //   -omega^2 ur = L ur  - beta*sigma*omega M ui
      //   -omega^2 ui = L ui  + beta*sigma*omega M ur      

      // OV_ABORT(" *** FINISH ME ***");
  
      // int i1m,i2m,i3m;
      // int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];

            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const bool isRectangular = mg.isRectangular();
                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

                realMappedGridFunction & coeffg = coeff[grid];
                MappedGridOperators & mgop = *coeffg.getOperators();
                
                // SparseRepForMGF & sparse = *coeffg.sparse;
                // const int numberOfComponentsForCoefficients = sparse.numberOfComponents;  // size of the system of equations
                // const int numberOfGhostLines                = sparse.numberOfGhostLines;
                // const int stencilSize                       = sparse.stencilSize;
                // const int stencilDim=stencilSize*numberOfComponentsForCoefficients; // number of coefficients per equation
                // const int equationOffset      = sparse.equationOffset;
                // IntegerArray & equationNumber = sparse.equationNumber;
                // IntegerArray & classify       = sparse.classify;
                // const int equationNumberBase1  =equationNumber.getBase(1);
                // const int equationNumberLength1=equationNumber.getLength(1);
                // const int equationNumberBase2  =equationNumber.getBase(2);
                // const int equationNumberLength2=equationNumber.getLength(2);
                // const int equationNumberBase3  =equationNumber.getBase(3);
                // const int orderOfAccuracy=mgop.getOrderOfAccuracy(); assert( orderOfAccuracy==2 );
                // // stencil width's and half-width's : 
                // // const int width = orderOfAccuracy+1;
                // const int width = stencilWidth;
                // const int halfWidth1 = (width-1)/2;
                // const int halfWidth2 = numberOfDimensions>1 ? halfWidth1 : 0;
                // const int halfWidth3 = numberOfDimensions>2 ? halfWidth1 : 0;
                    assert( coeffg.sparse!=NULL );
                    SparseRepForMGF & sparse = *coeffg.sparse;
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
                    Range M0 = stencilSize;  
                    Range M = coeffg.dimension(0);
                    Range M1 = M0 + CE(1,1); // for second equation

                OV_GET_SERIAL_ARRAY(Real,coeffg,coeffLocal);

                    Real dx[3]={1.,1.,1.};
                    Real dr[3]={1.,1.,1.};
                    if( isRectangular )
                    { // rectangular grid grid-spacings: 
                        mg.getDeltaX(dx);
                    }
                    else
                    {
            // unit square grid spacings: 
                        for( int dir=0; dir<3; dir++ )
                            dr[dir]=mg.gridSpacing(dir);           
                    }        

        // get some constants in the upwinding:
          // Here is the upwind prefactor : 
          //    In 2D: betaUpwind = (sigma*omega*c)/( sqrt(2)* 8 )
                    const Real dtUpwind=1.;             // dt to use when computing upwind coeff
                    const bool adjustForTimeStep=false; // do not adjust for dt
                    Real upwindDissipationCoefficient = cgWave.getUpwindDissipationCoefficient( grid,dtUpwind,adjustForTimeStep );
                    if( true )
                    {
                        Real adSosup = (c)/( sqrt(1.*numberOfDimensions) * pow(2.,orderOfAccuracy+1) );
                        printF(">>>>>>solveHelmholtzDirect: upwindDissipationCoefficient=%g, old-way=%g, omega=%g\n",upwindDissipationCoefficient,adSosup,omega);
                    }
          // Real betaUpwind = (omegaSign*omega*c)/( sqrt(1.*numberOfDimensions) * pow(2.,orderOfAccuracy+1) ); // OLD
                    Real betaUpwind = (omegaSign*omega)*upwindDissipationCoefficient;
          // betaUpwind *= 10; // ** TEST ***
          // betaUpwind =0.; // *** TEMP TEST
                    if( isRectangular )
                    {
            // --- Here is the upwind coefficient for Cartesian grids ---
            // const Real beta= (1./(8.*sqrt(2.)))*c;
                        upwindCoefficient[0] = betaUpwind/dx[0];  
                        upwindCoefficient[1] = betaUpwind/dx[1];
                    } 

                if( !isRectangular )
                    mg.update(MappedGrid::THEinverseVertexDerivative );

                    bool useOriginalMetrics = useSuperGrid && superGrid(grid) && !isRectangular;
          // OV_GET_SERIAL_ARRAY_CONDITIONAL(Real,mg.inverseVertexDerivative(),rsxyLocal2,!useOriginalMetrics);
                    OV_GET_SERIAL_ARRAY_CONDITIONAL(Real,mg.inverseVertexDerivative(),rsxyLocalNew,!useOriginalMetrics);
          // if( useSuperGrid && superGrid(grid) )
          // {
          //   // -- SuperGrid changes the metrics but upwinding uses the original metrics: 
          //   RealArray *& rxOriginal = cgWave.dbase.get<RealArray*>("rxOriginal");
          //   // rsxyLocal.reference(rxOriginal[grid]);
          // }
                    RealArray & rsxyLocal = useOriginalMetrics ? cgWave.dbase.get<RealArray*>("rxOriginal")[grid] : rsxyLocalNew;

        // bool useOriginalMetrics = superGrid(grid);
        // OV_GET_SERIAL_ARRAY_CONDITIONAL(Real,mg.inverseVertexDerivative(),rsxyLocal,!useOriginalMetrics);
        // if( superGrid(grid) )
        // {
        //   // -- SuperGrid changes the metrics but upwinding uses the original metrics: 
        //   RealArray *& rxOriginal = cgWave.dbase.get<RealArray*>("rxOriginal");
        //   rsxyLocal.reference(rxOriginal[grid]);
        // }

        // macro to turn 4D array rxsyLocal(:,:,:,*) into 5D array RXLocal(:,:,:,m,n)
                #define RXLocal(i1,i2,i3,m,n) rsxyLocal(i1,i2,i3,(m)+numberOfDimensions*(n))

                
                if( numberOfDimensions==3 )  
                {
                    OV_ABORT("finish me: 3D upwinding");
                }   

        // int extra=-1; 
                int extra=0; 
                getIndex(mg.gridIndexRange(),I1,I2,I3,extra); // add dissipation here  

                if( isRectangular )
                {
                        FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the grid
                        {
                            if( maskLocal(i1,i2,i3)>0 )
                            {
                                for( int eq = eq1; eq<=eq2; eq++ )
                                {
                  // The upwinding has a sqrt(-1) multiplying it and so it fits in the other component: 
                                    const int uc = eq==eq1 ? uic : urc;
                  // Real and Imaginary equations have different signs for the dissipation term: 
                                    const Real eqnSign = eq==eq1 ? -1. : 1.; 
                  // const Real eqnSign = eq==eq1 ? 1. : -1.; 
                  // Range MM = eq==eq1 ? M0 : M1;
                                    ForStencil(m1,m2,m3)
                                    {
                                        if( m1!=0 && m2!=0 ) // upwind stencil is "thin"
                                            continue;
                                        int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                        int mm = M123CE(m1,m2,m3,uc,eq);  // the system coeff-index 
                    // // This could be cleaned up using an index array -- do this for now 
                    // if( (m1==0 && abs(m2) > halfWidth2)  || ( m2==0 && abs(m1) >halfWidth1 ) )
                    // {
                    //   // --- upwind stencil lies outside usual stencil ---
                    //   // Choose an index past the last normal index 
                    //   // There are 4 extra coefficients in 2D for EACH equation
                    //   //                    X 4
                    //   //                    |
                    //   //                    O
                    //   //                    |
                    //   //            X---O---O---O---X         
                    //   //            1       |       2
                    //   //                    O
                    //   //                    |
                    //   //                    X 3
                    //   // 
                    //   mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq);  // index of last entry of normal stencil
                    //   mm += 4*(eq-eq1);   // offset for equation 
                    //   if( m1==-2 ) 
                    //     mm += 1;
                    //   else if( m1==2 ) 
                    //     mm += 2;
                    //   else if( m2==-2 ) 
                    //     mm += 3;
                    //   else if( m2==2 ) 
                    //     mm += 4;  
                    //   else
                    //   {
                    //     OV_ABORT("unexpected upwind stencil");
                    //   }                    
                    // }
                                        if( m2==0 )
                                            coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction 
                                        if( m1==0 )
                                            coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction
                    // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
                    // NOTE: This is probably not needed: 
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber
                                    }
                  // add extra entries for wider upwind stencil:
                  //
                  // --- upwind stencil lies outside usual stencil ---
                  // Choose an index past the last normal index 
                  // There are 4 extra coefficients in 2D for EACH equation
                  //                    X 4
                  //                    |
                  //                    O
                  //                    |
                  //            X---O---O---O---X         
                  //            1       |       2
                  //                    O
                  //                    |
                  //                    X 3
                  // 
                                    if( upwindWidth==5 && isRectangular )
                                    {
                    // Put extra entries in UNUSED SPOTS
                    // ** This works for Cartesian Grids ***
                                        int m1=-2, m2=0, m3=0;
                                        int mm = M123CE(-1,-1,0,uc,eq);  
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction 
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber    
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=2; m2=0; 
                                        mm = M123CE(1,-1,0,uc,eq);
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,-1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(-1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (-1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,-1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=-2;
                                        mm = M123CE(-1,1,0,uc,eq);
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=2;
                                        mm = M123CE(1,1,0,uc,eq);
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction  
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber              
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,-1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(-1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,-1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,-1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                    if( upwindWidth==5 && !isRectangular )
                                    {
                    // ---- Curvilinear grid : we need extra points in the stencil ---
                                        int m1=-2, m2=0, m3=0;
                                        int mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+1;                   // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction 
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=2; m2=0; 
                                        mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+2;                       // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,-1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(-1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (-1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,-1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=-2;
                                        mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+3;                       // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=2;
                                        mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+4;                       // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction  
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber              
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,-1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(-1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,-1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,-1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                    }      
                                } 
                            } // end if maskLocal
                        } // end for3d
                }
                else
                {
                        FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the grid
                        {
                            if( maskLocal(i1,i2,i3)>0 )
                            {
                  // eval upwind coeff for curvilinear grids 
                                    for( int dir=0; dir<numberOfDimensions; dir++ )
                                        upwindCoefficient[dir] = betaUpwind*sqrt( SQR(RXLocal(i1,i2,i3,dir,0)) + SQR(RXLocal(i1,i2,i3,dir,1)) )/dr[dir];
                                for( int eq = eq1; eq<=eq2; eq++ )
                                {
                  // The upwinding has a sqrt(-1) multiplying it and so it fits in the other component: 
                                    const int uc = eq==eq1 ? uic : urc;
                  // Real and Imaginary equations have different signs for the dissipation term: 
                                    const Real eqnSign = eq==eq1 ? -1. : 1.; 
                  // const Real eqnSign = eq==eq1 ? 1. : -1.; 
                  // Range MM = eq==eq1 ? M0 : M1;
                                    ForStencil(m1,m2,m3)
                                    {
                                        if( m1!=0 && m2!=0 ) // upwind stencil is "thin"
                                            continue;
                                        int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                        int mm = M123CE(m1,m2,m3,uc,eq);  // the system coeff-index 
                    // // This could be cleaned up using an index array -- do this for now 
                    // if( (m1==0 && abs(m2) > halfWidth2)  || ( m2==0 && abs(m1) >halfWidth1 ) )
                    // {
                    //   // --- upwind stencil lies outside usual stencil ---
                    //   // Choose an index past the last normal index 
                    //   // There are 4 extra coefficients in 2D for EACH equation
                    //   //                    X 4
                    //   //                    |
                    //   //                    O
                    //   //                    |
                    //   //            X---O---O---O---X         
                    //   //            1       |       2
                    //   //                    O
                    //   //                    |
                    //   //                    X 3
                    //   // 
                    //   mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq);  // index of last entry of normal stencil
                    //   mm += 4*(eq-eq1);   // offset for equation 
                    //   if( m1==-2 ) 
                    //     mm += 1;
                    //   else if( m1==2 ) 
                    //     mm += 2;
                    //   else if( m2==-2 ) 
                    //     mm += 3;
                    //   else if( m2==2 ) 
                    //     mm += 4;  
                    //   else
                    //   {
                    //     OV_ABORT("unexpected upwind stencil");
                    //   }                    
                    // }
                                        if( m2==0 )
                                            coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction 
                                        if( m1==0 )
                                            coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction
                    // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
                    // NOTE: This is probably not needed: 
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber
                                    }
                  // add extra entries for wider upwind stencil:
                  //
                  // --- upwind stencil lies outside usual stencil ---
                  // Choose an index past the last normal index 
                  // There are 4 extra coefficients in 2D for EACH equation
                  //                    X 4
                  //                    |
                  //                    O
                  //                    |
                  //            X---O---O---O---X         
                  //            1       |       2
                  //                    O
                  //                    |
                  //                    X 3
                  // 
                                    if( upwindWidth==5 && isRectangular )
                                    {
                    // Put extra entries in UNUSED SPOTS
                    // ** This works for Cartesian Grids ***
                                        int m1=-2, m2=0, m3=0;
                                        int mm = M123CE(-1,-1,0,uc,eq);  
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction 
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber    
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=2; m2=0; 
                                        mm = M123CE(1,-1,0,uc,eq);
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,-1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(-1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (-1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,-1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=-2;
                                        mm = M123CE(-1,1,0,uc,eq);
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=2;
                                        mm = M123CE(1,1,0,uc,eq);
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction  
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber              
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,-1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(-1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,-1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,-1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                    if( upwindWidth==5 && !isRectangular )
                                    {
                    // ---- Curvilinear grid : we need extra points in the stencil ---
                                        int m1=-2, m2=0, m3=0;
                                        int mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+1;                   // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction 
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=2; m2=0; 
                                        mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+2;                       // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m1)*upwindCoefficient[0];   // x-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,-1,0,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(-1), j2e=j2 + ie*(0), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (-1,0,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,-1,0,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=-2;
                                        mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+3;                       // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber          
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                        m1=0; m2=2;
                                        mm = M123CE(halfWidth1,halfWidth2,halfWidth3,uc,eq)+4;                       // *NOTE* special formula for location of coeff
                                        coeffLocal(mm,i1,i2,i3) += eqnSign*upwindWeights(m2)*upwindCoefficient[1];   // y-direction  
                                        j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                        setEquationNumber(mm, eq,i1,i2,i3,  uc,j1,j2,j3 );  // macro to set equationNumber              
                                        if( maskLocal(j1,j2,j3)==0 )  
                                            {
                                                assert( maskLocal(j1,j2,j3)==0 );
                                                if( true || classify(j1,j2,j3,eq) !=  SparseRepForMGF::active ) // this point may have already been extrapolated
                                                {
                                                    const int ice = eq; // interpolate this component
                                                    if( false )
                                                        printF("solveHelmholtzDirect: EXTRAPOLATE upwind stencil point: grid=%d, EQ=%d, (j1,j2,j3)=(%d,%d,%d) direction=(%d,%d,%d)\n",grid,eq,j1,j2,j3,0,-1,0);
                          // turn the unused point into an active point (new option added to SparseRep, July 10, 2023)
                                                    setClassify( eq,j1,j2,j3, SparseRepForMGF::active );   
                                                    Range Me = M0 + CE(0,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                                                    Me = M0 + CE(1,eq);
                                                    coeffLocal(Me,j1,j2,j3)=0.; 
                          // int mStart = M0.getBase() + CE(0,eq);
                          // Range Me(mStart,mStart+stencilSize -1); 
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // Range Me = MM;
                          // coeffLocal(Me,j1,j2,j3) = 0.0;         // zero all coeff to start
                          // --- fill in the coefficients of the extrapolation formula ---
                                                    for( int ie=0; ie<=extrapOrder; ie++ )
                                                    {
                                                        int m0 = ie + CE(ice,eq);  // start location in matrix for equation=eq and component=ice is CE(ice,eq)
                                                        coeffLocal(m0,j1,j2,j3) = extrapCoeff[ie];
                                                        int j1e=j1 + ie*(0), j2e=j2 + ie*(-1), j3e=j3 + ie*(0);     // index of point "m" in extrap formula is shifted in the direction (0,-1,0)
                                                        setEquationNumber(m0, eq,j1,j2,j3,  ice,j1e,j2e,j3e );      // macro to set equationNumber
                                                        if( ie>0 && maskLocal(j1e,j2e,j3e)==0 )
                                                        {
                                                            printF("solveHelmHoltzDirect: ERROR: extrapolation formula invovlves an unused point -- FIX ME: reduce order of extrap\n");
                                                            printF(" extrap unused point (j1,j2,j3)=(%d,%d,%d), direction=(%d,%d,%d), grid=%d, ie=%d, point (j1,j2,j3)=(%d,%d,%d)\n",
                                                                        j1,j2,j3,0,-1,0,grid,ie,j1e,j2e,j3e);
                                                            OV_ABORT("error");
                                                        }
                                                    }
                                                }
                                            }
                                    }      
                                } 
                            } // end if maskLocal
                        } // end for3d
                }


            }     
        } // end if upwind

    // --------------- FILL BOUNDARY CONDITIONS ----------------

        if( true )
        {
      // *** NEW WAY : July 6, 2023 ****

      // OV_ABORT("FINISH ME");

      // --- This mostly duplicates code from implicit.bC ---
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const bool isRectangular = mg.isRectangular();

                realMappedGridFunction & coeffg = coeff[grid];
                MappedGridOperators & mgop = *coeffg.getOperators();
                
                // SparseRepForMGF & sparse = *coeffg.sparse;
                // const int numberOfComponentsForCoefficients = sparse.numberOfComponents;  // size of the system of equations
                // const int numberOfGhostLines                = sparse.numberOfGhostLines;
                // const int stencilSize                       = sparse.stencilSize;
                // const int stencilDim=stencilSize*numberOfComponentsForCoefficients; // number of coefficients per equation
                // const int equationOffset      = sparse.equationOffset;
                // IntegerArray & equationNumber = sparse.equationNumber;
                // IntegerArray & classify       = sparse.classify;
                // const int equationNumberBase1  =equationNumber.getBase(1);
                // const int equationNumberLength1=equationNumber.getLength(1);
                // const int equationNumberBase2  =equationNumber.getBase(2);
                // const int equationNumberLength2=equationNumber.getLength(2);
                // const int equationNumberBase3  =equationNumber.getBase(3);
                // const int orderOfAccuracy=mgop.getOrderOfAccuracy(); assert( orderOfAccuracy==2 );
                // // stencil width's and half-width's : 
                // // const int width = orderOfAccuracy+1;
                // const int width = stencilWidth;
                // const int halfWidth1 = (width-1)/2;
                // const int halfWidth2 = numberOfDimensions>1 ? halfWidth1 : 0;
                // const int halfWidth3 = numberOfDimensions>2 ? halfWidth1 : 0;
                    assert( coeffg.sparse!=NULL );
                    SparseRepForMGF & sparse = *coeffg.sparse;
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
                    Range M0 = stencilSize;  
                    Range M = coeffg.dimension(0);
                    Range M1 = M0 + CE(1,1); // for second equation
                              
        // assert( coeffg.sparse!=NULL );

        // SparseRepForMGF & sparse = *coeffg.sparse;
        // const int numberOfComponentsForCoefficients = sparse.numberOfComponents;  // size of the system of equations
        // const int numberOfGhostLines                = sparse.numberOfGhostLines;
        // const int stencilSize                       = sparse.stencilSize;
        // const int stencilDim=stencilSize*numberOfComponentsForCoefficients; // number of coefficients per equation

            
        // const int equationOffset      = sparse.equationOffset;
        // IntegerArray & equationNumber = sparse.equationNumber;
        // IntegerArray & classify       = sparse.classify;

        // const int equationNumberBase1  =equationNumber.getBase(1);
        // const int equationNumberLength1=equationNumber.getLength(1);
        // const int equationNumberBase2  =equationNumber.getBase(2);
        // const int equationNumberLength2=equationNumber.getLength(2);
        // const int equationNumberBase3  =equationNumber.getBase(3);

        // const int orderOfAccuracy=mgop.getOrderOfAccuracy(); assert( orderOfAccuracy==2 );
                
        // // stencil width's and half-width's : 
        // // const int width = orderOfAccuracy+1;
        // const int width = stencilWidth;

        // const int halfWidth1 = (width-1)/2;
        // const int halfWidth2 = numberOfDimensions>1 ? halfWidth1 : 0;
        // const int halfWidth3 = numberOfDimensions>2 ? halfWidth1 : 0;


        // Range M0 = stencilSize;  // for first equation
        // Range M1 = M0 + CE(1,1); // for second equation
        // Range M = coeffg.dimension(0);

                OV_GET_SERIAL_ARRAY(Real,coeffg,coeffLocal);

                    Real dx[3]={1.,1.,1.};
                    Real dr[3]={1.,1.,1.};
                    if( isRectangular )
                    { // rectangular grid grid-spacings: 
                        mg.getDeltaX(dx);
                    }
                    else
                    {
            // unit square grid spacings: 
                        for( int dir=0; dir<3; dir++ )
                            dr[dir]=mg.gridSpacing(dir);           
                    }        



                const int e=0, ic=0; // equation number and component number 

                const int mDiag0 = M123CE(0,0,0,ic1,eq1);   // M123(0,0,0);              // index of diagonal entry, equation 0 
                const int mDiag1 = M123CE(0,0,0,ic2,eq2);    // index of diagonal entry, equation 1 



                ForBoundary(side,axis)
                {

                    Index Ib1,Ib2,Ib3;
                    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);

          // Set the index-shift for this side
                    is1=is2=is3=0;
                    isv[axis]=1-2*side;   // +1 on left and -1 on right  

                    const int axisp1 = (axis+1) % numberOfDimensions;    

                    if( mg.boundaryCondition(side,axis)==CgWave::dirichlet         ||
              // mg.boundaryCondition(side,axis)==CgWave::absorbing ||  // ** DO THIS FOR NOW : absorbing terminated with Dirichlet
                            mg.boundaryCondition(side,axis)==CgWave::exactBC )
                    {
            // ------------ FILL DIRICHLET BC ------------

                        printF("+++++ solveHelmholtzDirect BC: FILL MATRIX BC=%d FOR (grid,side,axis)=(%d,%d,%d) DIRICHLET/ABSORBING/EXACT\n",mg.boundaryCondition(side,axis),grid,side,axis);
                        coeffLocal(     M,Ib1,Ib2,Ib3) = 0.0;  // zero out any existing equations
                        coeffLocal(mDiag0,Ib1,Ib2,Ib3) = 1.0;            
                        coeffLocal(mDiag1,Ib1,Ib2,Ib3) = 1.0;

            // --- EXTRAPOLATE GHOST LINES ---

                        for( int ghost=1; ghost<=numberOfGhostLines; ghost++ )
                        {
                            {
                                getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                                coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0; // zero all coeff to start
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                  // --- fill in the coefficients of the extrapolation formula ---
                                    for( int me=0; me<=extrapOrder; me++ )
                                    {
                                        int m0 = me + M0.getBase(); 
                                        coeffLocal(m0,i1m,i2m,i3m) = extrapCoeff[me];
                                        int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                        setEquationNumber(m0, eq1,i1m,i2m,i3m,  ic1,j1,j2,j3 );      // macro to set equationNumber
                                    }
                                } // end FOR_3D
                            }
                            {
                                getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                                coeffLocal(M1,Ig1,Ig2,Ig3) = 0.0; // zero all coeff to start
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                  // --- fill in the coefficients of the extrapolation formula ---
                                    for( int me=0; me<=extrapOrder; me++ )
                                    {
                                        int m0 = me + M1.getBase(); 
                                        coeffLocal(m0,i1m,i2m,i3m) = extrapCoeff[me];
                                        int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                        setEquationNumber(m0, eq2,i1m,i2m,i3m,  ic2,j1,j2,j3 );      // macro to set equationNumber
                                    }
                                } // end FOR_3D
                            }
                        } // end for ghost

                    }
                    else if( mg.boundaryCondition(side,axis)==CgWave::neumann )
                    {
            // ------------ FILL NEUMANN BC ------------

                        printF("+++++ solveHelmholtzDirect BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) NEUMANN\n",grid,side,axis);
                        OV_ABORT("FINISH ME");

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

                // Specify that the above coeff value is the coefficient of component ic at the grid point (j1,j2,j3).
                                int j1=i1+m1, j2=i2+m2, j3=i3+m3;                       // the stencil is centred on the boundary pt (i1,i2,i3)
                                setEquationNumber(m, e,i1m,i2m,i3m,  ic,j1,j2,j3 );      // macro to set equationNumber
                            }

                        } // end FOR_3D

            // fill ghost 2 with extrapolation
                        for( int ghost=2; ghost<=numberOfGhostLines; ghost++ )
                        {
                            {
                                getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                                coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0; // zero all coeff to start
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                  // --- fill in the coefficients of the extrapolation formula ---
                                    for( int me=0; me<=extrapOrder; me++ )
                                    {
                                        int m0 = me + M0.getBase(); 
                                        coeffLocal(m0,i1m,i2m,i3m) = extrapCoeff[me];
                                        int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                        setEquationNumber(m0, eq1,i1m,i2m,i3m,  ic1,j1,j2,j3 );      // macro to set equationNumber
                                    }
                                } // end FOR_3D
                            }
                            {
                                getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                                coeffLocal(M1,Ig1,Ig2,Ig3) = 0.0; // zero all coeff to start
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                  // --- fill in the coefficients of the extrapolation formula ---
                                    for( int me=0; me<=extrapOrder; me++ )
                                    {
                                        int m0 = me + M1.getBase(); 
                                        coeffLocal(m0,i1m,i2m,i3m) = extrapCoeff[me];
                                        int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                        setEquationNumber(m0, eq2,i1m,i2m,i3m,  ic2,j1,j2,j3 );      // macro to set equationNumber
                                    }
                                } // end FOR_3D
                            }
                        } // end for ghost


                    }
                    else if( mg.boundaryCondition(side,axis)==CgWave::absorbing || 
                                      mg.boundaryCondition(side,axis)==CgWave::abcEM2 )
                    {
                        printF("+++++ solveHelmholtzDirect: BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) ABSORBING/EM2, c=%e\n",grid,side,axis,c);

            // EM2 
            //     c D+xD-x + .5*c D+yD-y - is*I*sgn*omega* D0x = 0 
            // 
            // Engquist-Majda order2 scheme
            //  D+t (-D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : left 
            //  D+t ( D0x W^n) + D+xD-x .5*(W^{n+1} + W^n ) + .5* D+yD-y .5*(W^{n+1} + W^n ) = 0 : right                
            //
            //   - D0x + c D+xD-x + .5*c D+yD-y  : left 
            //     D0x + c D+xD-x + .5*c D+yD-y  : right

             // *** FINISH ME ***
                        Real ca = c; // c * 0.;
            // if( par.adjustOmega )
            //   ca = c*tan(par.frequencyArray(1)*dt/2)/(par.frequencyArraySaved(1)*dt/2); % adjust c for dt errors 
            // end


            // res = -is*(unx-ucx)/dt + (.5*c)*( unxx + ucxx) + (.25*c)*( unyy + ucyy );

                        Range Rw(-halfWidth1,halfWidth1);
                        RealArray abcCoeffr(Rw,Rw,Rw), abcCoeffi(Rw,Rw,Rw);
                        abcCoeffr=0.; abcCoeffi=0.;

            // om = frequencyArray(freq); 

            // A(ie,ie)=  c/dx^2   + par.isign*1i*om/(2*dx);
            // A(ie,ib)= -2*c/dx^2 -c/dy^2;
            // A(ie,je)=  c/dx^2   - par.isign*1i*om/(2*dx);

            // if( dir==1 )
            //   A(ie,eqn(ix,iy-1)) = .5*c/dy^2; 
            //   A(ie,eqn(ix,iy+1)) = .5*c/dy^2; 
            // else
            //   A(ie,eqn(ix-1,iy)) = .5*c/dy^2; 
            //   A(ie,eqn(ix+1,iy)) = .5*c/dy^2;                  
            // end            


            // EM2: BC: 
            //     -is*(u).xt +  c( Dxx u + .5*Dyy u ) = 0
            // => 
            //   - is* i * omega*omegaSign* Dx u + c( Dxx u + .5*Dyy u ) = 0
            //  
            // u = ur + i*ui : 
            //     is* omega*omegaSign* Dx ui + c( Dxx ur + .5*Dyy ur ) = 0 
            //    -is* omega*omegaSign* Dx ur + c( Dxx ui + .5*Dyy ui ) = 0 
            //             
                        if( isRectangular )
                        {
                            const Real dxa2 = dx[axis]*dx[axis];
                            const Real dya2 = dx[axisp1]*dx[axisp1];
                            abcCoeffr(-is1,-is2,0) = (    ca/dxa2           ); //   + 1./(2.*dx[axis]*dt);    // ghost
                            abcCoeffr(   0,   0,0) = (-2.*ca/dxa2 - ca/dya2 );                           // boundary 
                            abcCoeffr(+is1,+is2,0) = (    ca/dxa2           ); //   - 1./(2.*dx[axis]*dt);    // first line in

              // ** CHECK sign ***
                            abcCoeffi(-is1,-is2,0) =  + omegaSign*omega/(2.*dx[axis]);    // ghost
                            abcCoeffi(+is1,+is2,0) =  - omegaSign*omega/(2.*dx[axis]);    // first line in

                            if( axis==0 )
                            {
                                abcCoeffr(0,-1,0) = ( .5*ca/(dx[1]*dx[1]) ); 
                                abcCoeffr(0,+1,0) = ( .5*ca/(dx[1]*dx[1]) ); 
                            }
                            else
                            {
                                abcCoeffr(-1,0,0) = ( .5*ca/(dx[0]*dx[0]) ); 
                                abcCoeffr(+1,0,0) = ( .5*ca/(dx[0]*dx[0]) );                   
                            }     
                        }
                        else
                        {
                              OV_ABORT("solveHelmholtzDirect:: IMP BC ABORBING -- FINISH ME: curvilinear"); 
                        }

                        if( false )
                        {
                            ::display(abcCoeffr,"abcCoeffr");
                            ::display(abcCoeffi,"abcCoeffi");
                        }
                        FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                        {
                                
                            int i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)

              // Specify that this a "real" equation on the first ghost line: 
              // (A "real" equation has a possible non-zero right-hand-side)
                            setClassify(eq1,i1m,i2m,i3m, SparseRepForMGF::ghost1);              
                            setClassify(eq2,i1m,i2m,i3m, SparseRepForMGF::ghost1);              

                            ForStencil(m1,m2,m3)
                            {
                // -- Equation 1 : Real part 
                //     is* omega*omegaSign* Dx ui + c( Dxx ur + .5*Dyy ur ) = 0 
                            
                                int mr  = M123CE(m1,m2,m3,ic1,eq1);                        // equation 1, component 1
                                coeffLocal(mr,i1m,i2m,i3m) = abcCoeffr(m1,m2,m3);
                                int j1=i1+m1, j2=i2+m2, j3=i3+m3;                          // the stencil is centred on the boundary pt (i1,i2,i3)
                                setEquationNumber(mr, eq1,i1m,i2m,i3m,  ic1,j1,j2,j3 );    

                                int mi  = M123CE(m1,m2,m3,ic2,eq1);                        // equation 1, component 2
                                coeffLocal(mi,i1m,i2m,i3m) = -abcCoeffi(m1,m2,m3);         // *CHECK SIGN 
                                setEquationNumber(mi, eq1,i1m,i2m,i3m,  ic2,j1,j2,j3 );    

                // -- Equation 2 : imaginary part: 
                //    -is* omega*omegaSign* Dx ur + c( Dxx ui + .5*Dyy ui ) = 0   
                                mi  = M123CE(m1,m2,m3,ic2,eq2);                            // equation 2, component 2
                                coeffLocal(mi,i1m,i2m,i3m) =  abcCoeffr(m1,m2,m3);
                                setEquationNumber(mi, eq2,i1m,i2m,i3m,  ic2,j1,j2,j3 );    

                                mr  = M123CE(m1,m2,m3,ic1,eq2);                            // equation 2, component 1
                                coeffLocal(mr,i1m,i2m,i3m) = +abcCoeffi(m1,m2,m3);
                                setEquationNumber(mr, eq2,i1m,i2m,i3m,  ic1,j1,j2,j3 );    

                            }



                        } // end FOR_3D

            // fill ghost 2 with extrapolation
                        for( int ghost=2; ghost<=numberOfGhostLines; ghost++ )
                        {
                            {
                                getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                                coeffLocal(M0,Ig1,Ig2,Ig3) = 0.0; // zero all coeff to start
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                  // --- fill in the coefficients of the extrapolation formula ---
                                    for( int me=0; me<=extrapOrder; me++ )
                                    {
                                        int m0 = me + M0.getBase(); 
                                        coeffLocal(m0,i1m,i2m,i3m) = extrapCoeff[me];
                                        int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                        setEquationNumber(m0, eq1,i1m,i2m,i3m,  ic1,j1,j2,j3 );      // macro to set equationNumber
                                    }
                                } // end FOR_3D
                            }
                            {
                                getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,ghost);
                                coeffLocal(M1,Ig1,Ig2,Ig3) = 0.0; // zero all coeff to start
                                FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                                {
                                    int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                  // --- fill in the coefficients of the extrapolation formula ---
                                    for( int me=0; me<=extrapOrder; me++ )
                                    {
                                        int m0 = me + M1.getBase(); 
                                        coeffLocal(m0,i1m,i2m,i3m) = extrapCoeff[me];
                                        int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                        setEquationNumber(m0, eq2,i1m,i2m,i3m,  ic2,j1,j2,j3 );      // macro to set equationNumber
                                    }
                                } // end FOR_3D
                            }
                        } // end for ghost

                    }
                    else if(  mg.boundaryCondition(side,axis)> 0 )
                    {
                        printF("solveHelmholtzDirect::ERROR: unknown boundaryCondition=%d \n",mg.boundaryCondition(side,axis));
                        OV_ABORT("error");

                    }
                }
            }  // end for grid

        }
        else
        {
      // -- OLD WAY: 

      // Here are the boundary conditions we can impose with the high-level operators: 
            coeff.applyBoundaryConditionCoefficients(eq1,urc,BCTypes::dirichlet,  BCTypes::allBoundaries);
            coeff.applyBoundaryConditionCoefficients(eq1,urc,BCTypes::extrapolate,BCTypes::allBoundaries);

            coeff.applyBoundaryConditionCoefficients(eq2,uic,BCTypes::dirichlet,  BCTypes::allBoundaries);
            coeff.applyBoundaryConditionCoefficients(eq2,uic,BCTypes::extrapolate,BCTypes::allBoundaries);    

        }


        coeff.finishBoundaryConditions();    

    // ---- FILL IN RHS ---- 
    // OV_ABORT("**FINISH ME");

        Real amp,kx,ky,kz,phi; // for plane wave incident field
        if( solveForScatteredField && knownSolutionOption==planeWave )
        {
      // Get planeWave parameters: 
            amp   = cgWave.dbase.get<Real>("ampPlaneWave");
            kx    = cgWave.dbase.get<Real>("kxPlaneWave");
            ky    = cgWave.dbase.get<Real>("kyPlaneWave");
            kz    = cgWave.dbase.get<Real>("kzPlaneWave");
            phi   = cgWave.dbase.get<Real>("phiPlaneWave");

            Real omegaPlaneWave = cgWave.dbase.get<Real>("omegaPlaneWave");
            printF("\n solveHelmholtzDirect: omega=%14.6e, omegaPlaneWave=%14.6e\n",omega,omegaPlaneWave);
        }

    // Index I1,I2,I3; 
    // Index Ib1,Ib2,Ib3; 
    // Index Ig1,Ig2,Ig3; 
        for( int grid=0; grid<numberOfComponentGrids; grid++ )
        {
            MappedGrid & mg = cg[grid];
            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);
            OV_GET_SERIAL_ARRAY(Real,fc[grid],fcLocal);

            getIndex(mg.dimension(),I1,I2,I3);

      // -- make sure the RHS is zero at active points (unused points which are extrpolated for upwinding)
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
        // if( isnan(fcLocal(i1,i2,i3,0)) )
        // {
        //   printF("solveHelmholtzDirect: ERROR: fc(%d,%d,%d,0) isnan\n",i1,i2,i3);
        //   OV_ABORT("error");
        // }
                if( maskLocal(i1,i2,i3) > 0 )
                    fcLocal(i1,i2,i3,0)=fLocal(i1,i2,i3);
                else
                    fcLocal(i1,i2,i3,0)=0.;
            }
      // fcLocal(I1,I2,I3,0) = fLocal(I1,I2,I3);
            fcLocal(I1,I2,I3,1) = 0.;

      // --- ASSIGN SPECIAL boundary conditions --- 
      //  Boundary conditions will be filled in by getHelmholtzForcing, but these are not all correct for the complex case

            ForBoundary(side,axis)
            {
                if( mg.boundaryCondition(side,axis)==CgWave::dirichlet )
                {
                    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
          // OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);

                    if( solveForScatteredField && knownSolutionOption==planeWave )
                    {
                        OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);

            // The INCIDENT field is 
            //    exp( i ( k.x - omega*t + phi ))
            //    phi = -pi/2 to change cos() into sin(), and sin into -cos() 
            //      --> The incident field for the wave solver is sin( k.x - omega*t) for some reason

            // we substract off the incident field of amp*[ sin( k.x+ phi ) - I cos( k.x + phi ) ]
                        if( mg.numberOfDimensions()==2 )
                        { // *check me*
                            fcLocal(Ib1,Ib2,Ib3,0)= -amp*sin( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + phi );
                            fcLocal(Ib1,Ib2,Ib3,1)=  amp*cos( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + phi );
                        }
                        else
                        {
                            fcLocal(Ib1,Ib2,Ib3,0)= -amp*sin( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + kz*xLocal(Ib1,Ib2,Ib3,2) + phi );
                            fcLocal(Ib1,Ib2,Ib3,1)= +amp*cos( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + kz*xLocal(Ib1,Ib2,Ib3,2) + phi );                  
                        }
                    }
                }
            } // end ForBoundary


      // 
      // *** CHECK that active points have a zero RHS ****
      // 
            bool checkActive=true; 
            if( checkActive && upwind )
            {
                SparseRepForMGF & sparse = *coeff[grid].sparse;
                intArray & classify = sparse.classify;
                OV_GET_SERIAL_ARRAY(int,classify,classifyLocal);
                int numActive=0;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {  
                    if( classifyLocal(i1,i2,i3,0)==SparseRepForMGF::active  )
                    {
                        numActive++;
                        if( fcLocal(i1,i2,i3,0)!=0. )  
                        {
                            printF("solveHelmDirect:ERROR: RHS fcLocal is NOT ZERO at an active point (unused points which are extrapolated for upwinding)\n");
                            OV_ABORT("ERROR");
                        } 
                    }   
                }
                printF("INFO: There are %d active points on grid=%d. (active points = unused points which are extrapolated for upwinding)\n",numActive,grid);

            }

      // if( solveForScatteredField )
      //   ::display(fcLocal,"solveHelmholtzDirect: fcLocal","%9.2e ");

        }


        Oges solver( cg );                     // create a solver
        solver.setCoefficientArray( coeff );   // supply coefficients

        if( true )
        {
                if( dbase.has_key("helmholtzOgesParameters") )
                {  
                    printf("solveHelmholtzDirect: OgesParameters found. Changing the direct Helmholtz solver parameters. \n");
                    OgesParameters & par = dbase.get<OgesParameters>("helmholtzOgesParameters");
                    solver.setOgesParameters(par);
                }
                else
                {
                    int solverType=OgesParameters::yale; 
          // solverType=OgesParameters::PETSc;
          // solverType=OgesParameters::PETScNew; // parallel
                    solver.set(OgesParameters::THEsolverType,solverType); 
                    if( solverType==OgesParameters::PETSc )
                      solver.set(OgesParameters::THEsolverMethod,OgesParameters::biConjugateGradientStabilized);
          // solver.set(OgesParameters::THEparallelSolverMethod,OgesParameters::gmres);
          // solver.set(OgesParameters::THErelativeTolerance,max(tol,REAL_EPSILON*10.));
          // solver.set(OgesParameters::THEmaximumNumberOfIterations,10000);
          // if( iluLevels>=0 )
          //   solver.set(OgesParameters::THEnumberOfIncompleteLULevels,iluLevels);
                }
                if( upwind )
                {
                    Real fillinRatio=0;
          //solver.get(OgesParameters::THEfillinRatio,fillinRatio );
          // printF("Current fillinRatio = %g\n",fillinRatio);
                    int stencilWidth = orderOfAccuracy+1;
                    int stencilSize=int( pow(stencilWidth,cg.numberOfDimensions())+1 );  // add 1 for interpolation equations
                    fillinRatio = (stencilSize+2)*2 + 20.; // what should this be ?
                    printF("Set new fillinRatio = %g\n",fillinRatio);
                    solver.set(OgesParameters::THEfillinRatio,fillinRatio );
                }
        }

        helmholtzSolverName = solver.parameters.getSolverName();
        printF("\n === Direct Helmholtz Solver (Complex) :\n %s\n =====\n",(const char*)helmholtzSolverName);    

        u=0.;  // initial guess for iterative solvers
        solver.solve( u,fc ); 

    // OV_ABORT("Solve the complex Helmholtz : FINISH ME");


        bool checkResidual=true;
        if( checkResidual )
        {
            realCompositeGridFunction & res = dbase.get<realCompositeGridFunction>("residual");
            Index Ib1,Ib2,Ib3;
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const bool isRectangular = mg.isRectangular();

                MappedGridOperators & mgop = op[grid];

                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                OV_GET_SERIAL_ARRAY(real,fc[grid],fLocal);
                OV_GET_SERIAL_ARRAY(real,res[grid],resLocal);

        // OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rsxyLocal); // for upwinding 

                    bool useOriginalMetrics = useSuperGrid && superGrid(grid) && !isRectangular;
          // OV_GET_SERIAL_ARRAY_CONDITIONAL(Real,mg.inverseVertexDerivative(),rsxyLocal2,!useOriginalMetrics);
                    OV_GET_SERIAL_ARRAY_CONDITIONAL(Real,mg.inverseVertexDerivative(),rsxyLocalNew,!useOriginalMetrics);
          // if( useSuperGrid && superGrid(grid) )
          // {
          //   // -- SuperGrid changes the metrics but upwinding uses the original metrics: 
          //   RealArray *& rxOriginal = cgWave.dbase.get<RealArray*>("rxOriginal");
          //   // rsxyLocal.reference(rxOriginal[grid]);
          // }
                    RealArray & rsxyLocal = useOriginalMetrics ? cgWave.dbase.get<RealArray*>("rxOriginal")[grid] : rsxyLocalNew;

                int extra=0; // -1;
                getIndex(cg[grid].indexRange(),I1,I2,I3,extra); 
                
                bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
                if( ok )
                {
                    resLocal=0.;
                    
                    RealArray lap(I1,I2,I3,2);
                    op[grid].derivative(MappedGridOperators::laplacianOperator,uLocal,lap,I1,I2,I3);

          // extra=-1; 
          // extra=0; 
          // getIndex(cg[grid].gridIndexRange(),I1,I2,I3,extra); 

                    if( true && useSuperGrid )
                    {
                        Real lapNorm = max(fabs(lap(I1,I2,I3,all)));
                        printF("checkResidual: grid=%d, lapNorm=%9.2e\n",grid,lapNorm);
                    }

                    if( !superGrid(grid)  || !isRectangular )
                    {
                        resLocal(I1,I2,I3,0) = cSq*lap(I1,I2,I3,0) + (omega*omega)*uLocal(I1,I2,I3,0) + (sigma*omega*damp)*uLocal(I1,I2,I3,1) - fLocal(I1,I2,I3,0); 
                        resLocal(I1,I2,I3,1) = cSq*lap(I1,I2,I3,1) + (omega*omega)*uLocal(I1,I2,I3,1) - (sigma*omega*damp)*uLocal(I1,I2,I3,0); 
                    }
                    else
                    {
            // ----- rectangular grid + superGrid ------
                        printF("Check residual for Cartesian grid + superGrid : grid=%d\n",grid);

                        assert( superGrid(grid) && isRectangular );

                        resLocal=0.;

            // useAbsorbingLayer(axis,grid) = 1 if this axis has a superGridLayer 
                        IntegerArray & useAbsorbingLayer = cgWave.dbase.get<IntegerArray>("useAbsorbingLayer");

                        getIndex(mg.dimension(),I1,I2,I3);
                        Range R2 = numberOfComponents;
                        RealArray ddDeriv(I1,I2,I3,R2), dDeriv(I1,I2,I3,R2);

                        getIndex(mg.indexRange(),I1,I2,I3,extra);
            // --- Eval xx and x derivatives ---
                        mgop.derivative( MappedGridOperators::xxDerivative,uLocal,ddDeriv,I1,I2,I3,R2);
                        mgop.derivative( MappedGridOperators::xDerivative, uLocal, dDeriv,I1,I2,I3,R2);

                        if( useAbsorbingLayer(0,grid) )
                        {
              // -- scale coefficients using superGrid functions --
                            RealArray *& etaxSuperGrid = cgWave.dbase.get<RealArray*>("etaxSuperGrid" );
                            RealArray & etax = etaxSuperGrid[grid];
                            for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )
                            {
                                ddDeriv(i1,I2,I3,R2) *= etax(i1,0);  // scale by "(r.x)^2"
                                  dDeriv(i1,I2,I3,R2) *= etax(i1,1);  // scale by "r.xx"
                            }
                        }

                        resLocal(I1,I2,I3,R2) = cSq*( ddDeriv(I1,I2,I3,R2) + dDeriv(I1,I2,I3,R2) );    // transformed xx derivative

            // --- Eval yy and y derivatives ---
                        mgop.derivative( MappedGridOperators::yyDerivative,uLocal,ddDeriv,I1,I2,I3,R2);
                        mgop.derivative( MappedGridOperators::yDerivative, uLocal, dDeriv,I1,I2,I3,R2);

                        if( useAbsorbingLayer(1,grid) )
                        {
              // -- scale coefficients using superGrid functions --  
                            RealArray *& etaySuperGrid = cgWave.dbase.get<RealArray*>("etaySuperGrid" );
                            RealArray & etay = etaySuperGrid[grid];
                            for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )
                            {
                                ddDeriv(I1,i2,I3,R2) *= etay(i2,0);  
                                  dDeriv(I1,i2,I3,R2) *= etay(i2,1);  
                            }
                        }
                        resLocal(I1,I2,I3,R2) += cSq*( ddDeriv(I1,I2,I3,R2) + dDeriv(I1,I2,I3,R2) );   

                        resLocal(I1,I2,I3,R2) += (omega*omega)*uLocal(I1,I2,I3,R2);  // omega^2 u 

                        if( damp!=0 )
                        { // ---- damping terms ---
                            resLocal(I1,I2,I3,0) += b1*uLocal(I1,I2,I3,1);  // b1*ui in ur eqn
                            resLocal(I1,I2,I3,1) += b2*uLocal(I1,I2,I3,0);  // b2*ur in ui eqn
                        }

                        resLocal(I1,I2,I3,0) -= fLocal(I1,I2,I3,0); 

                    } 

                    if( upwind ) 
                    {
                            Real dx[3]={1.,1.,1.};
                            Real dr[3]={1.,1.,1.};
                            if( isRectangular )
                            { // rectangular grid grid-spacings: 
                                mg.getDeltaX(dx);
                            }
                            else
                            {
                // unit square grid spacings: 
                                for( int dir=0; dir<3; dir++ )
                                    dr[dir]=mg.gridSpacing(dir);           
                            }        

              // Here is the upwind prefactor : 
              //    In 2D: betaUpwind = (sigma*omega*c)/( sqrt(2)* 8 )
                            const Real dtUpwind=1.;             // dt to use when computing upwind coeff
                            const bool adjustForTimeStep=false; // do not adjust for dt
                            Real upwindDissipationCoefficient = cgWave.getUpwindDissipationCoefficient( grid,dtUpwind,adjustForTimeStep );
                            if( true )
                            {
                                Real adSosup = (c)/( sqrt(1.*numberOfDimensions) * pow(2.,orderOfAccuracy+1) );
                                printF(">>>>>>solveHelmholtzDirect: upwindDissipationCoefficient=%g, old-way=%g, omega=%g\n",upwindDissipationCoefficient,adSosup,omega);
                            }
              // Real betaUpwind = (omegaSign*omega*c)/( sqrt(1.*numberOfDimensions) * pow(2.,orderOfAccuracy+1) ); // OLD
                            Real betaUpwind = (omegaSign*omega)*upwindDissipationCoefficient;
              // betaUpwind *= 10; // ** TEST ***
              // betaUpwind =0.; // *** TEMP TEST
                            if( isRectangular )
                            {
                // --- Here is the upwind coefficient for Cartesian grids ---
                // const Real beta= (1./(8.*sqrt(2.)))*c;
                                upwindCoefficient[0] = betaUpwind/dx[0];  
                                upwindCoefficient[1] = betaUpwind/dx[1];
                            } 

                        FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the grid
                        {
                            if( !isRectangular )
                            {
                                for( int dir=0; dir<numberOfDimensions; dir++ )
                                    upwindCoefficient[dir] = betaUpwind*sqrt( SQR(RXLocal(i1,i2,i3,dir,0)) + SQR(RXLocal(i1,i2,i3,dir,1)) )/dr[dir];
                            }
                            for( int m1=-upwindHalfWidth; m1<=upwindHalfWidth; m1++ )
                            {
                                resLocal(i1,i2,i3,0) += -upwindCoefficient[0]*upwindWeights(m1)*uLocal(i1+m1,i2,i3,1);  // x direction
                                resLocal(i1,i2,i3,0) += -upwindCoefficient[1]*upwindWeights(m1)*uLocal(i1,i2+m1,i3,1);  // y direction

                                resLocal(i1,i2,i3,1) +=  upwindCoefficient[0]*upwindWeights(m1)*uLocal(i1+m1,i2,i3,0);
                                resLocal(i1,i2,i3,1) +=  upwindCoefficient[1]*upwindWeights(m1)*uLocal(i1,i2+m1,i3,0);                
                            }
                        }
                    } 

                    Range all;
                    ForBoundary(side,axis)
                    {
                        int is = 1 -2*side;
            // set residual to zero on dirichlet boundaries 
                        if( mg.boundaryCondition(side,axis) == CgWave::dirichlet )
                        {
                            getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3);
                            if( solveForScatteredField==0 )
                            {
                              resLocal(Ib1,Ib2,Ib3,all)=0.;
                            }
                            else
                            {
                // --- scattering from a Dirichlet BC ----
                                OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
                                if( mg.numberOfDimensions()==2 )
                                { 
                                    resLocal(Ib1,Ib2,Ib3,0) = uLocal(Ib1,Ib2,Ib3,0) - (-amp*sin( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + phi ));
                                    resLocal(Ib1,Ib2,Ib3,1) = uLocal(Ib1,Ib2,Ib3,1) - (+amp*cos( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + phi ));
                                }
                                else
                                {
                                    resLocal(Ib1,Ib2,Ib3,0) = uLocal(Ib1,Ib2,Ib3,0) - (-amp*sin( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + kz*xLocal(Ib1,Ib2,Ib3,2) + phi ));
                                    resLocal(Ib1,Ib2,Ib3,1) = uLocal(Ib1,Ib2,Ib3,1) - (+amp*cos( kx*xLocal(Ib1,Ib2,Ib3,0) + ky*xLocal(Ib1,Ib2,Ib3,1) + kz*xLocal(Ib1,Ib2,Ib3,2) + phi ));
                                }
                                where( maskLocal(Ib1,Ib2,Ib3)<=0 )
                                {
                                    resLocal(Ib1,Ib2,Ib3,0)=0.;
                                    resLocal(Ib1,Ib2,Ib3,1)=0.;
                                }

                                Range R2=2;
                                Real maxResBC1 = max(abs(resLocal(Ib1,Ib2,Ib3,0)));
                                Real maxResBC2 = max(abs(resLocal(Ib1,Ib2,Ib3,1)));

                                printF(" Scattering Dirichlet BC: (side,axis,grid)={%d,%d,%d) max-res=[%9.2e,%9.2e]\n",side,axis,grid,maxResBC1,maxResBC2);
                            }

                          }
                          else if( mg.boundaryCondition(side,axis) == CgWave::absorbing ||
                                            mg.boundaryCondition(side,axis) == CgWave::abcEM2 )
                          {
                                getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3);
                                getGhostIndex(mg.indexRange(),side,axis,Ig1,Ig2,Ig3,1);

                // EM2: BC: 
                //     -is*(u).xt +  c( Dxx u + .5*Dyy u ) = 0
                // => 
                //   - is* i * omega*omegaSign* Dx u + c( Dxx u + .5*Dyy u ) = 0
                //  
                // u = ur + i*ui : 
                //     is* omega*omegaSign* Dx ui + c( Dxx ur + .5*Dyy ur ) = 0 
                //    -is* omega*omegaSign* Dx ur + c( Dxx ui + .5*Dyy ui ) = 0 
                // 
                                Range R2 = numberOfComponents;
                                RealArray ux(Ib1,Ib2,Ib3,R2), uy(Ib1,Ib2,Ib3,R2), uxx(Ib1,Ib2,Ib3,R2), uyy(Ib1,Ib2,Ib3,R2);
                                mgop.derivative(MappedGridOperators::xDerivative,uLocal,ux,Ib1,Ib2,Ib3,R2);
                                mgop.derivative(MappedGridOperators::yDerivative,uLocal,uy,Ib1,Ib2,Ib3,R2);

                                mgop.derivative(MappedGridOperators::xxDerivative,uLocal,uxx,Ib1,Ib2,Ib3,R2);
                                mgop.derivative(MappedGridOperators::yyDerivative,uLocal,uyy,Ib1,Ib2,Ib3,R2);

                                if( axis==0 )
                                {
                                    resLocal(Ig1,Ig2,Ig3,0) = +1.*(-is*omegaSign*omega*ux(Ib1,Ib2,Ib3,1)) - c*( uxx(Ib1,Ib2,Ib3,0) + .5*uyy(Ib1,Ib2,Ib3,0) );
                                    resLocal(Ig1,Ig2,Ig3,1) =     ( is*omegaSign*omega*ux(Ib1,Ib2,Ib3,0)) - c*( uxx(Ib1,Ib2,Ib3,1) + .5*uyy(Ib1,Ib2,Ib3,1) );
                                }
                                else
                                {
                                    resLocal(Ig1,Ig2,Ig3,0) =  1.*( -is*omegaSign*omega*uy(Ib1,Ib2,Ib3,1) ) - c*( uyy(Ib1,Ib2,Ib3,0) + .5*uxx(Ib1,Ib2,Ib3,0) );
                                    resLocal(Ig1,Ig2,Ig3,1) =  1.*(  is*omegaSign*omega*uy(Ib1,Ib2,Ib3,0) ) - c*( uyy(Ib1,Ib2,Ib3,1) + .5*uxx(Ib1,Ib2,Ib3,1) );
                                }

                                where( maskLocal(Ib1,Ib2,Ib3)<=0 )
                                {
                                    resLocal(Ig1,Ig2,Ig3,0)=0.;
                                    resLocal(Ig1,Ig2,Ig3,1)=0.;
                                }
                                Real maxResBC1 = max(abs(resLocal(Ig1,Ig2,Ig3,0)));
                                Real maxResBC2 = max(abs(resLocal(Ig1,Ig2,Ig3,1)));
                                Real maxRHS = max(abs(fLocal(Ig1,Ig2,Ig3,R2)));

                                printF(" Absorbing BC: (side,axis,grid)={%d,%d,%d) max-res=[%9.2e,%9.2e], max|rhs|=%9.2e\n",side,axis,grid,maxResBC1,maxResBC2,maxRHS);

                // ::display(resLocal(Ig1,Ig2,Ig3,0),"resLocal(Ig1,Ig2,Ig3,0)","%8.1e ");
                // ::display(resLocal(Ig1,Ig2,Ig3,1),"resLocal(Ig1,Ig2,Ig3,0)","%8.1e ");

                          }
                    } // end for boundary

          // -- zero residual at unused points ---
                    getIndex(mg.dimension(),I1,I2,I3);
                    where( maskLocal(I1,I2,I3) <=0  )   
                    {
                        resLocal(I1,I2,I3,0)=0.;
                        resLocal(I1,I2,I3,1)=0.;
                    }        

                }

                if( true )
                {
          // ::display(resLocal,"resLocal","%7.0e ");
                    Real maxRes = max(abs(resLocal));
                    if( false && maxRes > .01 )
                    {
                        ::display(resLocal,"resLocal","%7.0e ");
                    }
                    printF("grid=%d : maxRes=%9.2e\n",grid,maxRes);
                }     

            }  // end for grid

            const int maskOption=1;  // maskOption=1 : check points with mask>0
            Real urRes = maxNorm(res,0,maskOption)/(omega*omega);
            Real uiRes = maxNorm(res,1,maskOption)/(omega*omega);

            printF("\n ^^^^^^ solveHelmholtzDirect: max-rel-res [re,Im]=[%8.2e,%8.2e] ^^^^^^^^\n\n",urRes,uiRes);       

        } // end check residual

        bool plotSolution = false; // true; // for debugging we plot the solution

        if( plotSolution )
        {
            GL_GraphicsInterface & gi = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface());

            if( useSuperGrid && adjustPlotsForSuperGrid )
            {
                printF("solveHelmholtzDirect: adjustPlotsForSuperGrid : WARNING: This will change the computed solution !!\n");

                cgWave.adjustSolutionForSuperGrid( u );
            }

            PlotStuffParameters psp;
            char buff[100];
            psp.set(GI_TOP_LABEL,sPrintF(buff,"Helmholtz (complex) omega=%g, d=%g, SG=%d",omega,damp,useSuperGrid));

            if( cg.numberOfDimensions()==2 )
                psp.set(GI_PLOT_CONTOUR_LINES,false);

      //  aString title=sPrintF("Super-grid functions: width=%g p=%g q=%g",superGridWidth,pSG,qSG);
            gi.erase();
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            PlotIt::contour( gi, u, psp );

        }
        



    }
    return 0;
}