// This file automatically generated from solveHelmholtz.bC with bpp.
// ----Form and solve the Helmholtz equation using a direct or iterative solver ----

#include "CgWaveHoltz.h"
#include "CgWave.h"
// #include "Overture.h" 
// #include "MappedGridOperators.h"
#include "Oges.h"
#include "ParallelUtility.h"
#include "CompositeGridOperators.h"



#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )

// ============================================================================================
/// \brief Form and solve the Helmholtz equation using a direct or iterative solver.
// ============================================================================================
int CgWaveHoltz::solveHelmholtz( realCompositeGridFunction & u, realCompositeGridFunction & f  )
{
    int debug=1; 

    const real & omega     = dbase.get<real>("omega");

    CgWave & cgWave        = *dbase.get<CgWave*>("cgWave");
    const real & c         = cgWave.dbase.get<real>("c");

    const int & orderOfAccuracy       = cgWave.dbase.get<int>("orderOfAccuracy");
    const int & numberOfFrequencies   = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray  = dbase.get<RealArray>("frequencyArray");

    printF("\n ============== DIRECT SOLVE OF THE HELMHOLTZ EQUATION =============\n");
    printF("    c=%.4g, numberOfFrequencies=%d, orderOfAccuracy=%d\n", c,numberOfFrequencies,orderOfAccuracy);
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
        printF(" freq(%d) = %.6g\n",freq,frequencyArray(freq));
    }
    printF(" ===================================================================\n");

    ::display(frequencyArray,"frequencyArray");

    if( omega != frequencyArray(0) )
    {
        printF("solveHelmholtz: ERROR omega=%g, != frequencyArray(0)=%g (TEMP SANITY CHECK -- FIX ME)\n",omega,frequencyArray(0));
        OV_ABORT("error"); 
    }

    CompositeGrid & cg = *u.getCompositeGrid();
    const int & numberOfComponentGrids = cg.numberOfComponentGrids(); 

  // Fill in the forcing and boundary conditions: 
    cgWave.getHelmholtzForcing( f );

  // --- DO THIS FOR NOW : **FIX ME** -----
    const int boxHelmholtz=3; 
    Real kxBoxHelmholtz=1., kyBoxHelmholtz=1., kzBoxHelmholtz=1.;
    int knownSolutionOption=0; // no known solution
    if( cgWave.dbase.has_key("userDefinedKnownSolutionData") )
    {
    // printF("++++ solveHelmholtz: userDefinedKnownSolutionData is found\n");
        DataBase & db =  cgWave.dbase.get<DataBase>("userDefinedKnownSolutionData");
        const aString & userKnownSolution = db.get<aString>("userKnownSolution");
        if( userKnownSolution=="planeWave"  )
        {
            knownSolutionOption=1;                   // this number must match in bcOptWave.bf90
        }
        else if( userKnownSolution=="gaussianPlaneWave"  ) 
        {
            knownSolutionOption=2;                   // this number must match in bcOptWave.bf90
        }    
        else if( userKnownSolution=="boxHelmholtz"  ) 
        {
            knownSolutionOption=boxHelmholtz;                   // this number must match in bcOptWave.bf90
            kxBoxHelmholtz = cgWave.dbase.get<Real>("kxBoxHelmholtz");
            kyBoxHelmholtz = cgWave.dbase.get<Real>("kyBoxHelmholtz");
            kzBoxHelmholtz = cgWave.dbase.get<Real>("kzBoxHelmholtz"); 
            printF("\n @@@@ solveHelmholtz: boxHelmholtz: kxBoxHelmholtz=%g\n\n",kxBoxHelmholtz);     
        }
        else if( userKnownSolution=="polyPeriodic"  ) 
        {
            knownSolutionOption=4;                   // this number must match in bcOptWave.bf90
        } 

    } 
    else
    {
     // printF("+++ solveHelmholtz: userDefinedKnownSolutionData is NOT found\n");
    } 

    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
        const Real omega = frequencyArray(freq);

        Oges solver( cg );                     // create a solver

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


        IntegerArray boundaryConditions(2,3,numberOfComponentGrids);
        boundaryConditions = OgesParameters::dirichlet; // default
        RealArray bcData(2,2,3,numberOfComponentGrids);
        Range all; 

    // Solve constCoeff(0,grid)*I +constCoeff(1,grid)*Laplacian 
    // We solve:  omega^2 I + c^2 Delta = f 
        RealArray constantCoeff(2,numberOfComponentGrids);
    // constantCoeff(0,all) = -SQR(omega);    // FIX SIGN IN CGWAVE 
    // constantCoeff(1,all) = -SQR(c); 

        constantCoeff(0,all) = SQR(omega);    
        constantCoeff(1,all) = SQR(c); 

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
                  else if( mg.boundaryCondition(side,axis) > 0 )
                  {
                      printF("CgWave::solveHelmholtz:ERROR: unknown boundaryCondition=%d for (side,axis,grid)=(%d,%d,%d)\n", mg.boundaryCondition(side,axis),side,axis,grid);
                      OV_ABORT("ERROR");
                  }

            }

        }

    // boundaryConditions=OgesParameters::dirichlet;   // ** FIX ME ****
        bcData=0.;

        CompositeGridOperators op(cg);
        op.setOrderOfAccuracy(orderOfAccuracy);

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
                          OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
                          OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);

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
                          OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
             // fLocal(Ig1,Ig2,Ig3)=0.;              // Neumann BC -- fill ghost values
                    }
                    else if( mg.boundaryCondition(side,axis)>0 )
                    {
                        printF("CgWave::solveHelmholtz:ERROR: unknown boundaryCondition=%d for (side,axis,grid)=(%d,%d,%d)\n", mg.boundaryCondition(side,axis),side,axis,grid);
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
                OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
                OV_GET_SERIAL_ARRAY(real,fTemp[grid],fTempLocal);
                getIndex(cg[grid].dimension(),I1,I2,I3);
                fTempLocal(I1,I2,I3)=fLocal(I1,I2,I3,freq);
            }

            solver.solve( uTemp,fTemp ); 

            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                OV_GET_SERIAL_ARRAY(real,uTemp[grid],uTempLocal);
                getIndex(cg[grid].dimension(),I1,I2,I3);
                uLocal(I1,I2,I3,freq)=uTempLocal(I1,I2,I3);
            }
        }  


        real time= ParallelUtility::getMaxValue(getCPU()-time0);
        printF("\n*** freq=%d: omega=%.5g, max residual=%8.2e, time for direct Helmholtz solve = %8.2e (s) (iterations=%i) ***\n",
                      freq,omega,solver.getMaximumResidual(),time,solver.getNumberOfIterations());

    } // end for freq 
    return 0;
}