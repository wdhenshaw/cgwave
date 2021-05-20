// ----Form and solve the Helmholtz equation using a direct or iterative solver ----

#include "CgWaveHoltz.h"
#include "CgWave.h"
// #include "Overture.h" 
// #include "MappedGridOperators.h"
#include "Oges.h"
#include "ParallelUtility.h"
#include "CompositeGridOperators.h"

#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )

// ============================================================================================
/// \brief Form and solve the Helmholtz equation using a direct or iterative solver.
// ============================================================================================
int CgWaveHoltz::solveHelmholtz( realCompositeGridFunction & u, realCompositeGridFunction & f  )
{
  int debug=1; 

  const real & omega     = dbase.get<real>("omega");

  CgWave & cgWave        = *dbase.get<CgWave*>("cgWave");
  const real & c         = cgWave.dbase.get<real>("c");

  const int & orderOfAccuracy = cgWave.dbase.get<int>("orderOfAccuracy");

  printF("\n ============== DIRECT SOLVE OF THE HELMHOLTZ EQUATION =============\n");
  printF("    c=%.4g, omega=%.6g, orderOfAccuracy=%d\n", c,omega,orderOfAccuracy);
  printF(" ===================================================================\n");

  CompositeGrid & cg = *u.getCompositeGrid();
  const int & numberOfComponentGrids = cg.numberOfComponentGrids(); 

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

  // Fill in the forcing and boundary conditions: 
  cgWave.getHelmholtzForcing( f );

  // Index Ib1,Ib2,Ib3; 
  // Index Ig1,Ig2,Ig3; 
  // for( int grid=0; grid<numberOfComponentGrids; grid++ )
  // {
  //   MappedGrid & mg = cg[grid];
  //   for( int axis=0; axis<cg.numberOfDimensions(); axis++ )
  //   {
  //     for( int side=0; side<=1; side++ )
  //     {
  //       if( mg.boundaryCondition(side,axis)==CgWave::dirichlet )
  //       {
  //          getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
  //          OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
  //          fLocal(Ib1,Ib2,Ib3)=0.; // Dirichlet BC 

  //       }
  //       else if( mg.boundaryCondition(side,axis)==CgWave::neumann )
  //       {
  //          getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
  //          OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
  //          fLocal(Ig1,Ig2,Ig3)=0.;              // Neumann BC -- fill ghost values
  //       }
  //       else
  //       {
  //         printF("CgWave::solveHelmholtz:ERROR: unknown boundaryCondition=%d for (side,axis,grid)=(%d,%d,%d)\n", mg.boundaryCondition(side,axis),side,axis,grid);
  //         OV_ABORT("ERROR");
  //       }   
  //     }
  //   }
  // }

  u=0.;  // initial guess for iterative solvers
  real time0=getCPU();

  // ------- SOLVE THE HELMHOLTZ EQUATIONS -----
  solver.solve( u,f );   


  real time= ParallelUtility::getMaxValue(getCPU()-time0);
  printF("\n*** max residual=%8.2e, time for 1st solve of the Dirichlet problem = %8.2e (iterations=%i) ***\n",
         solver.getMaximumResidual(),time,solver.getNumberOfIterations());


  return 0;
}