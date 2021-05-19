// This file automatically generated from tcmWideStencil.bC with bpp.
//===============================================================================
//
//  Coefficient Matrix Example for a WIDE STENCIL (e.g. fitting a wider upwind dissipation into a regular stencil )
//
// 
//
// Usage: `tcmWideStencil [-g=<gridName>] [-solver=[yale][harwell][slap][petsc][mg]] [-debug=<value>] [-outputMatrix] ...
//              [-noTiming] [-check] [-plot] [-trig] [-tol=<value>] [-freq=<value>] [-dirichlet] [-neumann]' 
//
//   The -check option is used for regression testing -- it will test various solvers on a few grids
//
// NOTE:
// To get PETSc log info, compile PETScEquationSolver with the destructor calling PetscFinalize()
//   and use the command line arg -log_summary
//   memory usage: -trmalloc_log 
//
// Parallel examples:
//   mpirun -np 2 tcmWideStencil square20.hdf -solver=petsc
//   mpirun -np 2 tcmWideStencil cic.hdf -solver=petsc 
//   srun -N1 -n1 -ppdebug tcmWideStencil square20.hdf -solver=petsc 
//   srun -N1 -n2 -ppdebug tcmWideStencil sibe2.order2.hdf -solver=petsc 
//==============================================================================
#include "Overture.h"  
#include "MappedGridOperators.h"
#include "Oges.h"
#include "CompositeGridOperators.h"
#include "SquareMapping.h"
#include "AnnulusMapping.h"
#include "OGTrigFunction.h"
#include "OGPolyFunction.h"
#include "SparseRep.h"
#include "display.h"
#include "Ogmg.h"
#include "Checker.h"
#include "PlotStuff.h"
#include "ParallelUtility.h"
#include "LoadBalancer.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )



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
///   This follows the function in CgWave  
// ============================================================================================
int formImplicitTimeSteppingMatrix( CompositeGrid & cg, realCompositeGridFunction & impCoeff, Oges & impSolver, CompositeGridOperators & op   )
{
    real cpu0=getCPU();


 // make some shorter names for readability
    BCTypes::BCNames
        dirichlet           = BCTypes::dirichlet,
        neumann             = BCTypes::neumann,
        mixed               = BCTypes::mixed,
        extrapolate         = BCTypes::extrapolate,
        allBoundaries       = BCTypes::allBoundaries; 

    DataBase dbase;

    if( true)
    {
        dbase.put<int>("debug") =3 ;
        dbase.put<real>("dt") =.01;
        dbase.put<real>("c") = 1.;
        dbase.put<int>("orderOfAccuracy") = 2;
        dbase.put<int>("orderOfAccuracyInTime") = 2 ;
        IntegerArray & gridIsImplicit = dbase.put<IntegerArray>("gridIsImplicit");
        gridIsImplicit.redim(cg.numberOfComponentGrids());
        gridIsImplicit=1;
        dbase.put<real>("ad4"); // coeff of the artificial dissipation.
    
        RealArray & cImp= dbase.put<RealArray>("cImp");  
        cImp.redim(Range(-1,1));
        cImp(-1)=.25; 
        cImp( 0)=.5; 
        cImp(+1)=.25; 
    }


    int & debug                          = dbase.get<int>("debug");
    real & dt                            = dbase.get<real>("dt");
    const real & c                       = dbase.get<real>("c");
    const int & orderOfAccuracy          = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime    = dbase.get<int>("orderOfAccuracyInTime");
    const IntegerArray & gridIsImplicit  = dbase.get<IntegerArray>("gridIsImplicit");
    const real & ad4                     = dbase.get<real>("ad4"); // coeff of the artificial dissipation.

    bool addUpwinding = ad4>0.;


    addUpwinding=true; // 



    printF("\n ==================== FORM MATRIX FOR IMPLICI TIME-STEPPING ===================\n");
    printF("   c=%.4g, dt=%9.3e, orderOfAccuracy=%d, orderOfAccuracyInTime=%d addUpwinding=%d\n", c,dt,orderOfAccuracy,orderOfAccuracyInTime,addUpwinding);
    printF(" ================================================================================\n");

    const int & numberOfComponentGrids = cg.numberOfComponentGrids(); 
    const int & numberOfDimensions = cg.numberOfDimensions(); 

  // coefficients in implicit time-stepping  
  //  D+t D-t u = c^2 Delta( cImp(1) *u^{n+1} + cImp(0) *u^n + cImp(-1)* u^{n-1} )
    RealArray & cImp              = dbase.get<RealArray>("cImp");  

  // if( !dbase.has_key("impSolver") )
  // {
  //   dbase.put<Oges>("impSolver");
  // }
  // Oges & impSolver = dbase.get<Oges>("impSolver");
    impSolver.updateToMatchGrid( cg );                     

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


  // CompositeGridOperators & op = dbase.get<CompositeGridOperators>("operators");
    op.setOrderOfAccuracy(orderOfAccuracy);

    bool usePredefined= false; // true = old way

    if( usePredefined )
    {
    // // ---- use Oges predefined equations ***OLD WAY*** ----
    
    // IntegerArray boundaryConditions(2,3,numberOfComponentGrids);
    // RealArray bcData(2,2,3,numberOfComponentGrids);
    // bcData=0.;

    // Range all; 

    // // Solve constCoeff(0,grid)*I +constCoeff(1,grid)*Laplacian 
    // // We solve:  I - alpha*(c^2*dt^2)* Delta = ...
    // Real alpha=cImp(-1);
    // RealArray constantCoeff(2,numberOfComponentGrids);

    // constantCoeff(0,all) = 1.;
    // constantCoeff(1,all) = - alpha*SQR(c*dt);

    // // Assign boundary conditions for Oges
    // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    // {
    //   MappedGrid & mg = cg[grid];
    //   ForBoundary(side,axis)
    //   {
    //      if( mg.boundaryCondition(side,axis)==CgWave::dirichlet )
    //      {
    //        boundaryConditions(side,axis,grid) = OgesParameters::dirichlet;
    //      }
    //      else if( mg.boundaryCondition(side,axis)==CgWave::neumann )
    //      { 
    //        boundaryConditions(side,axis,grid) = OgesParameters::neumann;
    //      }
    //      else if( mg.boundaryCondition(side,axis) <= 0 )
    //      { 
    //        boundaryConditions(side,axis,grid) = mg.boundaryCondition(side,axis);
    //      }       
    //      else if( mg.boundaryCondition(side,axis) > 0 )
    //      {
    //        printF("CgWave::formImplicitTimeSteppingMatrix:ERROR: unknown boundaryCondition=%d for (side,axis,grid)=(%d,%d,%d)\n", mg.boundaryCondition(side,axis),side,axis,grid);
    //        OV_ABORT("ERROR");
    //      }

    //   }

    // }



    // impSolver.setEquationAndBoundaryConditions( OgesParameters::heatEquationOperator ,op,boundaryConditions,bcData,constantCoeff );

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
        Real adSosup;
        const int upwindHalfStencilWidth = orderOfAccuracy; 
        if( orderOfAccuracy==2 )
        {
            adSosup=-c*dt*1./8.;
            for( int m=0; m<4; m++ )
                upwindCoeff[m] =upwindCoeff4[m];
        }
        else if( orderOfAccuracy==4 )
        {
            adSosup=c*dt*5./288.;
            for( int m=0; m<4; m++ )
                upwindCoeff[m] =upwindCoeff6[m];      
        }
        else if( orderOfAccuracy==6 )
        {
            adSosup=-c*dt*31./8640.;
      // upwindDissCoeff=upwindDissCoeff8;
        }
        else
        {
          OV_ABORT("ERROR orderOfAccuracy");
        }



        Range all;
        int stencilWidth = orderOfAccuracy + 1;
        int numberOfGhostLines= orderOfAccuracy/2;  // fix me for UPWIND

        int extraEntries = 1;  // we add 1 extra entry for interpolation equations

    // extraEntries = 3;  // ******** TEMP *******************

        if( addUpwinding )
        {
      // -- Note: we do not always have to add extra entries for upwinding
      //          e.g. on Cartesian grids there are zeros in the existing stencil that could be used. 
            extraEntries = 2*numberOfDimensions; // for upwinding equations
      // stencilWidth += 2;
      // numberOfGhostLines +=1; 
        }

        const int baseStencilSize = pow(stencilWidth,cg.numberOfDimensions());   // number of entries in default stencil 
        const int stencilSize=int( baseStencilSize + extraEntries );                      // add extra for interpolation and upwind equations

        const int numberOfComponentsForCoefficients=1;
        const int stencilDimension=stencilSize*SQR(numberOfComponentsForCoefficients);

        const int baseStencilDimension=baseStencilSize*SQR(numberOfComponentsForCoefficients);



        printF(">>>> stenciWidth=%d, stencilSize=%d, numberOfGhostLines=%d, addUpwinding=%d\n",stencilWidth,stencilSize,numberOfGhostLines,addUpwinding);

    // use this coeff matrix: **FIX ME**
        if( !dbase.has_key("impCoeff") )
        {
            dbase.put<realCompositeGridFunction>("impCoeff");
        }
    // realCompositeGridFunction & impCoeff = dbase.get<realCompositeGridFunction>("impCoeff"); 

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

        int cc = 0; // component number 

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
      // We solve:  I - cImp(-1) * (c^2*dt^2)* Delta = ...

            const int mDiag = M123(0,0,0);              // index of diagonal entry

            if( gridIsImplicit(grid)==1 )
            {
        // ----- this grid is adavnced with IMPLICIT time-stepping ----
                getIndex(mg.gridIndexRange(),I1,I2,I3);
                RealArray lapCoeff(M0,I1,I2,I3);
                mgop.assignCoefficients(MappedGridOperators::laplacianOperator,lapCoeff,I1,I2,I3,0,0); // 
                

                if( true )
                {
                    coeffLocal(M0,I1,I2,I3)  = lapCoeff;  // Just put in the Laplacian for now 
                }
                else 
                {
                    Real ccLap = - cImp(-1)*SQR(c*dt); 
                    coeffLocal(M0,I1,I2,I3)  = ccLap*lapCoeff;

          // set diagonal entry
                    coeffLocal(mDiag,I1,I2,I3) += 1.0;
                }
            }
            else
            {
        // ----- this grid is adavnced with EXPLICIT time-stepping ----
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

                assert( isRectangular );

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
                 // ---Upwind coefficients for a curvlinear grid ---

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
                                const Real upwStencilValue = upwindCoeff[upwCase][iStencil+upwindHalfStencilWidth]; 
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
                                        assert( extraStencilLocation<stencilDimension );

                                        m = extraStencilLocation;
                                        extraStencilLocation++;

                                    }

                                    int ii  = indexToEquation( cc,i1,i2,i3 ); 
                                    int iis = indexToEquation( cc,i1s,i2s,i3s ); 
                                    printF("UPWIND: (i1,i2)=(%3d,%3d) : eqn=%d,  dir=%d ->  ADD entry m=%3d, (i1s,i2s)=(%3d,%3d) : eqn=%d, adxSosup=%9.3e, upwStencilValue = %9.3e\n",
                                                      i1,i2,ii, dir,m,i1s,i2s,iis, adxSosup[dir],upwStencilValue);

                  // ***TEST*** Just add difference coefficients 
                                    if( testUPW )
                                        coeffLocal(m,i1,i2,i3) += upwStencilValue;  // ************ TEMP *************
                                    else
                                        coeffLocal(m,i1,i2,i3) += adxSosup[dir] * upwStencilValue;


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

        // if( mg.boundaryCondition(side,axis)
                if( mg.boundaryCondition(side,axis) > 0 ) // ***************************** FORCE DIRICHLET FOR NOW
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

            if( addUpwinding )
            {
                ::display(impCoeff[grid],sPrintF("implicit time-stepping matrix on grid=%d",grid));
                displayCoeff(impCoeff[grid],sPrintF("implicit time-stepping matrix on grid=%d",grid));
        // OV_ABORT("stop here for now");
            }

        }

        impCoeff.finishBoundaryConditions(); 

        impSolver.setCoefficientArray( impCoeff );   // supply coefficients to Oges

    }

  // timing(timeForInitialize) += getCPU()-cpu0;

  
    return 0;
}







bool measureCPU=TRUE;
real
CPU()
// In this version of getCPU we can turn off the timing
{
    if( measureCPU )
        return getCPU();
    else
        return 0;
}

void
plotResults( PlotStuff & ps, Oges & solver, realCompositeGridFunction & u, realCompositeGridFunction & err )
// ==============================================================================================
// Plot results from Oges
// ==============================================================================================
{
            
    GraphicsParameters psp;

    aString answer;
    aString menu[]=
    {
        "solution",
        "error",
        "grid",
        "erase",
        "exit",
        ""
    };
        
    for( ;; )
    {
        ps.getMenuItem(menu,answer,"choose an option");
        if( answer=="exit" )
        {
            break;
        }
        else if( answer=="solution" )
        {
            psp.set(GI_TOP_LABEL,"Solution u"); 
            PlotIt::contour(ps,u,psp);
        }
        else if( answer=="error" )
        {
            psp.set(GI_TOP_LABEL,"error"); 
            PlotIt::contour(ps,err,psp);
        }
        else if( answer=="grid" )
        {
            psp.set(GI_TOP_LABEL,"grid"); 
            PlotIt::plot(ps,*u.getCompositeGrid(),psp);
        }
        else if( answer=="erase" )
        {
            ps.erase();
        }
            
    }

}


int 
assignForcing(int option, CompositeGrid & cg, realCompositeGridFunction & f, OGFunction & exact,
                            RealArray *varCoeff=NULL )
// ================================================================================================
//
/// \brief  Assign the right-hand-side. 
///
/// \param option (input) : 0 = dirichlet, 1=neumann, 2=mixed (variable coefficients)
// 
// ================================================================================================
{
    const int numberOfDimensions = cg.numberOfDimensions();
    
    Index I1,I2,I3;
    Index Ib1,Ib2,Ib3;
    Index Ig1,Ig2,Ig3;


    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
    // mg.mapping().getMapping().getGrid();
    // printF(" signForJacobian=%e\n",mg.mapping().getMapping().getSignForJacobian());
#ifdef USE_PPP
        realSerialArray fLocal; getLocalArrayWithGhostBoundaries(f[grid],fLocal);
#else
        realSerialArray & fLocal = f[grid]; 
#endif

        getIndex(mg.indexRange(),I1,I2,I3);  
        int includeGhost=1; // include parallel ghost pts in fLocal:
        bool ok = ParallelUtility::getLocalArrayBounds(f[grid],fLocal,I1,I2,I3,includeGhost);
        if( !ok ) continue; // there are no points on this processor.

        realArray & x= mg.center();
#ifdef USE_PPP
        realSerialArray xLocal; getLocalArrayWithGhostBoundaries(x,xLocal);
#else
        const realSerialArray & xLocal = x;
#endif

    // Assign the forcing : f = e.xx + e.yy + e.zz (e=exact solution)
        RealArray ed(I1,I2,I3);
        const int rectangularForTZ=0;
        fLocal=0.;
        for( int axis=0; axis<cg.numberOfDimensions(); axis++ )
        {
            int ntd=0, nxd[3]={0,0,0}; //
            nxd[axis]=2;  // compute e.xx (axis=0), e.yy (axis=1), ...
            exact.gd( ed,xLocal,mg.numberOfDimensions(),rectangularForTZ,ntd,nxd[0],nxd[1],nxd[2],I1,I2,I3,0,0.);
            fLocal(I1,I2,I3)+=ed(I1,I2,I3);
        }
                

        ForBoundary(side,axis)
        {
            if( mg.boundaryCondition(side,axis) > 0 )
            {
                #ifdef USE_PPP
                    const realSerialArray & normal  = mg.vertexBoundaryNormalArray(side,axis);
                #else
                    const realSerialArray & normal  = mg.vertexBoundaryNormal(side,axis);
                #endif

                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);

                bool ok = ParallelUtility::getLocalArrayBounds(f[grid],fLocal,Ib1,Ib2,Ib3,includeGhost);
                if( !ok ) continue; // there are no points on this processor.
                if( option==0 )
                {
          // Dirichlet BC's : assign the value of f on the boundary: 
                    RealArray ue(Ib1,Ib2,Ib3);
                    exact.gd( ue,xLocal,mg.numberOfDimensions(),rectangularForTZ,0,0,0,0,Ib1,Ib2,Ib3,0,0.);
                    fLocal(Ib1,Ib2,Ib3)=ue;
                }
                else
                {
          // Neumann or mixed BC's : assign the value of f on the ghost line: 

                    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
                    bool ok = ParallelUtility::getLocalArrayBounds(f[grid],fLocal,Ig1,Ig2,Ig3,includeGhost);
                    if( !ok ) continue; // there are no points on this processor.

                    realSerialArray uex(Ib1,Ib2,Ib3), uey(Ib1,Ib2,Ib3);
                    exact.gd( uex,xLocal,numberOfDimensions,rectangularForTZ,0,1,0,0,Ib1,Ib2,Ib3,0,0.);
                    exact.gd( uey,xLocal,numberOfDimensions,rectangularForTZ,0,0,1,0,Ib1,Ib2,Ib3,0,0.);

                    fLocal(Ig1,Ig2,Ig3) = normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey;
                    if( numberOfDimensions==3 )
                    {
                        exact.gd( uex,xLocal,numberOfDimensions,rectangularForTZ,0,0,0,1,Ib1,Ib2,Ib3,0,0.);  // uex = T.z
                        fLocal(Ig1,Ig2,Ig3) +=normal(Ib1,Ib2,Ib3,2)*uex;
                    }
                    if( option==2 )
                    { // -- Mixed BC's ---
                        RealArray ue(Ib1,Ib2,Ib3);
                        exact.gd( ue,xLocal,mg.numberOfDimensions(),rectangularForTZ,0,0,0,0,Ib1,Ib2,Ib3,0,0.);
                        fLocal(Ig1,Ig2,Ig3) = (  varCoeff[grid](Ib1,Ib2,Ib3,0)*ue(Ib1,Ib2,Ib3)
                                                                        +varCoeff[grid](Ib1,Ib2,Ib3,1)*fLocal(Ig1,Ig2,Ig3) );
                    }
                    
                }
                
            }
        }
    }
  // f.applyBoundaryCondition(0,BCTypes::dirichlet,BCTypes::allBoundaries,0.);   
  // f.display("Here is f");

    return 0;
}

int 
computeTheError( int option, CompositeGrid & cg, realCompositeGridFunction & u,
                                  realCompositeGridFunction & err, OGFunction & exact, real & error ) 
// ================================================================================================
//
//  Compute the error in the solution.
//
// /option (input) : 0 = dirichlet, 1=neumann, 2=mixed
// 
// ================================================================================================
{

    err=0.;
    error=0.;
    real errorWithGhostPoints=0;
    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        realArray & x= mg.center();
#ifdef USE_PPP
        realSerialArray xLocal; getLocalArrayWithGhostBoundaries(x,xLocal);
        realSerialArray uLocal; getLocalArrayWithGhostBoundaries(u[grid],uLocal);
        realSerialArray errLocal; getLocalArrayWithGhostBoundaries(err[grid],errLocal);
        intSerialArray maskLocal; getLocalArrayWithGhostBoundaries(cg[grid].mask(),maskLocal);
#else
        const realSerialArray & xLocal = x;
        realSerialArray & uLocal = u[grid]; 
        realSerialArray & errLocal = err[grid]; 
        const intSerialArray & maskLocal = cg[grid].mask();
#endif

        getIndex(cg[grid].indexRange(),I1,I2,I3,1);  
        int includeGhost=1; // include parallel ghost pts in uLocal
        bool ok = ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);

        real ueMax=0.; // holds the max value of the exact soln on this grid
        RealArray ue;
        if( ok )
        { // evaluate the exact solution
            ue.redim(I1,I2,I3);
            const int rectangularForTZ=0;
            exact.gd( ue,xLocal,mg.numberOfDimensions(),rectangularForTZ,0,0,0,0,I1,I2,I3,0,0.);
            ueMax=max(fabs(ue));
        }
        ueMax=ParallelUtility::getMaxValue(ueMax); // max value over all procs

        real gridErrWithGhost=0., gridErr=0.;
        if( ok )
        {
            where( maskLocal(I1,I2,I3)!=0 )
                errLocal(I1,I2,I3)=abs(uLocal(I1,I2,I3)-ue);

            gridErrWithGhost=max(errLocal(I1,I2,I3))/ueMax;

            getIndex(cg[grid].indexRange(),I1,I2,I3);  
            bool ok = ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
            if( !ok ) continue; // there are no points on this processor.

            where( maskLocal(I1,I2,I3)!=0 )
                errLocal(I1,I2,I3)=abs(uLocal(I1,I2,I3)-ue(I1,I2,I3));

            gridErr=max(errLocal(I1,I2,I3))/ueMax;
        }
        gridErr         =ParallelUtility::getMaxValue(gridErr); // max value over all procs
        gridErrWithGhost=ParallelUtility::getMaxValue(gridErrWithGhost); // max value over all procs

        error=max(error, gridErr );
        errorWithGhostPoints=max(errorWithGhostPoints, gridErrWithGhost);

        printF(" grid=%i (%s) max. rel. err=%e (%e with ghost)\n",grid,(const char*)cg[grid].getName(),
                      gridErr,gridErrWithGhost);

        if( Oges::debug & 4 )
        {
            display(u[grid],"solution u");
            display(err[grid],"abs(error on indexRange +1)");
      // abs(u[grid](I1,I2,I3)-exact(cg[grid],I1,I2,I3,0)).display("abs(error)");
        }
    }
    if( option==0 )
        printF("Maximum relative error with dirichlet bc's= %e (%e with ghost)\n",error,errorWithGhostPoints);  
    else if(option==1 )
        printF("Maximum relative error with neumann bc's= %e\n",error);  
    else if(option==1 )
        printF("Maximum relative error with mixed bc's= %e\n",error);  

    return 0;
}


int 
main(int argc, char *argv[])
{
    Overture::start(argc,argv);  // initialize Overture

   // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  // *** INIT_PETSC_SOLVER(); // ** TURN OFF FOR NOW -- May 17, 2021

    printF("Usage: tcmWideStencil -g=gridName [-solver=[yale][harwell][slap][petsc][mg]] [-debug=<value>][-outputMatrix]\n" 
                                          "[-noTiming] [-check] [-trig] [-tol=<value>] [-order=<value>] [-plot] [-ilu=] [-gmres] \n"
                                          "[-freq=<value>] [-dirichlet] [-neumann] [-mixed] [-testCommunicator] [-hypre] [-predefined]\n");

    const int maxNumberOfGridsToTest=3;
    int numberOfGridsToTest=maxNumberOfGridsToTest;
    aString gridName[maxNumberOfGridsToTest] =   { "square5", "cic", "sib" };
  // here are upper bounds on the errors we expect for each grid. This seems the only reliable
  // way to compare results from different machines, especially for iterative solvers.
    const real errorBound[maxNumberOfGridsToTest][2][2]=
        { 5.e-8,4.e-8,    5.e-7,9.e-7,  // square, dirichlet/neuman(DP) dir/neu(SP)
            7.e-4,2.e-3,    7.e-4,2.e-3, // cic
            6.e-3,7.e-3,    6.e-3,7.e-3  // sib
        };
    const int precision = REAL_EPSILON==DBL_EPSILON ? 0 : 1;
    int twilightZoneOption=0;
    
    int solverType=OgesParameters::yale; 
    aString solverName="yale";
    aString iterativeSolverType="bi-cg";
    bool check=false;
    real tol=1.e-8;
    int orderOfAccuracy=2;
    int plot=0;
    int iluLevels=-1; // -1 : use default
    int problemsToSolve=1+2;  // solve dirichlet=1 and neumann=2
    bool outputMatrix=false;
    bool testCommunicator=false;  // set to true to test PETSc when using only a subset of the processors.
    bool usePredefined=false;
    
    real fx=2., fy=2., fz=2.; // frequencies for trig TZ
    
    int len=0;
    if( argc >= 1 )
    { 
        for( int i=1; i<argc; i++ )
        {
            aString arg = argv[i];
            if( arg=="-noTiming" )
                measureCPU=FALSE;
            else if( (len=arg.matches("-debug=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&Oges::debug);
                printF("Setting Oges::debug=%i\n",Oges::debug);
            }
            else if( (len=arg.matches("-tol=")) )
            {
                sScanF(arg(len,arg.length()-1),"%e",&tol);
                printF("Setting tol=%e\n",tol);
            }
            else if( (len=arg.matches("-freq=")) )
            {
                sScanF(arg(len,arg.length()-1),"%e",&fx);
                fy=fx; fz=fx;
                printF("Setting fx=fy=fz=%e\n",fx);
            }
            else if( (len=arg.matches("-ilu=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&iluLevels);
                printF("Setting ilu levels =%i\n",iluLevels);
            }
            else if( (len=arg.matches("-gmres")) )
            {
                iterativeSolverType="gmres";
            }
            else if( (len=arg.matches("-hypre")) )
            {
                iterativeSolverType="hypre";
            }
            else if( (len=arg.matches("-testCommunicator")) )
            {
                testCommunicator=true;
                printF("Test the parallel PETSc solver using a subset of the processors\n");
            }
            else if( (len=arg.matches("-outputMatrix")) )
            {
                outputMatrix=true;
            }
            else if( (len=arg.matches("-dirichlet")) )
            {
                problemsToSolve=1; // just solve dirichlet problem
            }
            else if( (len=arg.matches("-neumann")) )
            {
                problemsToSolve=2; // just solve neumann problem
            }
            else if( (len=arg.matches("-mixed")) )
            {
                problemsToSolve=4; // just solve with mixed BC's (variable coeff)
            }
            else if( (len=arg.matches("-order=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&orderOfAccuracy);
                if( orderOfAccuracy!=2 && orderOfAccuracy!=4 && orderOfAccuracy!=6 )
                {
                    printF("ERROR: orderOfAccuracy should be 2, 4 or 6!\n");
                    OV_ABORT("ERROR");
                }
                printF("Setting orderOfAccuracy=%i\n",orderOfAccuracy);
            }
            else if( arg(0,7)=="-solver=" )
            {
                solverName=arg(8,arg.length()-1);
                if( solverName=="yale" )
                    solverType=OgesParameters::yale;
                else if( solverName=="harwell" )
                    solverType=OgesParameters::harwell;
                else if( solverName=="petsc" || solverName=="PETSc" )
#ifdef USE_PPP 
                {
                    solverType=OgesParameters::PETScNew;
          // printF("tcmWideStencil: Setting solverType=PETScNew = %i\n",(int)solverType);
                }
#else
                solverType=OgesParameters::PETSc;
#endif
                else if( solverName=="slap" || solverName=="SLAP" )
                    solverType=OgesParameters::SLAP;
                else if( solverName=="mg" || solverName=="multigrid" )
                    solverType=OgesParameters::multigrid;
                else
                {
                    printF("Unknown solver=%s \n",(const char*)solverName);
                    throw "error";
                }
                
                printF("Setting solverType=%i\n",solverType);
            }
            else if( arg=="-plot" )
            {
                plot=true;
            }
            else if( arg=="-check" )
            {
                check=true;
            }
            else if( arg=="-trig" )
            {
                twilightZoneOption=1;
            }
            else if( (len=arg.matches("-predefined")) )
            {
                usePredefined=true; // use predefined equations
                printF("Setting usePredfined=true: define problems from Oges predefined equations\n");
            }
            else if( (len=arg.matches("-g=")) ) // *new way*
            {
                numberOfGridsToTest=1;
                gridName[0]=arg(len,arg.length()-1);
                printF("Use grid=[%s]\n",(const char*)gridName[0]);
            }
            else
            {
                numberOfGridsToTest=1;
                gridName[0]=argv[1];
            }
        }
    }

    if( Oges::debug > 3 )
        SparseRepForMGF::debug=3;  

    aString checkFileName;
    if( REAL_EPSILON == DBL_EPSILON )
        checkFileName="tcmWideStencil.dp.check.new";  // double precision
    else  
        checkFileName="tcmWideStencil.sp.check.new";
    Checker checker(checkFileName);  // for saving a check file.

    printF("=================================================================================\n"
                  " --- tcmWideStencil --- test coefficient matrices: scalar problem on an overlapping grid   \n"
                  " \n"
                  "  Equation: Poisson.\n"
                  "  order of accuracy=%i.\n"
                  "  usePredefined=%i.\n"
                  ,orderOfAccuracy,(int)usePredefined);
    if( twilightZoneOption==0 )
        printF(" TwilightZone: polynomial, degree=%i.\n",orderOfAccuracy);
    else
        printF(" TwilightZone: trigonometric, fx=fy=fz=%e.\n",fx);
        
    printF("=================================================================================\n");

    PlotStuff ps(false,"tcmWideStencil");

  // make some shorter names for readability
    BCTypes::BCNames
        dirichlet           = BCTypes::dirichlet,
        neumann             = BCTypes::neumann,
        mixed               = BCTypes::mixed,
        extrapolate         = BCTypes::extrapolate,
        allBoundaries       = BCTypes::allBoundaries; 

    bool communicatorWasNewed=false;
    int pStart=-1, pEnd=-1;  // Use to distribute the grids in the CompositeGrid when we define a communicator
#ifdef USE_PPP
    MPI_Comm myComm=MPI_COMM_WORLD, myCommWorld=MPI_COMM_WORLD;  // Have PETSc use these communicators
    if( testCommunicator && solverType==OgesParameters::PETScNew )
    {
    // Here we test the PETSc solver when only some processors are involved in the parallel solve.
    // We build an MPI communicator that only includes some processors.
        communicatorWasNewed=true;

        const int np=max(1,Communication_Manager::Number_Of_Processors);
        const int myid=max(0,Communication_Manager::My_Process_Number);

        MPI_Group worldGroup, myGroup;
        MPI_Comm_group( MPI_COMM_WORLD,&worldGroup ); //get world group

        const int numRanks=min(2,np);
        int ranks[4]={1,2,3,4};  // include these ranks in myGroup
        printF("--TCM3-- Active processors = [%i",ranks[0]);
        for( int r=1; r<numRanks; r++ ){ printF(",%i",ranks[r]); }  // 
        printF("]\n");

        MPI_Group_incl(worldGroup, numRanks, ranks, &myGroup );   // myGgroup includes some ranks 
        MPI_Comm_create( MPI_COMM_WORLD, myGroup, &myComm ); // construct myComm
        MPI_Comm_create( MPI_COMM_WORLD, myGroup, &myCommWorld ); // construct myCommWorld
        MPI_Group_free(&worldGroup);
        MPI_Group_free(&myGroup);

        pStart=ranks[0]; pEnd=ranks[numRanks-1];

    // int colour = myid==(np-1);  // communicator includes this processor
    // // int colour = myid==(np-1);  // communicator includes this processor
    // // int colour = 1; 
    // int key=0;
    // MPI_Comm_split(MPI_COMM_WORLD,colour,key,&myComm );

    // if( myComm!=MPI_COMM_NULL )
    // {
    //   int size=-1; 
    //   MPI_Comm_size(myComm,&size);
    //   printf(" --TCM3-- myComm: myid=%i size=%i.\n",myid,size);
    // }
    // else
    // {
    //   printf(" --TCM3-- myComm: myid=%i myComm==NULL.\n",myid);
    // }
    }
    else
    {

    }
#endif

    int numberOfSolvers = check ? 2 : 1;
    real worstError=0.;
    for( int sparseSolver=0; sparseSolver<numberOfSolvers; sparseSolver++ )
    {
        if( check )
        {
            if( sparseSolver==0 )
            {
                solverName="yale";
                solverType=OgesParameters::yale;
            }
            else
            {
                solverName="slap";
                solverType=OgesParameters::SLAP;
            }
        }

        checker.setLabel(solverName,0);


        for( int it=0; it<numberOfGridsToTest; it++ )
        {
            aString nameOfOGFile=gridName[it];
            checker.setLabel(nameOfOGFile,1);

            printF("\n *****************************************************************\n"
                            " ******** Checking grid: %s                 ************ \n"
                            " *****************************************************************\n\n",(const char*)nameOfOGFile);

            CompositeGrid cg;
            if( pStart>=0  )
            {
                LoadBalancer loadBalancer;
                loadBalancer.setLoadBalancer(LoadBalancer::allToAll);
                printF("loadBalancer: pStart=%i, pEnd=%i\n",pStart,pEnd);
                
                loadBalancer.setProcessors(pStart, pEnd);

        // loadBalancer.setLoadBalancer(LoadBalancer::allToAll);
                getFromADataBase(cg,nameOfOGFile,loadBalancer);
            }
            else
            {
                getFromADataBase(cg,nameOfOGFile);
            }
            cg.displayDistribution("cg after reading.");
            
            cg.update(MappedGrid::THEmask | MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEvertexBoundaryNormal);

            if( Oges::debug >3 )
            {
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    displayMask(cg[grid].mask(),"mask");
            }


            const int inflow=1, outflow=2, wall=3;
    
      // create a twilight-zone function for checking the errors
            OGFunction *exactPointer;
            if( twilightZoneOption==1 ||
                    min(abs(cg[0].isPeriodic()(Range(0,cg.numberOfDimensions()-1))-Mapping::derivativePeriodic))==0 )
            {
        // this grid is probably periodic in space, use a trig function
                printF("TwilightZone: trigonometric polynomial, fx=%9.3e, fy=%9.3e, fz=%9.3e\n",fx,fy,fz);
                exactPointer = new OGTrigFunction(fx,fy,fz); 
            }
            else
            {
                printF("TwilightZone: algebraic polynomial\n");
        // cg.changeInterpolationWidth(2);

                int degreeOfSpacePolynomial = orderOfAccuracy; 
                int degreeOfTimePolynomial = 1;
                int numberOfComponents = cg.numberOfDimensions();
                exactPointer = new OGPolyFunction(degreeOfSpacePolynomial,cg.numberOfDimensions(),numberOfComponents,
                                                                                    degreeOfTimePolynomial);
        
            
            }
            OGFunction & exact = *exactPointer;

      // make a grid function to hold the coefficients
            Range all;
            Index I1,I2,I3, Ia1,Ia2,Ia3;

            const int width=orderOfAccuracy+1;
            int stencilSize=int(pow(width,cg.numberOfDimensions())+1);  // add 1 for interpolation equations

            realCompositeGridFunction coeff(cg,stencilSize,all,all,all); 

            const int numberOfGhostLines=orderOfAccuracy/2;
            coeff.setIsACoefficientMatrix(true,stencilSize,numberOfGhostLines);  
            coeff=0.;
        
      // create grid functions: 
            realCompositeGridFunction u(cg),f(cg);
            realCompositeGridFunction err(cg);

            real error;

            CompositeGridOperators op(cg);                            // create some differential operators
            op.setStencilSize(stencilSize);
            op.setOrderOfAccuracy(orderOfAccuracy);
      //   op.setTwilightZoneFlow(TRUE);
      // op.setTwilightZoneFlowFunction(exact);

            f.setOperators(op); // for apply the BC
            coeff.setOperators(op);
    
            BoundaryConditionParameters bcParams;

      // cout << "op.laplacianCoefficients().className: " << (op.laplacianCoefficients()).getClassName() << endl;
      // cout << "-op.laplacianCoefficients().className: " << (-op.laplacianCoefficients()).getClassName() << endl;
        
            Oges solver( cg );                     // create a solver
            
            #ifdef USE_PPP
            if( solverType==OgesParameters::PETScNew )
            {
                Oges::OGES_COMM_WORLD=myCommWorld;
                solver.setCommunicator( myComm );
            }
            #endif 
            
            solver.set(OgesParameters::THEsolverType,solverType); 


            if( true )
            {
        // Build a wide stencil
                usePredefined=false;

                formImplicitTimeSteppingMatrix( cg, coeff, solver, op   );

        // solver.setCoefficientArray( coeff );   // supply coefficients
            }
            else if( usePredefined )
            {
        // ------- use predefined equations ------

                IntegerArray boundaryConditions(2,3,cg.numberOfComponentGrids());
                RealArray bcData(2,2,3,cg.numberOfComponentGrids());
                boundaryConditions=OgesParameters::dirichlet;
                bcData=0.;
                solver.setEquationAndBoundaryConditions( OgesParameters::laplaceEquation,op,boundaryConditions,bcData);
            }
            else
            {
        // ----------- Define equations -----------

                if( false  )
                {
          // Do this for now -- optimized coeff's currently only for 4th order  *wdh* 2016/08/17
                    coeff=op.laplacianCoefficients();       // get the coefficients for the Laplace operator
                }
                else
                { // new way for parallel -- this avoids all communication
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                        op[grid].coefficients(MappedGridOperators::laplacianOperator,coeff[grid],I1,I2,I3);
                    }
                }

        // fill in the coefficients for the boundary conditions
                coeff.applyBoundaryConditionCoefficients(0,0,dirichlet,  allBoundaries);
                coeff.applyBoundaryConditionCoefficients(0,0,extrapolate,allBoundaries);
                bcParams.orderOfExtrapolation=orderOfAccuracy+1; // *wdh* 2016/08/17 -- check me
                if( orderOfAccuracy>=4 )
                {
                    bcParams.ghostLineToAssign=2;
                    coeff.applyBoundaryConditionCoefficients(0,0,extrapolate,allBoundaries,bcParams); // extrap 2nd ghost line
                    bcParams.ghostLineToAssign=1;  // reset *wdh* 2016/08/17
                }
                if( orderOfAccuracy>=6 )
                {
          // do this for now: *wdh* 2016/08/17
                    bcParams.ghostLineToAssign=3;
                    coeff.applyBoundaryConditionCoefficients(0,0,extrapolate,allBoundaries,bcParams); // extrap 3rd ghost line
                    bcParams.ghostLineToAssign=1;  // reset *wdh* 2016/08/17
                }
                coeff.finishBoundaryConditions();
        // coeff.display("Here is coeff after finishBoundaryConditions");

                if( false )
                {
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        displayCoeff(coeff[grid],sPrintF("Coeff matrix for grid %i",grid));
                    
            // coeff[grid].sparse->classify.display("the classify matrix after applying finishBoundaryConditions()");
            //    coeff[grid].display("this is the coefficient matrix");
                    }
                }
            } 
            

            
            if( outputMatrix )
                solver.set(OgesParameters::THEkeepSparseMatrix,true);
            
            if( solver.isSolverIterative() ) 
            {
                solver.setCommandLineArguments( argc,argv );

                if( iterativeSolverType=="gmres" )
                {
                    solver.set(OgesParameters::THEsolverMethod,OgesParameters::generalizedMinimalResidual);
                    solver.set(OgesParameters::THEpreconditioner,OgesParameters::incompleteLUPreconditioner);
                }
                else if( iterativeSolverType=="hypre" )
                {
          // NOTE: hypre is called through PETSc
          // NOTE: Hypre AMG is a PC type within a Kyrlov solver such as gmres or bcgs, 
                    solver.set(OgesParameters::THEparallelSolverMethod,OgesParameters::gmres);
                    solver.set(OgesParameters::THEparallelPreconditioner,OgesParameters::hyprePreconditioner);
                    solver.set(OgesParameters::THEpreconditioner,OgesParameters::hyprePreconditioner);
                    solver.set(OgesParameters::THEparallelExternalSolver,OgesParameters::hypre);

                    solver.parameters.setPetscOption("-ksp_type","gmres");
                    solver.parameters.setPetscOption("-pc_type","hypre");
                    solver.parameters.setPetscOption("-pc_hypre_type","boomeramg");
                    solver.parameters.setPetscOption("-pc_hypre_boomeramg_strong_threshold",".5");
                    solver.parameters.setPetscOption("-pc_hypre_boomeramg_max_levels","20");
                    solver.parameters.setPetscOption("-pc_hypre_boomeramg_coarsen_type","Falgout");

                }
                else
                {
                    if( solverType==OgesParameters::PETSc )
                        solver.set(OgesParameters::THEsolverMethod,OgesParameters::biConjugateGradientStabilized);
                    else if( solverType==OgesParameters::PETScNew )
                    { // parallel: -- NOTE: in parallel the solveMethod should be preonly and the parallelSolverMethod bicgs etc.
                        solver.set(OgesParameters::THEbestIterativeSolver);
            // solver.set(OgesParameters::THEparallelSolverMethod,OgesParameters::biConjugateGradient);
            // solver.set(OgesParameters::THEparallelSolverMethod,OgesParameters::gmres);
            // Use an LU solver on each processor:
            // solver.set(OgesParameters::THEpreconditioner,OgesParameters::luPreconditioner);
            // This also works: Use an LU on each processor:
            // solver.parameters.setPetscOption("-sub_pc_type","lu");
                    }
                    else
                        solver.set(OgesParameters::THEsolverMethod,OgesParameters::biConjugateGradient);
                }
                
                solver.set(OgesParameters::THErelativeTolerance,max(tol,REAL_EPSILON*10.));
                solver.set(OgesParameters::THEmaximumNumberOfIterations,10000);
                if( iluLevels>=0 )
                    solver.set(OgesParameters::THEnumberOfIncompleteLULevels,iluLevels);
            }    

            printF("\n === Solver:\n %s\n =====\n",(const char*)solver.parameters.getSolverName());

            if( false )
                solver.parameters.display();
            

      // ---------------------------------
      // --------- Dirichlet BC's --------
      // ---------------------------------
            if( problemsToSolve % 2 ==1 )
            {

                if( !usePredefined )
                    solver.setCoefficientArray( coeff );   // supply coefficients

        // Assign the right-hand-side f  
                assignForcing( 0,cg,f,exact );
              
                f.display("RHS f","%6.2f ");

        // f=0.; // *** TEST
    
                u=0.;  // initial guess for iterative solvers
                real time0=CPU();
                solver.solve( u,f );   // solve the equations
                real time= ParallelUtility::getMaxValue(CPU()-time0);
                printF("\n*** max residual=%8.2e, time for 1st solve of the Dirichlet problem = %8.2e (iterations=%i) ***\n",
                              solver.getMaximumResidual(),time,solver.getNumberOfIterations());

        // solve again
                if( true )
                {
          // u=0.;
                    time0=CPU();
                    solver.solve( u,f );   // solve the equations
                    time= ParallelUtility::getMaxValue(CPU()-time0);
                    printF("\n*** max residual=%8.2e, time for 2nd solve of the Dirichlet problem = %8.2e (iterations=%i) ***\n\n",
                                  solver.getMaximumResidual(),time,solver.getNumberOfIterations());
                }
                if( outputMatrix )
                {
                    printF("tcmWideStencil:INFO: save the matrix to file tcmWideStencilMatrix.out (using writeMatrixToFile). \n");
                    solver.writeMatrixToFile("tcmWideStencilMatrix.out");

                    aString fileName = "sparseMatrix.dat";
                    printF("tcmWideStencil:INFO: save the matrix to file %s (using outputSparseMatrix)\n",(const char*)fileName);
                    solver.outputSparseMatrix( fileName );
                }
            
            
        // ---- check the errors in the solution ---

                const int numberOfGridsPoints=max(1,cg.numberOfGridPoints());
                const real solverSize=solver.sizeOf();
                printF(".....solver: size = %8.2e (bytes), grid-pts=%i, reals/grid-pt=%5.2f \n",
                              solverSize,numberOfGridsPoints,solverSize/(numberOfGridsPoints*sizeof(real)));

        // u.display("Here is the solution to Laplacian(u)=f");
                computeTheError( 0,cg,u,err,exact, error );
                worstError=max(worstError,error);

                checker.setCutOff(errorBound[it][precision][0]); checker.printMessage("dirichlet: error",error,time);
    
                if( plot )
                {
                    ps.createWindow("tcmWideStencil");
                    plotResults( ps,solver,u,err );
                }
                
            }


      // ------------------------------------------
      // --------- Neumann or Mixed BC's ----------
      // ------------------------------------------
      // if( (problemsToSolve/2) % 2 ==1 ||
      //     (problemsToSolve/4) % 2 ==1 )
      // {
      //   bool neumannBCs =(problemsToSolve/2) % 2 ==1;
      //   bool mixedBCs   =(problemsToSolve/4) % 2 ==1;
      //   aString optionName = neumannBCs ? "neumann" : "mixed";
      //   RealArray *varCoeff=NULL;  // holds variable coefficients

      //   if( usePredefined )
      //   {
      //   // ------- use predefined equations ------

      //     IntegerArray boundaryConditions(2,3,cg.numberOfComponentGrids());
      //     RealArray bcData(2,2,3,cg.numberOfComponentGrids());
      //     boundaryConditions=OgesParameters::neumann;
      //     bcData=0.;
      //     solver.setEquationAndBoundaryConditions( OgesParameters::laplaceEquation,op,boundaryConditions,bcData);

      //   }
      //   else
      //   {
      //     // ----------- Define equations -----------


      //     coeff=0.;
      //     if( false && orderOfAccuracy>= 6  )
      //     {
      //       // Do this for now -- optimized coeff's currently only for 4th order  *wdh* 2016/08/17
      //       coeff=op.laplacianCoefficients();       // get the coefficients for the Laplace operator
      //     }
      //     else
      //     {
      //       // optimized way
      //       for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      //       {
      //         getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
      //         op[grid].coefficients(MappedGridOperators::laplacianOperator,coeff[grid],I1,I2,I3);
      //       }
      //     }
                
      //     // fill in the coefficients for the boundary conditions
      //     if( neumannBCs )
      //     {
      //       coeff.applyBoundaryConditionCoefficients(0,0,neumann,allBoundaries);
      //     }
      //     else
      //     {
      //       // -- mixed BC's with variable coefficients --

      //       bcParams.setVariableCoefficientOption(  BoundaryConditionParameters::spatiallyVaryingCoefficients );
      //       varCoeff = new RealArray [cg.numberOfComponentGrids()];
      //       for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      //       {
      //         MappedGrid & mg = cg[grid];
      //         int numGhost=1;
      //         getIndex(mg.gridIndexRange(),I1,I2,I3,numGhost);
      //         realArray & vertex = mg.vertex();
      //         OV_GET_SERIAL_ARRAY_CONST(real,vertex,x);
      //         int includeGhost=1;
      //         bool ok = ParallelUtility::getLocalArrayBounds(vertex,x,I1,I2,I3,includeGhost);
      //         if( ok ) 
      //         {
      //           // varCoeff only needs to be allocated on the boundary allocate on entire grid 
      //           // so we can assign all boundaries in one call (below)
      //           RealArray & vc = varCoeff[grid];
      //           vc.redim(I1,I2,I3,2);  // holds variable coefficients
      //           bcParams.setVariableCoefficientsArray( &vc );        

      //           // coeff of u 
      //           vc(I1,I2,I3,0)=1.+ .025*SQR(x(I1,I2,I3,0)) + .03*SQR(x(I1,I2,I3,1));   
      //           // coeff of u.n : (this value must not be zero)
      //           vc(I1,I2,I3,1)=2. + .1*SQR(x(I1,I2,I3,0)) + .05*SQR(x(I1,I2,I3,1)); 
      //         }
                        
      //         coeff[grid].applyBoundaryConditionCoefficients(0,0,mixed,allBoundaries,bcParams);

      //         // reset:
      //         bcParams.setVariableCoefficientsArray( NULL ); 
      //       }
      //       // reset: 
      //       bcParams.setVariableCoefficientOption( BoundaryConditionParameters::spatiallyConstantCoefficients );
                    
      //     }
                
      //     if( orderOfAccuracy>=4 )
      //     {
      //       bcParams.ghostLineToAssign=2;
      //       coeff.applyBoundaryConditionCoefficients(0,0,extrapolate,allBoundaries,bcParams); // extrap 2nd ghost line
      //       bcParams.ghostLineToAssign=1;  // reset

      //     }
                
      //     if( orderOfAccuracy>=6 )
      //     {
      //       // do this for now *wdh* 2016/08/17
      //       bcParams.ghostLineToAssign=3;
      //       coeff.applyBoundaryConditionCoefficients(0,0,extrapolate,allBoundaries,bcParams); // extrap 2nd ghost line
      //       bcParams.ghostLineToAssign=1;  // reset

      //     }
                
      //     coeff.finishBoundaryConditions();
      //     if( false )
      //     {
      //       for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      //         ::displayCoeff(coeff[grid],sPrintF("coeff on grid=%i",grid));
      //     }
                

      //     solver.setCoefficientArray( coeff );   // supply coefficients

      //   } // end if !usePredefined 
                
      //   bool singularProblem=neumannBCs; 
      //   // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      //   // { // this loop does nothing for now 
      //   //   MappedGrid & mg = cg[grid];
      //   //   ForBoundary(side,axis)
      //   //   {
      //   //     if( mg.boundaryCondition(side,axis) > 0  )
      //   //     { 
      //   //     }
      //   //     else if( mg.boundaryCondition(side,axis) ==inflow ||  mg.boundaryCondition(side,axis) ==outflow )
      //   //     {
      //   //       singularProblem=false;
      //   //     }
      //   //   }
      //   // }

      //   // Assign the right-hand-side f  
      //   const int option = neumannBCs ? 1 : 2;
      //   assignForcing( option,cg,f,exact,varCoeff );

      //   delete [] varCoeff;

      //   // if the problem is singular Oges will add an extra constraint equation to make the system nonsingular
      //   if( singularProblem )
      //     solver.set(OgesParameters::THEcompatibilityConstraint,TRUE);

      //   // Tell the solver to refactor the matrix since the coefficients have changed
      //   solver.setRefactor(TRUE);
      //   // we need to reorder too because the matrix changes a lot for the singular case
      //   solver.setReorder(TRUE);
                
      //   if( singularProblem )
      //   {
      //     // we need to first initialize the solver before we can fill in the rhs for the compatibility equation
      //     solver.initialize();
      //     realCompositeGridFunction ue(cg);
      //     exact.assignGridFunction(ue,0.);
      //     real value=0.;
      //     solver.evaluateExtraEquation(ue,value);

      //     if( Oges::debug & 4 )
      //       printF(" Neumann: RHS for singular equation=%10.3e\n",value);
                    
      //     solver.setExtraEquationValues(f,&value );
      //   }
      //   if( Oges::debug & 4 )
      //   {
      //     f.display("RHS before solve");
      //   }
                
      //   u=0.;  // initial guess for iterative solvers
      //   real time0=CPU();
      //   solver.solve( u,f );   // solve the equations
      //   real time= ParallelUtility::getMaxValue(CPU()-time0);
      //   printF("\n*** residual=%8.2e, time for 1st solve of the %s problem = %8.2e (iterations=%i)\n",
      //          solver.getMaximumResidual(),(const char*)optionName,time,solver.getNumberOfIterations());

      //   // turn off refactor for the 2nd solve
      //   solver.setRefactor(FALSE);
      //   solver.setReorder(FALSE);
      //   // u=0.;  // initial guess for iterative solvers
      //   time0=CPU();
      //   solver.solve( u,f );   // solve the equations
      //   time= ParallelUtility::getMaxValue(CPU()-time0);
      //   printF("\n*** residual=%8.2e, time for 2nd solve of the %s problem = %8.2e (iterations=%i)\n\n",
      //          solver.getMaximumResidual(),(const char*)optionName, time,solver.getNumberOfIterations());

      //   if( outputMatrix )
      //   {
      //     printF("tcmWideStencil:INFO: save the matrix to file tcmWideStencilMatrix.out (using writeMatrixToFile). \n");
      //     solver.writeMatrixToFile("tcmWideStencilMatrixNeumann.out");

      //     aString fileName = "sparseMatrixNeumann.dat";
      //     printF("tcmWideStencil:INFO: save the matrix to file %s (using outputSparseMatrix)\n",(const char*)fileName);
      //     solver.outputSparseMatrix( fileName );
      //   }
            

      //   computeTheError( option,cg,u,err,exact, error );

      //   worstError=max(worstError,error);
            
      //   checker.setCutOff(errorBound[it][precision][1]);
      //   aString buff;
      //   checker.printMessage(sPrintF(buff,"%s: error",(const char*)optionName),error,time);
            
      // } // end Neumann
            
            if( plot )
            {
                if( !( problemsToSolve % 2 ==1 ))
                    ps.createWindow("tcmWideStencil");
                plotResults( ps,solver,u,err );
            }

            delete exactPointer; exactPointer=0;// kkc 090902, this was a memory leak making new OGFunction's for each grid w/o releasing the previous one

        }  // end it (number of grids)
        

    }  // end sparseSolver

#ifdef USE_PPP
    if( communicatorWasNewed )
    {
        if( myComm !=MPI_COMM_NULL )
            MPI_Comm_free(&myComm);
        if( myCommWorld !=MPI_COMM_NULL )
            MPI_Comm_free(&myCommWorld);
    }
#endif
    
    fflush(0);
    printF("\n\n ************************************************************************************************\n");
    if( worstError > .025 )
        printF(" ************** Warning, there is a large error somewhere, worst error =%e ******************\n",
                      worstError);
    else
        printF(" ************** Test apparently successful, worst error =%e ******************\n",worstError);
    printF(" **************************************************************************************************\n\n");

    fflush(0);

    Overture::finish();          

    return(0);
}




