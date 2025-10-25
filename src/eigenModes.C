// This file automatically generated from eigenModes.bC with bpp.
// =====================================================================================
// Routines to compute eigenvalues and eigenvectors
// =====================================================================================
#include "CgWave.h"
#include "ParallelUtility.h"
#include "gridFunctionNorms.h"
#include "Integrate.h"
#include "CompositeGridOperators.h"   

#include "PlotStuff.h"
#include "GL_GraphicsInterface.h" 

// lapack routines
#ifdef OV_USE_DOUBLE
    #define SYEV EXTERN_C_NAME(dsyev)
    #define GEEV EXTERN_C_NAME(dgeev)
#else
    #define SYEV EXTERN_C_NAME(ssyev)
    #define GEEV EXTERN_C_NAME(sgeev)
#endif

// lapack routines
#ifdef OV_USE_DOUBLE
    #define GETRF EXTERN_C_NAME(dgetrf)
    #define GETRS EXTERN_C_NAME(dgetrs)
#else
    #define GETRF EXTERN_C_NAME(sgetrf)
    #define GETRS EXTERN_C_NAME(sgetrs)
#endif

extern "C"
{
  // Eigenvalues/eigenvectors of a real symmetric matrix:
    void SYEV( const char *jobz, const char* uplo, int & n, real & a, const int & lda,
                          real & w, real & work, int & lwork, int & info );

    void GEEV( char *jobvl, char* jobvr, int & n, real & a, const int & lda,
                          real & wr, real & wi, real &vl, int & ldvl, real & vr, int & ldvr, real & work, int & lwork, int & info );

  // PA = LU factor
    void GETRF( int & m, int & n, Real & a, const int & lda, int & ipvt, int & info );
  // Solve given LU
    void GETRS( char *trans, int & n, int & nhrs, Real & a, const int & lda, int & ipvt, Real & b, const int & ldb, int & info );

}

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )


// --------------------------------------------------------------------------------------
//   Macro: return the index's for possible active points
//            
//  NOTE: This macro appears in solveSLEPc.bC and eigenModes.bC 
// --------------------------------------------------------------------------------------
    

// ===============================================================
// Macro: Compute the inner product of two grid functions          
// ===============================================================

// ===============================================================
// Macro: Gram-Schmidt orthogonalization
// ===============================================================

// =================================================================================
/// \brief Update the eigenmodes: normalize solution, apply Rayeigh-Ritz.
///
// =================================================================================
int CgWave::updateEigenmodes()
{
    const int & iteration = dbase.get<int>("iteration");

    const int & numberOfFrequencies     = dbase.get<int>("numberOfFrequencies");
    RealArray & frequencyArray          = dbase.get<RealArray>("frequencyArray");
    RealArray & frequencyArrayAdjusted  = dbase.get<RealArray>("frequencyArrayAdjusted");
    const int & computeEigenmodes       = dbase.get<int>("computeEigenmodes");
    RealArray & resVector               = dbase.get<RealArray>("resVector");  // WaveHoltz residual
    const int & assignRitzFrequency     = dbase.get<int>("assignRitzFrequency");
    const int & useAccurateInnerProduct = dbase.get<int>("useAccurateInnerProduct");

    const Real residual = iteration>0 ? resVector(iteration-1) : 1.;
    printF("Cgwave::updateEigenModes:: normalize eigenmodes, iteration=%d, WaveHoltz residual=%8.2e\n",iteration,residual);

    assert( computeEigenmodes );

  // save computed eigenvalues here 
    const int maxIterations=500; // ** FIX ME**
    if( !dbase.has_key("eigenValues") )
    {
        RealArray & eigenValues = dbase.put<RealArray>("eigenValues");
        eigenValues.redim(numberOfFrequencies);
        eigenValues=-1.;

    // -- keep iteration history: 

        RealArray & eigsRayleighQuotient = dbase.put<RealArray>("eigsRayleighQuotient");
        eigsRayleighQuotient.redim(maxIterations);
        eigsRayleighQuotient=-1.;

        RealArray & eigsRayleighRitz = dbase.put<RealArray>("eigsRayleighRitz");
        eigsRayleighRitz.redim(maxIterations);
        eigsRayleighRitz=-1.;    
    }

    if( !dbase.has_key("numEigenVectors" ))
        dbase.put<int>("numEigenVectors")=0;

    int & numEigenVectors            = dbase.get<int>("numEigenVectors");             // number of eigenPairs
    RealArray & eigenValues          = dbase.get<RealArray>("eigenValues");           // current guess at eigenvalues
    RealArray & eigsRayleighQuotient = dbase.get<RealArray>("eigsRayleighQuotient");  // eigenvalues by iteration 
    RealArray & eigsRayleighRitz     = dbase.get<RealArray>("eigsRayleighRitz");

    numEigenVectors=1;  // For now we just compute one eigenPair


    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");

  // ----- NORMALIZE v ----
    RealArray vNorm(numberOfFrequencies);
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
        vNorm(freq) = maxNorm(v,freq);
    }
    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
        getIndex(cg[grid].dimension(),I1,I2,I3);
        bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
        if( ok )
        {
            for( int freq=0; freq<numberOfFrequencies; freq++ )
            {
                vLocal(I1,I2,I3,freq) /= vNorm(freq);
            }
        }

    }

  // estimate the eigenvalue using a Rayleigh quotient 
  // lambda = (vk^T L vk)/( vk^T vk)


    printF("Cgwave::updateEigenModes:: estimate the eigenvalue using the Rayleigh quotient...\n");
    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
    

    if( !dbase.has_key("integrate") )
    {
        Integrate & integrate = dbase.put<Integrate>("integrate");
        integrate.updateToMatchGrid(cg);
    }
    Integrate & integrate = dbase.get<Integrate>("integrate");
    realCompositeGridFunction lap(cg);   // ***** do this for now ... is there a work space we can use instead?

  // **** FIX ME : use getRayleighQuotient ************************************************************
    if( 1==1 )
    {
        Real lambda = getRayleighQuotient( v );
        eigenValues(0)=lambda;

      if( iteration<maxIterations )
          eigsRayleighQuotient(iteration)=lambda;
    }
    else
    {
    // *** OLD WAY 

    // for( int freq=0; freq<numberOfFrequencies; freq++ )
    // {

    //   // Compute   lap = L v 
    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
    //     OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);

    //     getIndex(cg[grid].dimension(),I1,I2,I3);

    //     operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lapLocal,I1,I2,I3,freq);

    //   }

    //   // Compute ( v, vLap )
    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     MappedGrid & mg = cg[grid];
    //     getIndex(mg.dimension(),I1,I2,I3);
    //     OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
    //     OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);      

    //     lapLocal(I1,I2,I3) = vLocal(I1,I2,I3,freq)*lapLocal(I1,I2,I3);
    //   }

    //   // *** vLocal should have zero BC's -- so no need to change lap I think ****
    //   lap.periodicUpdate(); // is this needed ? interpolate may do this
    //   lap.interpolate();  


    //   Real vDotLap = integrate.volumeIntegral(lap);    

    //   // Compute ( v,v )
    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     MappedGrid & mg = cg[grid];
    //     getIndex(mg.dimension(),I1,I2,I3);
    //     OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
    //     OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);  

    //     lapLocal(I1,I2,I3) = vLocal(I1,I2,I3,freq)*vLocal(I1,I2,I3,freq); // store in lapLocal
    //   }

    //   Real vDotv = integrate.volumeIntegral(lap);  

    //   Real lambda = vDotLap/vDotv;
    //   if( lambda<0 )
    //   {
    //     lambda=sqrt(-lambda);
    //     printF("Rayleigh quotient: freq=%d: lambda(RQ)=%16.8e (omega=%16.8e) (Lap(phi) = -lambda^2 phi)\n",freq,lambda,frequencyArray(freq));
    //   }
    //   else
    //   {
    //     printF("Cgwave::updateEigenModes::ERROR: rayleigh quotient for Laplacian is postive!? lambda=%12.4e\n",lambda);
    //     OV_ABORT("error");
    //   }

    //   eigenValues(freq)=lambda;

    //   if( iteration<maxIterations )
    //     eigsRayleighQuotient(iteration)=lambda;
    // }

    } // end OLD WAY 



    
  // -----------------------------------------------------
  // --------- Rayleigh Ritz Orthogonal Projection -------
  // -----------------------------------------------------

  // Given a set of orthonormal vectors qj
  //       Qm = [ q1 q2 ... qm ]
  // Find the eigenvalues to the Galerkin projected equations
  //       Bm = Qm^* A Qm 
  //
    if( !dbase.has_key("iterationStartRR") ) 
        dbase.put<int>("iterationStartRR") = 5;

    int & iterationStartRR = dbase.get<int>("iterationStartRR");      // start at this WH iteration

  // max number of Rayleigh-Ritz vectors: 
    const int maxNumberOfVectorsRR = dbase.get<int>("numberOfRitzVectors");
  // const int maxNumberOfVectorsRR = 30; // 5; // 30;  // 5;     

    if( !dbase.has_key("iterationRR" ) )
    {
        dbase.put<int>("iterationRR")=0;
    }

    int & iterationRR = dbase.get<int>("iterationRR");


    if( !dbase.has_key("eigenVectorRayleighRitz" ) )
    {
        realCompositeGridFunction & eigenVectorRayleighRitz = dbase.put<realCompositeGridFunction>("eigenVectorRayleighRitz");
        eigenVectorRayleighRitz.updateToMatchGrid(cg);
    }
    realCompositeGridFunction & eigenVectorRayleighRitz = dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");
    

  // ---- Only start the Rayleigh Ritz when the residual goes below a given tolerance
    const Real residTolRR=1.e-2;
    if( residual>residTolRR )
        iterationStartRR=iteration+1;

  // if( iteration>= iterationStartRR && iterationRR<maxNumberOfVectorsRR )
    if( iteration>= iterationStartRR )
    {
        printF("---- updateEigenModes: Rayleigh-Ritz-iteration=%d (iteration=%d)----\n",iterationRR,iteration);

        int freq=0;
        assert( numberOfFrequencies==1 );

        Range all;
        if( !dbase.has_key("Q") )
        {
              dbase.put<int>("numVectorsSaved") = 0;
              dbase.put<int>("currentVector")   =-1;
              RealCompositeGridFunction & vSave = dbase.put<RealCompositeGridFunction>("vSave");  // save past solutions
              vSave.updateToMatchGrid(cg,all,all,all,maxNumberOfVectorsRR);

              RealCompositeGridFunction & Q     = dbase.put<RealCompositeGridFunction>("Q");  // holds orthonormal basis for Krylov space
              Q.updateToMatchGrid(cg,all,all,all,maxNumberOfVectorsRR);

              RealArray & Bm = dbase.put<RealArray>("Bm");  // holds B matrix 
              Bm.redim(maxNumberOfVectorsRR,maxNumberOfVectorsRR);
        }
        int & numVectorsSaved             = dbase.get<int>("numVectorsSaved");
        int & currentVector               = dbase.get<int>("currentVector");
        RealCompositeGridFunction & vSave = dbase.get<RealCompositeGridFunction>("vSave");
        RealCompositeGridFunction & Q     = dbase.get<RealCompositeGridFunction>("Q");
        RealArray & Bm                    = dbase.get<RealArray>("Bm");    

        RealCompositeGridFunction w(cg,all,all,all);  // work space
        w=0.;                                         // inner product ignores ghost so make sure these are zero

    // --------------------------------------------------
    // ----- START Gram-Schmidt orthogonalization -------
    // --------------------------------------------------

        currentVector = (currentVector+1) % maxNumberOfVectorsRR;
        numVectorsSaved = min(numVectorsSaved+1,maxNumberOfVectorsRR);

        int numberOfVectorsToUse = numVectorsSaved; // Keep this many Q vectors (eliminate those that are linear dependent)

    // ---------- save the current v in a circular list -------------
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            vSave[grid](all,all,all,currentVector) = v[grid](all,all,all,freq);

        bool useHouseholder=false; // true;
        bool useReverseOrder=true;
            if( useHouseholder )
            {
                printF("RR: Use HOUSEHOLDER ORTHOGONALIZATION\n");
                if( cg.numberOfComponentGrids()>1 )
                {
                    printF("HOUSEHOLDER: finish me for more than one component grid\n");
                    OV_ABORT("finish me");
                }
                if( !dbase.has_key("vHouse") )
                {
                      RealCompositeGridFunction & vHouse = dbase.put<RealCompositeGridFunction>("vHouse");  // save Householder vectors here 
                      vHouse.updateToMatchGrid(cg,all,all,all,maxNumberOfVectorsRR);
                }
                RealCompositeGridFunction & vHouse = dbase.get<RealCompositeGridFunction>("vHouse");  
        //   Householder computes
        //       Qn * ... *Q2 *Q1 A = R     (we don't use R)
        // using Householder reflectors: 
        //   Qj = I - 2 vj vj^*
        // 
        // To start Q holds our matrix A = [ a1 a2 .. an ]
                for( int ia=0; ia<numVectorsSaved; ia++ )
                {
                    int next;
                    if( useReverseOrder )
                        next = ( currentVector-ia + maxNumberOfVectorsRR) % maxNumberOfVectorsRR; // use in reverse order
                    else
                        next = ( ia + maxNumberOfVectorsRR) % maxNumberOfVectorsRR; // forward order      
          // int next = ( currentVector - ia + maxNumberOfVectorsRR) % maxNumberOfVectorsRR; // reverse order
          // int next = ( currentVector + ia + maxNumberOfVectorsRR) % maxNumberOfVectorsRR; // forward order
          // int next=ia; 
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        Q[grid](all,all,all,ia) = vSave[grid](all,all,all,next);
                    }
                }
        // Make a list of the first active grid points *** finish me***
        // for now use i1=1,2,3,..
        //             i2=0, i3=0; 
                int i1=0, i2=0, i3=0, grid0=0;
                int i1Start=0; 
        // ---- Compute the scalings for the unit vectors so that they have norm 1 ----
                RealArray unitVectorValue(numVectorsSaved);
                for( int ia=0; ia<numVectorsSaved; ia++ )
                {
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        vHouse[grid](all,all,all,ia) = 0.; 
                    i1 = i1Start + ia;
                    vHouse[grid0](i1,i2,i3,ia)=1.; // unit vector (before scaling)
                    Real eNorm2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ia);
                            }
                            eNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            eNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            eNorm2 += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ia);
                                    }
                                }
                            }
                            eNorm2 = ParallelUtility::getSum( eNorm2 );
                        }  
                    assert( eNorm2>=0. );
                    eNorm2=sqrt(eNorm2);
                    unitVectorValue(ia)= 1./eNorm2;
                    printF("Unit vector ia=%d has norm = %9.3e, unitVectorValue=%9.3e\n",eNorm2,unitVectorValue(ia));
                }
                for( int ia=0; ia<numVectorsSaved; ia++ )
                {
          // Form Householder vector vj
          //   x = aj
          //   vj = x + sign(x1)*norm(x)
          //   vj = vj/norm(vj) 
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        vHouse[grid](all,all,all,ia) = Q[grid](all,all,all,ia);  // current column in altered "A"  = "x" 
                    }
          // printF("ia=%d: xNorm2=%9.3e\n",ia,xNorm2);
          // The Householder vectors get shorter as the algorithm progresses 
          // Just zero out the top entries in the column to get the same effect.
                    for( int j1=0; j1<ia; j1++ )
                    {
                        i1 = i1Start + j1;
                        vHouse[grid0](i1,i2,i3,ia)=0.; // zero out entries in column above current digaonal 
                    }
                    Real xNorm2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ia);
                            }
                            xNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            xNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            xNorm2 += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ia);
                                    }
                                }
                            }
                            xNorm2 = ParallelUtility::getSum( xNorm2 );
                        }  
                    assert( xNorm2>=0. );
                    xNorm2=sqrt(xNorm2);
                    i1 = i1Start + ia;  // diagonal entry --  do this for now ******FIX ME ****************************************************************
                    Real x1 = vHouse[grid0](i1,i2,i3,ia); 
                    Real signx1 = x1>=0. ? 1. : -1;
                    vHouse[grid0](i1,i2,i3,ia) += signx1*xNorm2*unitVectorValue(ia); // (i1,i2,i3,grid0)  must be a valid point in the inner product
                    Real vNorm2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ia);
                            }
                            vNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            vNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            vNorm2 += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ia);
                                    }
                                }
                            }
                            vNorm2 = ParallelUtility::getSum( vNorm2 );
                        }  
                    assert( vNorm2>=0 );
                    vNorm2=sqrt(vNorm2);
          // printF("ia=%d: x1=%9.3e, vNorm2=%9.3e\n",ia,x1,vNorm2);
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        vHouse[grid](all,all,all,ia) *= (1./vNorm2);
          // Now multiply Qj =  I - 2 vj vj^* times each remaining column of A 
          // for( int ja=ia+1; ja<numVectorsSaved; ja++ )
                    for( int ja=ia; ja<numVectorsSaved; ja++ )  // try this ??
                    {
                        Real vDotQja;
                            if( useAccurateInnerProduct )
                            {
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                    OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ja);
                                }
                                vDotQja = integrate.volumeIntegral(w);
                            }
                            else
                            {
                                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                int iab[2]; 
                                int ii1,ii2,ii3;  
                                vDotQja=0.; 
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    MappedGrid & mg = cg[grid];
                                    const IntegerArray & gid = mg.gridIndexRange();
                                    OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                        Iv[2]=Range(0,0);
                                        for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                        {
                                            for( int side=0; side<=1; side++ )
                                            {
                                                int is = 1-2*side;
                                                iab[side]=gid(side,axis);
                                                const int bc = mg.boundaryCondition(side,axis);
                                                if( bc==CgWave::dirichlet )
                                                {
                                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                                }
                                                else if( bc==CgWave::neumann )
                                                {
                          // include boundary
                                                }
                                                else if( bc>0 )
                                                {
                                                    printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                    OV_ABORT("error");
                                                }
                                                else if( bc<0 )
                                                {
                          // periodic -- include left end
                                                    if( side==1 )
                                                        iab[side] += is; 
                                                }
                                                else
                                                {
                          // interpolation boundary : include end 
                                                }
                                            }
                                            Iv[axis] = Range(iab[0],iab[1]);
                                        }
                                    bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                    if( ok )
                                    {
                                        FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                        {
                                            if( maskLocal(ii1,ii2,ii3)>0 )
                                                vDotQja += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ja);
                                        }
                                    }
                                }
                                vDotQja = ParallelUtility::getSum( vDotQja );
                            }  
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            Q[grid](all,all,all,ja) -= (2.*vDotQja) * vHouse[grid](all,all,all,ia);
                        }
                    }
                }
        // Form the matrix of Q = Q1 * Q2 * ... * Qn
        // Form column j of Q from 
        //         Q1*Q2* ... Qn* * e_j , j=0,1,2,n-1
        //  where
        //       e_j = unit vector in direction j
                for( int ia=0; ia<numVectorsSaved; ia++ )
                {
          // Q[grid](all,all,all,ia) = Unit vector e[ia];
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        Q[grid](all,all,all,ia) = 0.;
                    i1 = i1Start + ia;                 // do this for now ******FIX ME ******************************************
          // Put "1" into Q for a unit vector but we must rescale below 
          // since the unit vect depends on the inner product
                    Q[grid0](i1,i2,i3,ia) = unitVectorValue(ia);  // this must be a valid point in the inner product
          // Multiply Qj in reverse order
                    for( int ja=numVectorsSaved-1; ja>=0; ja-- )
                    {
                        Real vDotQ;
                            if( useAccurateInnerProduct )
                            {
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                    OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ja)*vvLocal(I1,I2,I3,ia);
                                }
                                vDotQ = integrate.volumeIntegral(w);
                            }
                            else
                            {
                                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                int iab[2]; 
                                int ii1,ii2,ii3;  
                                vDotQ=0.; 
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    MappedGrid & mg = cg[grid];
                                    const IntegerArray & gid = mg.gridIndexRange();
                                    OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                        Iv[2]=Range(0,0);
                                        for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                        {
                                            for( int side=0; side<=1; side++ )
                                            {
                                                int is = 1-2*side;
                                                iab[side]=gid(side,axis);
                                                const int bc = mg.boundaryCondition(side,axis);
                                                if( bc==CgWave::dirichlet )
                                                {
                                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                                }
                                                else if( bc==CgWave::neumann )
                                                {
                          // include boundary
                                                }
                                                else if( bc>0 )
                                                {
                                                    printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                    OV_ABORT("error");
                                                }
                                                else if( bc<0 )
                                                {
                          // periodic -- include left end
                                                    if( side==1 )
                                                        iab[side] += is; 
                                                }
                                                else
                                                {
                          // interpolation boundary : include end 
                                                }
                                            }
                                            Iv[axis] = Range(iab[0],iab[1]);
                                        }
                                    bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                    if( ok )
                                    {
                                        FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                        {
                                            if( maskLocal(ii1,ii2,ii3)>0 )
                                                vDotQ += uuLocal(ii1,ii2,ii3,ja)*vvLocal(ii1,ii2,ii3,ia);
                                        }
                                    }
                                }
                                vDotQ = ParallelUtility::getSum( vDotQ );
                            }  
            // printF("House: form Q[%d]: ja=%d, vDotQ=%9.3e\n",ia,ja,vDotQ);
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            Q[grid](all,all,all,ia) -= (2.*vDotQ) * vHouse[grid](all,all,all,ja);
                    }
          // Check:
                    Real qNorm2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ia);
                            }
                            qNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            qNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(Q[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            qNorm2 += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ia);
                                    }
                                }
                            }
                            qNorm2 = ParallelUtility::getSum( qNorm2 );
                        }  
                    qNorm2 = sqrt(qNorm2); 
                    printF("ia=%d: qNorm2=%9.3e  (SHOULD BE 1 NOW)\n",ia,qNorm2);  
          // Normalize Q 
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        Q[grid](all,all,all,ia) *= (1./qNorm2);
                }
        // -- check that Q[0] = vSave[0]/norm(vSave[0])
                if( true )
                {
                    int ia=0;
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        vHouse[grid](all,all,all,ia) = vSave[grid](all,all,all,ia);
                    Real vNorm2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ia);
                            }
                            vNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            vNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            vNorm2 += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ia);
                                    }
                                }
                            }
                            vNorm2 = ParallelUtility::getSum( vNorm2 );
                        }  
                    assert( vNorm2>=0 );
                    vNorm2=sqrt(vNorm2);
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        vHouse[grid](all,all,all,ia) *= (1./vNorm2);
                    int ib=1; 
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        vHouse[grid](all,all,all,ib) = vHouse[grid](all,all,all,ia) - Q[grid](all,all,all,ia);
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ib)*vvLocal(I1,I2,I3,ib);
                            }
                            vNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            vNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            vNorm2 += uuLocal(ii1,ii2,ii3,ib)*vvLocal(ii1,ii2,ii3,ib);
                                    }
                                }
                            }
                            vNorm2 = ParallelUtility::getSum( vNorm2 );
                        }  
                    assert( vNorm2>=0 );
                    vNorm2=sqrt(vNorm2);  
                    printF(" norm( Q[0] - vSave[0]/norm(vSave[0]) = %8.2e\n",vNorm2);
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        vHouse[grid](all,all,all,ib) = vHouse[grid](all,all,all,ia) + Q[grid](all,all,all,ia);
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ib)*vvLocal(I1,I2,I3,ib);
                            }
                            vNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            vNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,vHouse[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(vHouse[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            vNorm2 += uuLocal(ii1,ii2,ii3,ib)*vvLocal(ii1,ii2,ii3,ib);
                                    }
                                }
                            }
                            vNorm2 = ParallelUtility::getSum( vNorm2 );
                        }  
                    assert( vNorm2>=0 );
                    vNorm2=sqrt(vNorm2);  
                    printF(" norm( Q[0] + vSave[0]/norm(vSave[0]) = %8.2e\n",vNorm2);
          // OV_ABORT("stop here for now");    
                }
            }
            else if( true )
            {
        // 
        //  --- Modified Gram-Schmidt starting with the newest vector first, (i.e. use in reverse order) ---
        // 
        //  Thus the best v is kept first as q1 
                for( int ia=0; ia<numVectorsSaved; ia++ )
                {
                    int next;
                    if( useReverseOrder )
                        next = ( currentVector-ia + maxNumberOfVectorsRR) % maxNumberOfVectorsRR; // use in reverse order
                    else
                        next = ( ia + maxNumberOfVectorsRR) % maxNumberOfVectorsRR; // forward order
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        Q[grid](all,all,all,ia) = vSave[grid](all,all,all,next);        
                    for( int ja=0; ja<ia; ja++ )
                    {
            // subtract off components of previous q vectors (modified Gram Schmidt -- use current Q[ia] instead of v)
                        Real dotProduct=0.;
                            if( useAccurateInnerProduct )
                            {
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ja)*vvLocal(I1,I2,I3,ia);
                                }
                                dotProduct = integrate.volumeIntegral(w);
                            }
                            else
                            {
                                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                int iab[2]; 
                                int ii1,ii2,ii3;  
                                dotProduct=0.; 
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    MappedGrid & mg = cg[grid];
                                    const IntegerArray & gid = mg.gridIndexRange();
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                        Iv[2]=Range(0,0);
                                        for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                        {
                                            for( int side=0; side<=1; side++ )
                                            {
                                                int is = 1-2*side;
                                                iab[side]=gid(side,axis);
                                                const int bc = mg.boundaryCondition(side,axis);
                                                if( bc==CgWave::dirichlet )
                                                {
                                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                                }
                                                else if( bc==CgWave::neumann )
                                                {
                          // include boundary
                                                }
                                                else if( bc>0 )
                                                {
                                                    printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                    OV_ABORT("error");
                                                }
                                                else if( bc<0 )
                                                {
                          // periodic -- include left end
                                                    if( side==1 )
                                                        iab[side] += is; 
                                                }
                                                else
                                                {
                          // interpolation boundary : include end 
                                                }
                                            }
                                            Iv[axis] = Range(iab[0],iab[1]);
                                        }
                                    bool ok=ParallelUtility::getLocalArrayBounds(Q[grid],uuLocal,I1,I2,I3);
                                    if( ok )
                                    {
                                        FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                        {
                                            if( maskLocal(ii1,ii2,ii3)>0 )
                                                dotProduct += uuLocal(ii1,ii2,ii3,ja)*vvLocal(ii1,ii2,ii3,ia);
                                        }
                                    }
                                }
                                dotProduct = ParallelUtility::getSum( dotProduct );
                            }  
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            Q[grid](all,all,all,ia) -= dotProduct*Q[grid](all,all,all,ja);
                    }
          // Normalize Q[ia]
                    Real qNorm2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ia);
                            }
                            qNorm2 = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            qNorm2=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(Q[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            qNorm2 += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ia);
                                    }
                                }
                            }
                            qNorm2 = ParallelUtility::getSum( qNorm2 );
                        }  
                    if( qNorm2 > REAL_EPSILON*10. )
                        qNorm2 = sqrt(qNorm2);
                    else
                    {
            // Q's are linearly dependent **WE SHOULD NOT USE THIS LATEST Q ***    **** FIX ME ****
                        printF("RayleighRitz: INFO: Q vector ia=%d has 2-norm=%8.2e !! (Ignoring remaing vectors)\n",ia);
                        qNorm2 = 1.;
                        numberOfVectorsToUse = min(numberOfVectorsToUse,ia);
                        break;
                    }
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        Q[grid](all,all,all,ia) *= (1./qNorm2);
                }
            }
            else
            {
        // --- Gram-Schmidt with last latest vector last ----
                const int ia = iterationRR;
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    Q[grid](all,all,all,ia) = v[grid](all,all,all,freq);
                for( int ja=0; ja<ia; ja++ )
                {
          // subtract off components of previous q vectors (modified Gram Schmidt -- use current Q[ia] instead of v)
                    Real dotProduct=0.;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ja)*vvLocal(I1,I2,I3,ia);
                            }
                            dotProduct = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            dotProduct=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(Q[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            dotProduct += uuLocal(ii1,ii2,ii3,ja)*vvLocal(ii1,ii2,ii3,ia);
                                    }
                                }
                            }
                            dotProduct = ParallelUtility::getSum( dotProduct );
                        }  
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        Q[grid](all,all,all,ia) -= dotProduct*Q[grid](all,all,all,ja);
                }
        // Normalize Q[ia]
                Real qNorm2;
                    if( useAccurateInnerProduct )
                    {
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ia);
                        }
                        qNorm2 = integrate.volumeIntegral(w);
                    }
                    else
                    {
                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                        int iab[2]; 
                        int ii1,ii2,ii3;  
                        qNorm2=0.; 
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            const IntegerArray & gid = mg.gridIndexRange();
                            OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                Iv[2]=Range(0,0);
                                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                {
                                    for( int side=0; side<=1; side++ )
                                    {
                                        int is = 1-2*side;
                                        iab[side]=gid(side,axis);
                                        const int bc = mg.boundaryCondition(side,axis);
                                        if( bc==CgWave::dirichlet )
                                        {
                                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                        }
                                        else if( bc==CgWave::neumann )
                                        {
                      // include boundary
                                        }
                                        else if( bc>0 )
                                        {
                                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                            OV_ABORT("error");
                                        }
                                        else if( bc<0 )
                                        {
                      // periodic -- include left end
                                            if( side==1 )
                                                iab[side] += is; 
                                        }
                                        else
                                        {
                      // interpolation boundary : include end 
                                        }
                                    }
                                    Iv[axis] = Range(iab[0],iab[1]);
                                }
                            bool ok=ParallelUtility::getLocalArrayBounds(Q[grid],uuLocal,I1,I2,I3);
                            if( ok )
                            {
                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                {
                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                        qNorm2 += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ia);
                                }
                            }
                        }
                        qNorm2 = ParallelUtility::getSum( qNorm2 );
                    }  
                if( qNorm2 > REAL_EPSILON*10. )
                    qNorm2 = sqrt(qNorm2);
                else
                {
                    printF("RayleighRitz: INFO: Q vector ia=%d has 2-norm=%8.2e !!\n",ia);
                    qNorm2 = 1.;
                }
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    Q[grid](all,all,all,ia) *= (1./qNorm2);
            }
      // **** check orthogonality:
            if( true )
            {
                printF("\n ***** Rayleigh-Ritz CHECK GRAM-SCHIMDT ORTHOGONALIZATION : numberOfVectorsToUse=%d *****\n",numberOfVectorsToUse);
                Real maxErr=0;
                for( int ia=0; ia<numberOfVectorsToUse; ia++ )
                {
                    for( int ja=ia; ja<numberOfVectorsToUse; ja++ )
                    {
                        Real dotProduct;
                            if( useAccurateInnerProduct )
                            {
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ia)*vvLocal(I1,I2,I3,ja);
                                }
                                dotProduct = integrate.volumeIntegral(w);
                            }
                            else
                            {
                                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                int iab[2]; 
                                int ii1,ii2,ii3;  
                                dotProduct=0.; 
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    MappedGrid & mg = cg[grid];
                                    const IntegerArray & gid = mg.gridIndexRange();
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],uuLocal);
                                    OV_GET_SERIAL_ARRAY(Real,Q[grid],vvLocal);
                                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                        Iv[2]=Range(0,0);
                                        for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                        {
                                            for( int side=0; side<=1; side++ )
                                            {
                                                int is = 1-2*side;
                                                iab[side]=gid(side,axis);
                                                const int bc = mg.boundaryCondition(side,axis);
                                                if( bc==CgWave::dirichlet )
                                                {
                                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                                }
                                                else if( bc==CgWave::neumann )
                                                {
                          // include boundary
                                                }
                                                else if( bc>0 )
                                                {
                                                    printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                    OV_ABORT("error");
                                                }
                                                else if( bc<0 )
                                                {
                          // periodic -- include left end
                                                    if( side==1 )
                                                        iab[side] += is; 
                                                }
                                                else
                                                {
                          // interpolation boundary : include end 
                                                }
                                            }
                                            Iv[axis] = Range(iab[0],iab[1]);
                                        }
                                    bool ok=ParallelUtility::getLocalArrayBounds(Q[grid],uuLocal,I1,I2,I3);
                                    if( ok )
                                    {
                                        FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                        {
                                            if( maskLocal(ii1,ii2,ii3)>0 )
                                                dotProduct += uuLocal(ii1,ii2,ii3,ia)*vvLocal(ii1,ii2,ii3,ja);
                                        }
                                    }
                                }
                                dotProduct = ParallelUtility::getSum( dotProduct );
                            }  
            // printF("Rayleigh-Ritz ( Q[%d],Q[%d] ) = %9.3e\n",ia,ja,dotProduct);
                        if( ia==ja ) 
                            dotProduct -= 1.;  // diagonal entries should be 1
                        maxErr=max(maxErr,fabs(dotProduct));
                    }
                }
                printF("*** MAX ERR in ORTHO-NORMALIZATION ( Qi,Qj )=delta_ij : %8.2e\n",maxErr);
                if( false && numVectorsSaved==3 )
                {
                    OV_ABORT("stop here for now");
                }
            }
      // ----- END Gram-Schmidt orthogonalization -------



    // ----- END Gram-Schmidt orthogonalization -------

    // Form entries in Bm (assume symmetric)
    //    B[i,j] = qi * L qj 

        for( int ia=0; ia<numberOfVectorsToUse; ia++ )
        {
      // Compute   lap = L q[ia]
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                OV_GET_SERIAL_ARRAY(Real,Q[grid],qLocal);
                OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);
                getIndex(cg[grid].dimension(),I1,I2,I3);
                operators[grid].derivative(MappedGridOperators::laplacianOperator,qLocal,lapLocal,I1,I2,I3,ia);
            }

      // Compute ( q[ja], Lap q[ia] )

      // *** NEED TO SET BOUNDARY POINTS TO ZERO ****
            w=0.; 

      // for( int ja=0; ja<=ia; ja++ )
            for( int ja=0; ja<numberOfVectorsToUse; ja++ )
            {
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];

          // const int extra=-1;
                    getIndex(mg.gridIndexRange(),I1,I2,I3);   

                    OV_GET_SERIAL_ARRAY(Real,Q[grid],qLocal);
                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);      
                    OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);      

                    wLocal(I1,I2,I3) = qLocal(I1,I2,I3,ja)*lapLocal(I1,I2,I3);
                }

                w.periodicUpdate();
                w.interpolate();  

                Bm(ia,ja) = integrate.volumeIntegral(w);  
        // if( ia!=ja )
        //   Bm(ja,ia) = Bm(ia,ja);  // symmetric entry
            }
        }


  
        int info, n=numberOfVectorsToUse, lda=n; 
    // int lwork = 3*n-1;
        int lwork = 5*n;   // at least 4*n for GEEV
        RealArray A(n,n), eigs(n), work(lwork);
  
        if( true )
        {
      // do NOT assume matrix is symmetric
  
      // jobvl : "N" : left eigenvectors are not computed
      // jobvr : "V" " right eigenvectors are computed
            RealArray wr(n), wi(n); // eigenvalue 
            RealArray vl(1,1), vr(n,n); 
            int ldvl=1, ldvr=n;
            Range N=n; 
            A(N,N)=Bm(N,N);      
            GEEV( "N", "V", n, A(0,0), lda, wr(0), wi(0), vl(0,0), ldvl, vr(0,0), ldvr, work(0), lwork, info );

            if( info!=0 )
            {
                printF("ERROR return from symmetric eigenvalue routine GEEV: info=%d\n",info);
                OV_ABORT("ERROR");
            }      

      // void GEEV( char *jobvl, char* jobvr, int & n, real & a, const int & lda,
      //         real & wr, real & wi, real &vl, int & ldvl, real & vr, int & ldvr, real & work, int & lwork, int & info );

            for( int i=0; i<n; i++ )
            {
                eigs(i)=wr(i);
                printF("GEEV: lambda(%d) = %12.4e + I* %12.4e\n",i,wr(i),wi(i)); 
            }
            A = vr;  // eigenvectors
        }
        else
        {
      // --- compute eigenvalues assuming a symmetric matrix ----

      // Compute the eigenvalues of Bm(0:ia,0:ia)
      //   jobz : "N" compute eigs, "V" compute eigs and EV's
      //   uplo : "U" or "L" matrix is provided in upper or lower part
      //   eigs : eigenvalues inb ascending order
      //   A  : matrix on input, eigenvectors on output

            Range N=n; 
            A(N,N)=Bm(N,N);
      // Real lamrq = sqrt(-A(n-1,n-1));
      // printF("**** B(n-1,n-1) => RQ = %12.4e\n",lamrq);
            SYEV( "V", "U", n, A(0,0), lda, eigs(0), work(0), lwork, info );
            if( info!=0 )
            {
                printF("ERROR return from symmetric eigenvalue routine SYEV: info=%d\n",info);
                OV_ABORT("ERROR");
            }

        }

    // --- FIND closest Rayleigh Ritz eigenvalue to the Rayeligh Quotient Eigenvalue ---
        Real lambdaRQ=eigenValues(0); 
        Real lambdaRR; // holds closest Ritz eigenvalue to the RQ eigenvalue
        int indexRR=-1;
    // int indexRR=0;
    // Real lambdaRR=eigs(indexRR);
    // if( lambdaRR<0 )
    //   lambdaRR = sqrt(-lambdaRR);
    // else
    // {
    //   printF(" Rayleigh-Ritz: eigs[%d]=%12.4e **WARNING: WRONG SIGN****\n",0,eigs(0));
    // }
    // printF(" Rayleigh-Ritz: lambda[%d]=%16.8e (lambdaRQ=%16.8e)\n",0,lambdaRR,lambdaRQ);
    // Real relErrLam = fabs(lambda-lambdaRQ)/lambdaRQ; 
        Real relErrLam = REAL_MAX;
        int numToKeep=1;
        for( int i=0; i<n; i++ )
        {
            Real lambda = eigs(i);
            if( lambda<=0 )
            {
                numToKeep=i+1; 
                lambda=sqrt(-lambda);
                printF(" Rayleigh-Ritz: lambda[%d]=%16.8e (lambdaRQ=%16.8e)\n",i,lambda,lambdaRQ);

                Real err = fabs(lambda-lambdaRQ)/lambdaRQ; 
                if( err < relErrLam )
                {
          // printF("closer: fabs(lambda-lambdaRQ)=%9.3e < fabs(lambda-lambdaRR)=%9.3e")
                    relErrLam = err;
                    lambdaRR=lambda;
                    indexRR=i;
                }
            }
            else
            {
                printF(" Rayleigh-Ritz: eigs[%d]=%12.4e **WARNING: WRONG SIGN****\n",i,eigs(i));

            }
        }
        assert( indexRR>=0 );
        printF(" Rayleigh-Ritz: iterationRR=%d: closest: eigs[%d]= lambdaRR=%16.8e (lambdaRQ=%16.8e)\n",iterationRR,indexRR,lambdaRR,lambdaRQ);

        eigsRayleighRitz(iterationRR)=lambdaRR;

    // -- compute the Rayleigh-Ritz eigenvector 
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            getIndex(mg.dimension(),I1,I2,I3);  

            OV_GET_SERIAL_ARRAY(Real,Q[grid],qLocal);
            OV_GET_SERIAL_ARRAY(Real,eigenVectorRayleighRitz[grid],eigenVectorRayleighRitzLocal);      

            eigenVectorRayleighRitzLocal(I1,I2,I3)=0.;
            if( grid==0 )
            {
                printF("RR: eigenvector[%d]=[",indexRR);
                for( int i=0; i<numToKeep; i++ )
                    printF("%12.4e, ",A(i,indexRR));
                printF("]\n");
            }
            for( int i=0; i<numToKeep; i++ )
            {
                eigenVectorRayleighRitzLocal(I1,I2,I3) += A(i,indexRR)*qLocal(I1,I2,I3,i);  // assumes eigenvector stored in a column
            }
            

        }

    // -- normalize by max-norm
        Real evNorm = maxNorm(eigenVectorRayleighRitz);
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,eigenVectorRayleighRitz[grid],eigenVectorRayleighRitzLocal);  
            eigenVectorRayleighRitzLocal *= (1./evNorm);
        }

        Real relErrEigenvalue, relErrEigenvector;
        bool checkRayleighRitz=true;
        getErrorsInEigenmodes( relErrEigenvalue, relErrEigenvector, checkRayleighRitz );
        printF("\n ++++ Rayleigh-Ritz: rel-err in eig=%8.2e, eigenvector=%8.2e (numToKeep=%d, numVectorsSaved=%d)\n\n",relErrEigenvalue, relErrEigenvector,numToKeep,numVectorsSaved);
    // OV_ABORT("FINISH ME : RQ");


    // ** CHECK eigenPair residual and save Ritz pair if small ********* FINISH ME *********

        iterationRR++;




    } // end Rayleigh-Quotient


    return 0;
}


// ===============================================================
/// \brief Correct the current guess at the eigenfunction
///  using the current Ritz vector.
// ===============================================================
int CgWave::
correctEigenfunction()
{

    const int & computeEigenmodes      = dbase.get<int>("computeEigenmodes");
    if( !computeEigenmodes )
        return 0;

    RealArray & resVector              = dbase.get<RealArray>("resVector");  // WaveHoltz residual
    const int & assignRitzFrequency    = dbase.get<int>("assignRitzFrequency");
    const int & iteration              = dbase.get<int>("iteration");
    const int & iterationRR            = dbase.get<int>("iterationRR");
    const int & iterationStartRR       = dbase.get<int>("iterationStartRR");      // start at this WH iteration
    realCompositeGridFunction & v      = dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & eigenVectorRayleighRitz = dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");

    int it = iteration-iterationStartRR; 
    if( it>0 && ( it % assignRitzFrequency == 0 ) )
    {
        printF("\n ZZZZZZZZZ  CORRECT EIGENFUNCTION : ASSIGN THE CURRENT RITZ VECTOR at iteration=%d, iterationStartRR=%d ZZZZZZZZZZZ\n",
                    iteration,iterationStartRR);
        const int freq=0;
        Index I1,I2,I3;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,eigenVectorRayleighRitz[grid],eigenVectorRayleighRitzLocal); 
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            getIndex(cg[grid].dimension(),I1,I2,I3);

            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
            if( ok )
                vLocal(I1,I2,I3,freq) = eigenVectorRayleighRitzLocal(I1,I2,I3);
        }      
    }

    return 0;
}


// ===============================================================
/// \brief Adjust the driving frequency for the EigenWave algorithm
// ===============================================================
int CgWave::
adjustEigenWaveFrequency()
{
    bool adjustEigenFrequency=true;

    const int & computeEigenmodes = dbase.get<int>("computeEigenmodes");
    if( !adjustEigenFrequency || !computeEigenmodes )
        return 0;

    const int & iteration = dbase.get<int>("iteration");

    RealArray & resVector              = dbase.get<RealArray>("resVector");  // WaveHoltz residual
    const Real residual = iteration>0 ? resVector(iteration-1) : 1.;

  // ---- Only start adjusting omega when residual goes below a given tolerance
  // const Real residTolRR=1.e-2;
    const Real residTolRR=1.e-2;
    if( (residual > residTolRR) || iteration<5  )
        return 0;


  // if( iteration<5 )  // ************ DO THIS FOR NOW ****************** FIX ME 
  //   return 0;

    RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
    RealArray & periodArray           = dbase.get<RealArray>("periodArray");

    IntegerArray & numPeriodsArray    = dbase.get<IntegerArray>("numPeriodsArray"); 

    RealArray & eigsRayleighQuotient = dbase.get<RealArray>("eigsRayleighQuotient");
    const Real omegaNew = eigsRayleighQuotient(iteration);

    frequencyArray(0) = omegaNew;
    periodArray(0)    = (twoPi/frequencyArray(0))*numPeriodsArray(0);

  // do this too: 
    real & omega                      = dbase.get<real>("omega");
    real & Tperiod                    = dbase.get<real>("Tperiod");
    real & tFinal                     = dbase.get<real>("tFinal");

    omega   = frequencyArray(0);
    Tperiod = periodArray(0);
    tFinal  = periodArray(0);    // This is important!

  // RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");
  // RealArray & periodArraySave       = dbase.get<RealArray>("periodArraySave");

  // frequencyArraySave(0)  = frequencyArray(0);
  // periodArraySave(0)     = periodArray(0);

  // frequencyArraySave(0) = omegaNew;
  // periodArraySave(0)    = (twoPi/frequencyArraySave(0))*numPeriodsArray(0);

  // frequencyArray(0)  = frequencyArraySave(0);
  // periodArray(0)     = periodArraySave(0);

  // // do this too??
  // real & omega                      = dbase.get<real>("omega");
  // real & Tperiod                    = dbase.get<real>("Tperiod");
  // real & omegaSave                  = dbase.get<real>("omegaSave");
  // omega     = omegaNew;
  // omegaSave = omegaNew;
  // Tperiod   = periodArraySave(0);

  // initialize();

    printF("\n ***** ADJUSTING FREQUENCY FOR EIGNWAVE: iteration=%d, omega=%14.6e, T=%14.6e (from RQ)\n\n",iteration,omegaNew,Tperiod);

    return 0;
}

// ===============================================================
// Macro: Add to grid functions
//     w = fact1*v + fact2*uev   
// ===============================================================

// ==========================================================================
/// \brief Estimate eigenvalues and eigenvectors
//    ***OLD VERSION ***  
/// \param relErrEigenvalue, relErrEigenvector (output) : relative errors
/// \param checkRayleighRitz (input) : if true, check the error in the Rayleigh-Ritz approximation
/// \Note: The error in the eigenvector is stored in the "error" grid function.
/// \Return value: 1 = success, 0=errors not computed 
///
/// \Details: The computed eigenvectors to check are: 
///  if( checkRayleighRitz )
///      dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz")
///   else
///     dbase.get<realCompositeGridFunction>("v")
// ===========================================================================
int CgWave::getErrorsInEigenmodes( Real & relErrEigenvalue, Real & relErrEigenvector, bool checkRayleighRitz /* =false */ )
{
  //    ***OLD VERSION ***  
  // if( true )
  //   OV_ABORT("error -- obsolete");

    const int & computeEigenmodes = dbase.get<int>("computeEigenmodes");
    if( !computeEigenmodes ) 
        return 0;

  // -- This next call will read in any known eigenmodes ---
  // cgWave.initializeDeflation();
    if( !dbase.has_key("uev") )  
        return 0;


    const int & iteration = dbase.get<int>("iteration");
    const int & useAccurateInnerProduct = dbase.get<int>("useAccurateInnerProduct");

  // ----- determine errors in eigenmodes if we know the true eigenmodes -----
    printF(" CgWave::computeErrorsInEigenmodes: COMPUTE ERRORS IN EIGENMODES iteration=%d (checkRayleighRitz=%d)------\n",
              iteration,(int)checkRayleighRitz);


    const int & numberOfFrequencies   = dbase.get<int>("numberOfFrequencies");
    RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");

  // Computed eigenvalues are stored here: 
    RealArray & eigenValues = dbase.get<RealArray>("eigenValues");

    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
    CompositeGrid & cg = *v.getCompositeGrid();

  // ----------------------------------------------
  // -------- eigenvector(s) to check: ------------
  // ----------------------------------------------
    realCompositeGridFunction & eigVector = checkRayleighRitz==0 ? v : dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");

  // Here are the "true" eigenvectors and eigenvalues:
    realCompositeGridFunction & uev = dbase.get<realCompositeGridFunction>("uev");
    RealArray & eig                 = dbase.get<RealArray>("eig");
    IntegerArray & eigMultiplicity  = dbase.get<IntegerArray>("eigMultiplicity");

    const int numberOfEigenvectors = uev.getComponentBound(0) - uev.getComponentBase(0) + 1;
    printF(">> There are %d eigenvectors (true, from file).\n",numberOfEigenvectors);

    const int maxIterations=500; // *** FIX ME ***
    if( !dbase.has_key("errEigenvector") )
    {
        RealArray & errEigenvalue = dbase.put<RealArray>("errEigenvalue");
        errEigenvalue.redim(maxIterations);
        errEigenvalue=1.;

        RealArray & errEigenvector = dbase.put<RealArray>("errEigenvector");
        errEigenvector.redim(maxIterations);
        errEigenvector=1.;

        dbase.put<int>("eigIndex");
    }
    RealArray & errEigenvalue  = dbase.get<RealArray>("errEigenvalue");
    RealArray & errEigenvector = dbase.get<RealArray>("errEigenvector");
    int & eigIndex             = dbase.get<int>("eigIndex");


    if( checkRayleighRitz && !dbase.has_key("errEigenvectorRR") )
    {
        RealArray & errEigenvalueRR = dbase.put<RealArray>("errEigenvalueRR");
        errEigenvalueRR.redim(maxIterations);
        errEigenvalueRR=1.;

        RealArray & errEigenvectorRR = dbase.put<RealArray>("errEigenvectorRR");
        errEigenvectorRR.redim(maxIterations);
        errEigenvectorRR=1.;
    }

    Range all;
    uev.updateToMatchGrid(cg,all,all,all,numberOfEigenvectors); // the cg in uev may have gone out of scope --- FIX ME ??

    if( !dbase.has_key("integrate") )
    {
        Integrate & integrate = dbase.put<Integrate>("integrate");
        integrate.updateToMatchGrid(cg);
    }
    Integrate & integrate = dbase.get<Integrate>("integrate");

    realCompositeGridFunction & error = dbase.get<realCompositeGridFunction>("error");
    error.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
    error=0.;
    for( int freq=0; freq<numberOfFrequencies; freq++ )
        error.setName(sPrintF("err%d",freq),freq);

    RealCompositeGridFunction w(cg);  // temp space 
    Index I1,I2,I3;
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
        Real lambda = eigenValues(freq);  // computed eigenvalue
        if( checkRayleighRitz )
        {
            RealArray & eigsRayleighRitz = dbase.get<RealArray>("eigsRayleighRitz");
            int & iterationRR            = dbase.get<int>("iterationRR");
            int ilamRR = min(iterationRR-1,eigsRayleighRitz.getBound(0));
            lambda = eigsRayleighRitz(ilamRR);
            if( lambda<0 )
            {
                printF("ERROR: lambda(Ritz)=%9.2e, ilamRR=%d, iterationRR=%d\n",lambda,ilamRR,iterationRR);
                ::display(eigsRayleighRitz,"eigsRayleighRitz");
                OV_ABORT("error");
            }
        }

        Real diffMin = REAL_MAX*.1;
    // int eigIndex=0;  
        for( int i=0; i<numberOfEigenvectors; i++ )
        {
      // Assume eigenvalues are sorted so we keep looking while the difference decreases
            Real diff = fabs( lambda - eig(0,i) ) ;
            if( diff < diffMin )
            {
                  eigIndex=i; diffMin=diff;
            }
            else
            {
                break;
            }
        }
        relErrEigenvalue = fabs(lambda-eig(0,eigIndex))/max(1.,fabs(eig(0,eigIndex)));


    // --- Compute the error in the eigenvector ---
    // const Real evNorm        = maxNorm( uev,eigIndex );    // may not have max-norm 1
    // if( fabs(evNorm-1.) > REAL_EPSILON*1000. )
    // {
    //   // normalize exact eigenvalue to match the computed : maxNorm =1
    //   printF("NORMALIZE EXACT EIGENVECTOR (norm was %9.2e)\n",evNorm); 
    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     uev[grid](all,all,all,eigIndex) *= (1./evNorm); 
    //   }
    // }

        const Real eigVectorNorm = maxNorm( eigVector,freq );  // should be 1   
        if( fabs(eigVectorNorm-1.) > REAL_EPSILON*1000. )
        {
            printF("INFO: computed eigVector does not have max-norm 1! eigVectorNorm=%12.5e (eigVectorNorm-1=%9.3e). Scaling...\n",eigVectorNorm,eigVectorNorm-1.);
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(Real,eigVector[grid],eigVectorLocal);

                bool ok=ParallelUtility::getLocalArrayBounds(eigVector[grid],eigVectorLocal,I1,I2,I3);
                if( ok )
                    eigVectorLocal(I1,I2,I3,freq) *= (1./eigVectorNorm); 

        // eigVector[grid](all,all,all,freq) *= (1./eigVectorNorm); 
            }
      // OV_ABORT("error");
        }

        relErrEigenvector = REAL_MAX*.1;
        printF("***INFO: eigenvalue eigIndex=%d has multiplicity=%d (norm=%9.3e)****\n",eigIndex,eigMultiplicity(eigIndex),eigVectorNorm);
        if( eigMultiplicity(eigIndex)==1 )
        {
      // ---- simple eigenvalue --
            if( true )
            {
        //  Project eigenVector v onto uev 
        //  v = alpha*uev
        //  alpha = (v,uev)/(uev,euv)
                Real dotProduct;
                    if( useAccurateInnerProduct )
                    {
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,freq)*vvLocal(I1,I2,I3,eigIndex);
                        }
                        dotProduct = integrate.volumeIntegral(w);
                    }
                    else
                    {
                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                        int iab[2]; 
                        int ii1,ii2,ii3;  
                        dotProduct=0.; 
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            const IntegerArray & gid = mg.gridIndexRange();
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                Iv[2]=Range(0,0);
                                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                {
                                    for( int side=0; side<=1; side++ )
                                    {
                                        int is = 1-2*side;
                                        iab[side]=gid(side,axis);
                                        const int bc = mg.boundaryCondition(side,axis);
                                        if( bc==CgWave::dirichlet )
                                        {
                                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                        }
                                        else if( bc==CgWave::neumann )
                                        {
                      // include boundary
                                        }
                                        else if( bc>0 )
                                        {
                                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                            OV_ABORT("error");
                                        }
                                        else if( bc<0 )
                                        {
                      // periodic -- include left end
                                            if( side==1 )
                                                iab[side] += is; 
                                        }
                                        else
                                        {
                      // interpolation boundary : include end 
                                        }
                                    }
                                    Iv[axis] = Range(iab[0],iab[1]);
                                }
                            bool ok=ParallelUtility::getLocalArrayBounds(eigVector[grid],uuLocal,I1,I2,I3);
                            if( ok )
                            {
                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                {
                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                        dotProduct += uuLocal(ii1,ii2,ii3,freq)*vvLocal(ii1,ii2,ii3,eigIndex);
                                }
                            }
                        }
                        dotProduct = ParallelUtility::getSum( dotProduct );
                    }  

                Real uevDot;
                    if( useAccurateInnerProduct )
                    {
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,eigIndex)*vvLocal(I1,I2,I3,eigIndex);
                        }
                        uevDot = integrate.volumeIntegral(w);
                    }
                    else
                    {
                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                        int iab[2]; 
                        int ii1,ii2,ii3;  
                        uevDot=0.; 
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            const IntegerArray & gid = mg.gridIndexRange();
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                Iv[2]=Range(0,0);
                                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                {
                                    for( int side=0; side<=1; side++ )
                                    {
                                        int is = 1-2*side;
                                        iab[side]=gid(side,axis);
                                        const int bc = mg.boundaryCondition(side,axis);
                                        if( bc==CgWave::dirichlet )
                                        {
                                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                        }
                                        else if( bc==CgWave::neumann )
                                        {
                      // include boundary
                                        }
                                        else if( bc>0 )
                                        {
                                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                            OV_ABORT("error");
                                        }
                                        else if( bc<0 )
                                        {
                      // periodic -- include left end
                                            if( side==1 )
                                                iab[side] += is; 
                                        }
                                        else
                                        {
                      // interpolation boundary : include end 
                                        }
                                    }
                                    Iv[axis] = Range(iab[0],iab[1]);
                                }
                            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                            if( ok )
                            {
                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                {
                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                        uevDot += uuLocal(ii1,ii2,ii3,eigIndex)*vvLocal(ii1,ii2,ii3,eigIndex);
                                }
                            }
                        }
                        uevDot = ParallelUtility::getSum( uevDot );
                    }  

                Real alpha =  dotProduct/uevDot;
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    error[grid] = alpha*uev[grid](all,all,all,eigIndex) - eigVector[grid](all,all,all,freq);  // holds error
                }                  

                relErrEigenvector = maxNorm( error )/eigVectorNorm;

            }
            else
            {
                Real dotProduct;
                    if( useAccurateInnerProduct )
                    {
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,freq)*vvLocal(I1,I2,I3,eigIndex);
                        }
                        dotProduct = integrate.volumeIntegral(w);
                    }
                    else
                    {
                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                        int iab[2]; 
                        int ii1,ii2,ii3;  
                        dotProduct=0.; 
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            const IntegerArray & gid = mg.gridIndexRange();
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                Iv[2]=Range(0,0);
                                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                {
                                    for( int side=0; side<=1; side++ )
                                    {
                                        int is = 1-2*side;
                                        iab[side]=gid(side,axis);
                                        const int bc = mg.boundaryCondition(side,axis);
                                        if( bc==CgWave::dirichlet )
                                        {
                                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                        }
                                        else if( bc==CgWave::neumann )
                                        {
                      // include boundary
                                        }
                                        else if( bc>0 )
                                        {
                                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                            OV_ABORT("error");
                                        }
                                        else if( bc<0 )
                                        {
                      // periodic -- include left end
                                            if( side==1 )
                                                iab[side] += is; 
                                        }
                                        else
                                        {
                      // interpolation boundary : include end 
                                        }
                                    }
                                    Iv[axis] = Range(iab[0],iab[1]);
                                }
                            bool ok=ParallelUtility::getLocalArrayBounds(eigVector[grid],uuLocal,I1,I2,I3);
                            if( ok )
                            {
                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                {
                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                        dotProduct += uuLocal(ii1,ii2,ii3,freq)*vvLocal(ii1,ii2,ii3,eigIndex);
                                }
                            }
                        }
                        dotProduct = ParallelUtility::getSum( dotProduct );
                    }  
                const Real eSign = dotProduct>0 ? -1. : +1.;
        // w = 1.*v + eSign*uev 
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        getIndex(cg[grid].dimension(),I1,I2,I3);
                        OV_GET_SERIAL_ARRAY(Real,eigVector[grid],vvLocal);
                        OV_GET_SERIAL_ARRAY(Real,uev[grid],uevLocal);
                        OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                        bool ok=ParallelUtility::getLocalArrayBounds(eigVector[grid],vvLocal,I1,I2,I3);
                        if( ok )
                            wLocal(I1,I2,I3) = 1.*vvLocal(I1,I2,I3,freq) + eSign*uevLocal(I1,I2,I3,eigIndex);
                    }

                relErrEigenvector = maxNorm( w )/eigVectorNorm;
            }

      // if( 1==1 )
      // {
      //   Real uevMaxNorm = maxNorm( uev,eigIndex );
      //   printF("TRUE: maxNorm=%9.3e\n",uevMaxNorm);
      // }

        }
        else if( eigMultiplicity(eigIndex)==2 )
        {
      // Project the eigenvector onto the eigenspace spanned by 
      //     ue1v = uev[eigIndex] 
      //  and 
      //     ue2v = uev[eigIndex+1]
      //
      //    v = alpha1*uev1 + alpha2*euv2
      //    A = [ uev1 uev2 ] 
      //      A [ alphav] = b 
      // Least Squares:
      //    [ A11 A12 ] [alpha1] = [ b1 ]

            int eigIndex1 = eigIndex;
            int eigIndex2 = eigIndex+1;

      // const Real eigTol=1.e-4;
            const Real eigTol = dbase.get<Real>("eigenValueTolForMultiplicity");

            if( fabs( eig(0,eigIndex1)-eig(0,eigIndex2) )/fabs( eig(0,eigIndex1) ) > eigTol )
            {
        // no match to eigIndex+1, try eigIndex-1 (maybe just find closest)
                eigIndex2=max(0,eigIndex-1); // try this one
                if( fabs( eig(0,eigIndex1)-eig(0,eigIndex2) )/fabs( eig(0,eigIndex1) ) > eigTol )
                {              
                    printF("ERROR: multiple eigenvalue but eigIndex1=%d is not the same a eigIndex2=%d\n",eigIndex1,eigIndex2);
                    printF("  eig(0,eigIndex1)=%12.4e, eig(0,eigIndex2)=%12.4e\n",eig(0,eigIndex1),eig(0,eigIndex2));
                    OV_ABORT("error");
                }
            }
            if( eigIndex1>eigIndex2 )
            { // flip these
                int temp=eigIndex1; eigIndex1=eigIndex2; eigIndex2=temp;
            }

            if( eigMultiplicity(eigIndex2)!=2  )
            {
                printF("\nERROR: multiple eigenvalue but eigIndex1=%d mult=%d, and eigIndex2=%d has mult=%d\n\n",
                    eigIndex1,eigMultiplicity(eigIndex1),eigIndex2,eigMultiplicity(eigIndex2));
                printF(" eig(0,eigIndex1  )=%16.8e\n",eig(0,eigIndex1));
                printF(" eig(0,eigIndex2  )=%16.8e\n",eig(0,eigIndex2));
                printF(" eig(0,eigIndex1-1)=%16.8e\n",eig(0,eigIndex1-1));
                printF(" eig(0,eigIndex1+1)=%16.8e\n",eig(0,eigIndex1+1));
            }


            
            RealArray a(2,2), ai(2,2), b(2);
            for( int i1=0; i1<2; i1++ )
            {
                int ii = i1==0 ? eigIndex1 : eigIndex2;
                    if( useAccurateInnerProduct )
                    {
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,freq);
                        }
                        b(i1) = integrate.volumeIntegral(w);
                    }
                    else
                    {
                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                        int iab[2]; 
                        int ii1,ii2,ii3;  
                        b(i1)=0.; 
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            const IntegerArray & gid = mg.gridIndexRange();
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                Iv[2]=Range(0,0);
                                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                {
                                    for( int side=0; side<=1; side++ )
                                    {
                                        int is = 1-2*side;
                                        iab[side]=gid(side,axis);
                                        const int bc = mg.boundaryCondition(side,axis);
                                        if( bc==CgWave::dirichlet )
                                        {
                                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                        }
                                        else if( bc==CgWave::neumann )
                                        {
                      // include boundary
                                        }
                                        else if( bc>0 )
                                        {
                                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                            OV_ABORT("error");
                                        }
                                        else if( bc<0 )
                                        {
                      // periodic -- include left end
                                            if( side==1 )
                                                iab[side] += is; 
                                        }
                                        else
                                        {
                      // interpolation boundary : include end 
                                        }
                                    }
                                    Iv[axis] = Range(iab[0],iab[1]);
                                }
                            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                            if( ok )
                            {
                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                {
                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                        b(i1) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,freq);
                                }
                            }
                        }
                        b(i1) = ParallelUtility::getSum( b(i1) );
                    }  

                for( int i2=0; i2<2; i2++ )
                {
                    int jj = i2==0 ? eigIndex1 : eigIndex2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,jj);
                            }
                            a(i1,i2) = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            a(i1,i2)=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            a(i1,i2) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,jj);
                                    }
                                }
                            }
                            a(i1,i2) = ParallelUtility::getSum( a(i1,i2) );
                        }  
                }

            }
      //  A  = ( a11 a12 )
      //       ( a21 a22 )
      //  A^{-1} = (1/det)*( a22  -a12 )
      //                   ( -a21  a11)
            Real det = a(0,0)*a(1,1) - a(0,1)*a(1,0);
            ai(0,0)=  a(1,1)/det;
            ai(0,1)= -a(0,1)/det;
            ai(1,0)= -a(1,0)/det;
            ai(1,1)=  a(0,0)/det;
            Real alpha1 = ai(0,0)*b(0) + ai(0,1)*b(1);
            Real alpha2 = ai(1,0)*b(0) + ai(1,1)*b(1);
            printF("\n -- double eig found: eigIndex1=%d, eigIndex2=%d, alpha1=%g, alpha2=%g (det=%9.2e)\n",eigIndex1,eigIndex2,alpha1,alpha2,det);

            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(Real,error[grid],errorLocal);
                OV_GET_SERIAL_ARRAY(Real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(Real,eigVector[grid],eigVectorLocal);

                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                    errorLocal = alpha1*uevLocal(I1,I2,I3,eigIndex1) + alpha2*uevLocal(I1,I2,I3,eigIndex2) - eigVectorLocal(I1,I2,I3,freq);  // holds error

        // error[grid] = alpha1*uev[grid](all,all,all,eigIndex1) + alpha2*uev[grid](all,all,all,eigIndex2) - eigVector[grid](all,all,all,freq);  // holds error
            }                  

            relErrEigenvector = maxNorm( error )/eigVectorNorm;


        }
        else
        {
      // **New way -- arbitrary multiplicity ---

      // -- this came from genEigs.bC ----
            const int multiplicity=eigMultiplicity(eigIndex);
            int md=multiplicity;

            const int i0=eigIndex; 
            RealArray a(md,md), ai(md,md), b(md), alpha(md);
            for( int i1=0; i1<md; i1++ )
            {
                int ii = i0 + i1;
                    if( useAccurateInnerProduct )
                    {
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,freq);
                        }
                        b(i1) = integrate.volumeIntegral(w);
                    }
                    else
                    {
                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                        int iab[2]; 
                        int ii1,ii2,ii3;  
                        b(i1)=0.; 
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            const IntegerArray & gid = mg.gridIndexRange();
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                Iv[2]=Range(0,0);
                                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                {
                                    for( int side=0; side<=1; side++ )
                                    {
                                        int is = 1-2*side;
                                        iab[side]=gid(side,axis);
                                        const int bc = mg.boundaryCondition(side,axis);
                                        if( bc==CgWave::dirichlet )
                                        {
                                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                        }
                                        else if( bc==CgWave::neumann )
                                        {
                      // include boundary
                                        }
                                        else if( bc>0 )
                                        {
                                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                            OV_ABORT("error");
                                        }
                                        else if( bc<0 )
                                        {
                      // periodic -- include left end
                                            if( side==1 )
                                                iab[side] += is; 
                                        }
                                        else
                                        {
                      // interpolation boundary : include end 
                                        }
                                    }
                                    Iv[axis] = Range(iab[0],iab[1]);
                                }
                            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                            if( ok )
                            {
                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                {
                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                        b(i1) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,freq);
                                }
                            }
                        }
                        b(i1) = ParallelUtility::getSum( b(i1) );
                    }  

                for( int i2=0; i2<md; i2++ )
                {
                    int jj = i0 + i2;
                        if( useAccurateInnerProduct )
                        {
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,jj);
                            }
                            a(i1,i2) = integrate.volumeIntegral(w);
                        }
                        else
                        {
                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                            int iab[2]; 
                            int ii1,ii2,ii3;  
                            a(i1,i2)=0.; 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                MappedGrid & mg = cg[grid];
                                const IntegerArray & gid = mg.gridIndexRange();
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                                OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                    Iv[2]=Range(0,0);
                                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                    {
                                        for( int side=0; side<=1; side++ )
                                        {
                                            int is = 1-2*side;
                                            iab[side]=gid(side,axis);
                                            const int bc = mg.boundaryCondition(side,axis);
                                            if( bc==CgWave::dirichlet )
                                            {
                                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                            }
                                            else if( bc==CgWave::neumann )
                                            {
                        // include boundary
                                            }
                                            else if( bc>0 )
                                            {
                                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                                OV_ABORT("error");
                                            }
                                            else if( bc<0 )
                                            {
                        // periodic -- include left end
                                                if( side==1 )
                                                    iab[side] += is; 
                                            }
                                            else
                                            {
                        // interpolation boundary : include end 
                                            }
                                        }
                                        Iv[axis] = Range(iab[0],iab[1]);
                                    }
                                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                                if( ok )
                                {
                                    FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                    {
                                        if( maskLocal(ii1,ii2,ii3)>0 )
                                            a(i1,i2) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,jj);
                                    }
                                }
                            }
                            a(i1,i2) = ParallelUtility::getSum( a(i1,i2) );
                        }  
                }

            }

      // PA = LU factor
            IntegerArray ipvt(md);
            int info;
            GETRF( md,md,a(0,0), md, ipvt(0), info );                
            if( info!=0 )
            {
                printF("ERROR return from GETRF, info=%d\n",info);
                OV_ABORT("error");
            }
      // Solve given LU
            int nrhs=1;
            GETRS( "N", md, nrhs, a(0,0), md, ipvt(0), b(0), md, info );
            alpha=b;
      // ::display(alpha,"alpha");
            if( info!=0 )
            {
                printF("ERROR return from GETRS, info=%d\n",info);
                OV_ABORT("error");
            }                


            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
        // WARNING: We cannot scale here since we may make both eigenvectors the same!
        // u[grid](all,all,all,i) = alpha1*u[grid](all,all,all,i) + alpha2*u[grid](all,all,all,j);
        // general case 
                OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                OV_GET_SERIAL_ARRAY(Real,eigVector[grid],eigVectorLocal);
                OV_GET_SERIAL_ARRAY(Real,uev[grid],uevLocal);

                getIndex(cg[grid].dimension(),I1,I2,I3);
                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                {
                    wLocal(I1,I2,I3) = -eigVectorLocal(I1,I2,I3,freq);
                    for( int i=0; i<multiplicity; i++ )
                        wLocal(I1,I2,I3) += alpha(i)*uevLocal(I1,I2,I3,i0+i);   
                }       


            }       

        }

    // // save errors for plotting 
    // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    // {
    //   error[grid]=w[grid];
    // }

        printF("freq=%d: lambda=%18.10e, true=%18.10e, eigenValue-rel-err=%8.2e, eigenVector-rel-err=%8.2e (checkRayleighRitz=%d)\n",
                      freq,lambda,eig(0,eigIndex),relErrEigenvalue,relErrEigenvector,(int)checkRayleighRitz);

        if( !checkRayleighRitz )
        {
            if( iteration<maxIterations )
            {
                errEigenvalue(iteration)  = relErrEigenvalue;
                errEigenvector(iteration) = relErrEigenvector;
            }
        }
        else
        {
            RealArray & errEigenvalueRR  = dbase.get<RealArray>("errEigenvalueRR");
            RealArray & errEigenvectorRR = dbase.get<RealArray>("errEigenvectorRR");
            int & iterationRR            = dbase.get<int>("iterationRR");
            errEigenvalueRR(iterationRR)  = relErrEigenvalue;
            errEigenvectorRR(iterationRR) = relErrEigenvector;

        }

    }

    printF("CgWave::computeErrorsInEigenmodes: iteration=%d: relative errors in eigenvalue=%8.2e, eigenvector=%8.2e\n",iteration,relErrEigenvalue, relErrEigenvector);

    return 1;
}


// ==========================================================================
/// \brief Compute the error in an eignavlue, eigenvector pair
/// \param lambda, eigVector, component (input): eigenvalue and eigenvector to check.
/// \param lambdaTrue  (output) : true eigenvalue
/// \param relErrEigenvalue, relErrEigenvector  (output) : relative errors
/// \param eigIndex (output) : index into array of true eigenPairs of the closest eigenpair
/// \param multipleEigIndex (output) : for multiple eigs, index of first eigenpair.
/// \param saveErrors (input) : if true, save the error in the grid function "error"
///
/// \Return value: 1 = success, 0=errors not computed 
///
// ===========================================================================
int CgWave::getErrorInEigenPair( const Real lambda, realCompositeGridFunction & eigVector, const int component, 
                                                                  Real & lambdaTrue, Real & relErrEigenvalue, Real & relErrEigenvector, 
                                                                  int & eigIndex, int & multipleEigIndex, 
                                                                  bool saveErrors /* = false */ )
{
    const int & computeEigenmodes = dbase.get<int>("computeEigenmodes");
    if( !computeEigenmodes ) 
        return 0;

    

  // -- This next call will read in any known eigenmodes ---
  // cgWave.initializeDeflation();
    if( !dbase.has_key("uev") )  
        return 0;

  // ----- determine errors in eigenmodes if we know the true eigenmodes -----
  // printF(" CgWave::getErrorInEigenPair, START...\n");


  // Here are the "true" eigenvectors and eigenvalues:
    realCompositeGridFunction & uev     = dbase.get<realCompositeGridFunction>("uev");
    RealArray & eig                     = dbase.get<RealArray>("eig");
    IntegerArray & eigMultiplicity      = dbase.get<IntegerArray>("eigMultiplicity");
    const int & useAccurateInnerProduct = dbase.get<int>("useAccurateInnerProduct");

    const int numberOfEigenvectors = uev.getComponentBound(0) - uev.getComponentBase(0) + 1;
  // printF(">> CgWave::getErrorInEigenPair: There are %d eigenvectors (true, from file).\n",numberOfEigenvectors);

    Range all;
  // uev.updateToMatchGrid(cg,all,all,all,numberOfEigenvectors); // the cg in uev may have gone out of scope --- FIX ME ??

  // int & eigIndex             = dbase.get<int>("eigIndex");

    Real diffMin = REAL_MAX*.1;
    eigIndex=0;  
    for( int i=0; i<numberOfEigenvectors; i++ )
    {
    // Note : Cannot assume sorted when searching, due to multiplicities
        Real diff = fabs( lambda - eig(0,i) ) ;
        if( diff < diffMin )
        {
              eigIndex=i; diffMin=diff;
        }
    }
  // printF("**** lambda=%g, eigIndex=%d, lambdaTrue=%g\n",lambda,eigIndex,eig(0,eigIndex));

    multipleEigIndex = eigIndex; 
    if( eigMultiplicity(eigIndex)>1 )
    {
    // Shift eigIndex to be the first of a multiple eig
        int multiplicity=eigMultiplicity(multipleEigIndex);
    // const Real eigTol=1.e-4;  // FIX ME 
        const Real eigTol = dbase.get<Real>("eigenValueTolForMultiplicity");
        for( int m=1; m<multiplicity; m++ )
        {
            if( multipleEigIndex>0 && eigMultiplicity(multipleEigIndex)==eigMultiplicity(multipleEigIndex-1 ) &&
                    fabs( eig(0,multipleEigIndex-1) - eig(0,multipleEigIndex) )/( 1.+fabs(eig(0,multipleEigIndex) )) < eigTol )
            {
                multipleEigIndex--;
            }
            else
            {
                break;
            }
        }
    }
    lambdaTrue = eig(0,eigIndex);
    relErrEigenvalue = fabs(lambda-lambdaTrue)/max(1.,fabs(lambdaTrue));




  // --- Compute the error in the eigenvector ---
    const Real eigVectorNorm = maxNorm( eigVector,component );  // should be 1   
    if( fabs(eigVectorNorm-1.) > REAL_EPSILON*1000. )
    {
        printF("INFO: computed eigVector does not have max-norm 1! eigVectorNorm=%12.5e (eigVectorNorm-1=%9.3e). Scaling...\n",eigVectorNorm,eigVectorNorm-1.);
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            eigVector[grid](all,all,all,component) *= (1./eigVectorNorm); 
        }
    // OV_ABORT("error");
    }

    if( !dbase.has_key("integrate") )
    {
        Integrate & integrate = dbase.put<Integrate>("integrate");
        if( useAccurateInnerProduct )
            integrate.updateToMatchGrid(cg);
    }
    Integrate & integrate = dbase.get<Integrate>("integrate"); 

    realCompositeGridFunction & error = dbase.get<realCompositeGridFunction>("error");
  // if( saveErrors )
  // {
  //   int numErrComp = error.getComponentBound(0)- error.getComponentBase(0)+1;
  //   if( numErrComp != numberOfEigenvectors )
  //   {
  //     printF("\n CgWave::getErrorInEigenPair ++++++ REDIMENSION ERROR GF ++++++++\n\n");
  //     error.updateToMatchGrid( cg,all,all,all,numberOfEigenvectors );
  //   }
  //   for( int ie=0; ie<numberOfEigenvectors; ie++ )
  //     error.setName(sPrintF("err%d",ie),ie);    
  // }

    realCompositeGridFunction w(cg,all,all,all); // temp space 
    Index I1,I2,I3;
    relErrEigenvector = REAL_MAX*.1;

  // printF("***INFO: eigenvalue eigIndex=%d has multiplicity=%d (norm=%9.3e)****\n",eigIndex,eigMultiplicity(eigIndex),eigVectorNorm);
    if( eigMultiplicity(eigIndex)==1 )
    {
    // ---- simple eigenvalue --

    //  Project eigenVector v onto uev 
    //  v = alpha*uev
    //  alpha = (v,uev)/(uev,euv)
        Real dotProduct;
            if( useAccurateInnerProduct )
            {
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(Real,eigVector[grid],uuLocal);
                    OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                    wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,component)*vvLocal(I1,I2,I3,eigIndex);
                }
                dotProduct = integrate.volumeIntegral(w);
            }
            else
            {
                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                int iab[2]; 
                int ii1,ii2,ii3;  
                dotProduct=0.; 
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    const IntegerArray & gid = mg.gridIndexRange();
                    OV_GET_SERIAL_ARRAY(Real,eigVector[grid],uuLocal);
                    OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                        Iv[2]=Range(0,0);
                        for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                        {
                            for( int side=0; side<=1; side++ )
                            {
                                int is = 1-2*side;
                                iab[side]=gid(side,axis);
                                const int bc = mg.boundaryCondition(side,axis);
                                if( bc==CgWave::dirichlet )
                                {
                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                }
                                else if( bc==CgWave::neumann )
                                {
                  // include boundary
                                }
                                else if( bc>0 )
                                {
                                    printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                    OV_ABORT("error");
                                }
                                else if( bc<0 )
                                {
                  // periodic -- include left end
                                    if( side==1 )
                                        iab[side] += is; 
                                }
                                else
                                {
                  // interpolation boundary : include end 
                                }
                            }
                            Iv[axis] = Range(iab[0],iab[1]);
                        }
                    bool ok=ParallelUtility::getLocalArrayBounds(eigVector[grid],uuLocal,I1,I2,I3);
                    if( ok )
                    {
                        FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                        {
                            if( maskLocal(ii1,ii2,ii3)>0 )
                                dotProduct += uuLocal(ii1,ii2,ii3,component)*vvLocal(ii1,ii2,ii3,eigIndex);
                        }
                    }
                }
                dotProduct = ParallelUtility::getSum( dotProduct );
            }  

        Real uevDot;
            if( useAccurateInnerProduct )
            {
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                    OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                    wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,eigIndex)*vvLocal(I1,I2,I3,eigIndex);
                }
                uevDot = integrate.volumeIntegral(w);
            }
            else
            {
                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                int iab[2]; 
                int ii1,ii2,ii3;  
                uevDot=0.; 
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    const IntegerArray & gid = mg.gridIndexRange();
                    OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                    OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                        Iv[2]=Range(0,0);
                        for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                        {
                            for( int side=0; side<=1; side++ )
                            {
                                int is = 1-2*side;
                                iab[side]=gid(side,axis);
                                const int bc = mg.boundaryCondition(side,axis);
                                if( bc==CgWave::dirichlet )
                                {
                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                }
                                else if( bc==CgWave::neumann )
                                {
                  // include boundary
                                }
                                else if( bc>0 )
                                {
                                    printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                    OV_ABORT("error");
                                }
                                else if( bc<0 )
                                {
                  // periodic -- include left end
                                    if( side==1 )
                                        iab[side] += is; 
                                }
                                else
                                {
                  // interpolation boundary : include end 
                                }
                            }
                            Iv[axis] = Range(iab[0],iab[1]);
                        }
                    bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                    if( ok )
                    {
                        FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                        {
                            if( maskLocal(ii1,ii2,ii3)>0 )
                                uevDot += uuLocal(ii1,ii2,ii3,eigIndex)*vvLocal(ii1,ii2,ii3,eigIndex);
                        }
                    }
                }
                uevDot = ParallelUtility::getSum( uevDot );
            }  

        Real alpha =  dotProduct/uevDot;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
      // general case 
            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],eigVectorLocal);
            OV_GET_SERIAL_ARRAY(Real,uev[grid],uevLocal);

            getIndex(cg[grid].dimension(),I1,I2,I3);
            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
            if( ok )
                wLocal(I1,I2,I3) = alpha*uevLocal(I1,I2,I3,eigIndex)  -eigVectorLocal(I1,I2,I3,component);

        }                  

        relErrEigenvector = maxNorm( w )/eigVectorNorm;

    }


    else
    {
    // -- arbitrary multiplicity ---

    // -- this came from genEigs.bC ----
        const int multiplicity=eigMultiplicity(eigIndex);
        int md=multiplicity;

    // printF("**** lambda=%g, multiplicity=%d, component=%d, for eigIndex=%d\n",lambda,multiplicity,component,eigIndex);

        const int i0=multipleEigIndex; 
        RealArray a(md,md), b(md), alpha(md);
        for( int i1=0; i1<md; i1++ )
        {
            int ii = i0 + i1;
                if( useAccurateInnerProduct )
                {
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                        OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                        OV_GET_SERIAL_ARRAY(Real,eigVector[grid],vvLocal);
                        OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                        wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,component);
                    }
                    b(i1) = integrate.volumeIntegral(w);
                }
                else
                {
                    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                    int iab[2]; 
                    int ii1,ii2,ii3;  
                    b(i1)=0.; 
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cg[grid];
                        const IntegerArray & gid = mg.gridIndexRange();
                        OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                        OV_GET_SERIAL_ARRAY(Real,eigVector[grid],vvLocal);
                        OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                            Iv[2]=Range(0,0);
                            for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                            {
                                for( int side=0; side<=1; side++ )
                                {
                                    int is = 1-2*side;
                                    iab[side]=gid(side,axis);
                                    const int bc = mg.boundaryCondition(side,axis);
                                    if( bc==CgWave::dirichlet )
                                    {
                                          iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                    }
                                    else if( bc==CgWave::neumann )
                                    {
                    // include boundary
                                    }
                                    else if( bc>0 )
                                    {
                                        printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                        OV_ABORT("error");
                                    }
                                    else if( bc<0 )
                                    {
                    // periodic -- include left end
                                        if( side==1 )
                                            iab[side] += is; 
                                    }
                                    else
                                    {
                    // interpolation boundary : include end 
                                    }
                                }
                                Iv[axis] = Range(iab[0],iab[1]);
                            }
                        bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                        if( ok )
                        {
                            FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                            {
                                if( maskLocal(ii1,ii2,ii3)>0 )
                                    b(i1) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,component);
                            }
                        }
                    }
                    b(i1) = ParallelUtility::getSum( b(i1) );
                }  

            for( int i2=0; i2<md; i2++ )
            {
                int jj = i0 + i2;
                    if( useAccurateInnerProduct )
                    {
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            wLocal(I1,I2,I3) = uuLocal(I1,I2,I3,ii)*vvLocal(I1,I2,I3,jj);
                        }
                        a(i1,i2) = integrate.volumeIntegral(w);
                    }
                    else
                    {
                        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                        int iab[2]; 
                        int ii1,ii2,ii3;  
                        a(i1,i2)=0.; 
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            const IntegerArray & gid = mg.gridIndexRange();
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],uuLocal);
                            OV_GET_SERIAL_ARRAY(Real,uev[grid],vvLocal);
                            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
                            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                Iv[2]=Range(0,0);
                                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                                {
                                    for( int side=0; side<=1; side++ )
                                    {
                                        int is = 1-2*side;
                                        iab[side]=gid(side,axis);
                                        const int bc = mg.boundaryCondition(side,axis);
                                        if( bc==CgWave::dirichlet )
                                        {
                                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                        }
                                        else if( bc==CgWave::neumann )
                                        {
                      // include boundary
                                        }
                                        else if( bc>0 )
                                        {
                                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                            OV_ABORT("error");
                                        }
                                        else if( bc<0 )
                                        {
                      // periodic -- include left end
                                            if( side==1 )
                                                iab[side] += is; 
                                        }
                                        else
                                        {
                      // interpolation boundary : include end 
                                        }
                                    }
                                    Iv[axis] = Range(iab[0],iab[1]);
                                }
                            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uuLocal,I1,I2,I3);
                            if( ok )
                            {
                                FOR_3D(ii1,ii2,ii3,I1,I2,I3)
                                {
                                    if( maskLocal(ii1,ii2,ii3)>0 )
                                        a(i1,i2) += uuLocal(ii1,ii2,ii3,ii)*vvLocal(ii1,ii2,ii3,jj);
                                }
                            }
                        }
                        a(i1,i2) = ParallelUtility::getSum( a(i1,i2) );
                    }  
            }

        }

    // PA = LU factor
        IntegerArray ipvt(md);
        int info;
        GETRF( md,md,a(0,0), md, ipvt(0), info );                
        if( info!=0 )
        {
            printF("getErrorInEigenPair:ERROR return from GETRF, info=%d (>0 means system is singular).\n",info);
      // OV_ABORT("error");
        }
        if( info==0 )
        {
      // Solve given LU
            int nrhs=1;
            GETRS( "N", md, nrhs, a(0,0), md, ipvt(0), b(0), md, info );
            alpha=b;
      // ::display(alpha,"alpha");
            if( info!=0 )
            {
                printF("getErrorInEigenPair:ERROR return from GETRS, info=%d\n",info);
                OV_ABORT("error");
            }
        }
        else
        {
            printF("getErrorInEigenPair:ERROR: system for computing the error in the eigenvector, multiplicity=%d, is singular.  \n",multiplicity);
            alpha=0;
            alpha(0)=1.; 
        }                


        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
      // WARNING: We cannot scale here since we may make both eigenvectors the same!
      // u[grid](all,all,all,i) = alpha1*u[grid](all,all,all,i) + alpha2*u[grid](all,all,all,j);
        // general case 

      // general case 
            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
            OV_GET_SERIAL_ARRAY(Real,eigVector[grid],eigVectorLocal);
            OV_GET_SERIAL_ARRAY(Real,uev[grid],uevLocal);

            getIndex(cg[grid].dimension(),I1,I2,I3);

            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
            if( ok )
            {
                wLocal(I1,I2,I3) = -eigVectorLocal(I1,I2,I3,component);
                for( int i=0; i<multiplicity; i++ )
                    wLocal(I1,I2,I3) += alpha(i)*uevLocal(I1,I2,I3,i0+i);          
            }

        } 

        Real wNorm = maxNorm( w );
        relErrEigenvector = wNorm/eigVectorNorm;

    // OV_ABORT("CgWave::getErrorInEigenPair: FINISH ME -- eigMultiplicity > 2 ")
    }

    if( saveErrors )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,error[grid],errorLocal);
            OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);
            getIndex(cg[grid].dimension(),I1,I2,I3);
            bool ok=ParallelUtility::getLocalArrayBounds(w[grid],wLocal,I1,I2,I3);
            if( ok )
                errorLocal(I1,I2,I3,component) = wLocal(I1,I2,I3);
        }  

    }


    if( false )
        printF("getErrorInEigenPair: lambda=%16.8e, eigIndex=%3d (true discrete): relErrEigenvalue=%8.2e, relErrEigenvector=%8.2e, multiplicity=%d (computed)\n",
                lambda,eigIndex,relErrEigenvalue,relErrEigenvector,eigMultiplicity(eigIndex));

    return 1;
}



// ==========================================================================================
///  \brief Evaluate the Rayleigh Quotient for a component of the grid function v
// ==========================================================================================
Real CgWave::getRayleighQuotient( realCompositeGridFunction & v, int component /* =0 */ )
{
    
  // printF("Cgwave::getRayleighQuotient:: estimate the eigenvalue using the Rayleigh quotient...\n");

    Real lambda=0.;

    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

    FILE *& debugFile  = dbase.get<FILE*>("debugFile");
    FILE *& pDebugFile = dbase.get<FILE*>("pDebugFile");  
    
  // useDiscreteInnerProduct=true  : use the Integrate class to compute integrals
  //                        =false : sum over active interior points 

  // bool useDiscreteInnerProduct = false; // true;  // These options need to be compared
    const int & useAccurateInnerProduct = dbase.get<int>("useAccurateInnerProduct");
    if( !useAccurateInnerProduct )
    {
    // ----- Compute the RQ by summing over active points ------
    // This should be ALMOST AS GOOD as using an inner product when V_j is close to an eigenvector
    //     L_h V_j = \lambda^2 V_j 


        realCompositeGridFunction lap(cg);   // ***** do this for now ... is there a work space we can use instead?

        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
        int iab[2];    

    // Compute   lap = L v 
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();

            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);

                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
          
            operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lapLocal,I1,I2,I3,component);

        }

    // Compute   ( v, vLap ) and (v, v )
        Real vDotLap=0., vDotv=0.;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();

            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);  

                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);

            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3) > 0 )
                {    
                    vDotLap += vLocal(i1,i2,i3,component)*lapLocal(i1,i2,i3);
                    vDotv   += vLocal(i1,i2,i3,component)*vLocal(i1,i2,i3,component);
                }
            }

        }
    // if( true )
    // {
    //   fprintf(pDebugFile,"getRQ: component=%3d: vDotLap=%20.13e vDotv=%20.13e (Before getSum)\n",component,vDotLap,vDotv);
    // }

    // if( false )
    // {
    //   Real sum;
    //   MPI_Allreduce(&vDotLap, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    //   vDotLap=sum;
    //   MPI_Allreduce(&vDotv,   &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    //   vDotv=sum;

    // }
        vDotLap=ParallelUtility::getSum(vDotLap);
        vDotv  =ParallelUtility::getSum(vDotv);

        lambda = vDotLap/vDotv;  // this is really -lambda^2 .. sqrt taken below 

    // if( true )
    // {
    //   fprintf(pDebugFile,"getRQ: component=%3d: vDotLap=%20.13e vDotv=%20.13e lamSq=%20.13e\n",component,vDotLap,vDotv,lambda);
    // }

    }
    else
    {
    // Compute the RQ by using discrete approximations to the L2 inner product

        if( !dbase.has_key("integrate") )
        {
            Integrate & integrate = dbase.put<Integrate>("integrate");
            integrate.updateToMatchGrid(cg);
        }
        Integrate & integrate = dbase.get<Integrate>("integrate");

        realCompositeGridFunction lap(cg);   // ***** do this for now ... is there a work space we can use instead?

        Index I1,I2,I3;

    // Compute   lap = L v 
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);

            getIndex(cg[grid].dimension(),I1,I2,I3);
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);

            operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lapLocal,I1,I2,I3,component);

        }

    // Compute ( v, vLap )
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            getIndex(mg.dimension(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);      
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
            if( ok )
                lapLocal(I1,I2,I3) = vLocal(I1,I2,I3,component)*lapLocal(I1,I2,I3);
        }

    // *** vLocal should have zero BC's -- so no need to change lap I think ****
        lap.periodicUpdate(); // is this needed ? interpolate may do this

        lap.interpolate();   // Probably do NOT want to do this if v is not smooth


        Real vDotLap = integrate.volumeIntegral(lap);    

    // Compute ( v,v )
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            getIndex(mg.dimension(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,lap[grid],lapLocal);  

            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
            if( ok )
                lapLocal(I1,I2,I3) = vLocal(I1,I2,I3,component)*vLocal(I1,I2,I3,component); // store in lapLocal
        }

        Real vDotv = integrate.volumeIntegral(lap);  

        lambda = vDotLap/vDotv;  // this is really -lambda^2 .. sqrt taken below 

    }

    if( lambda<0 )
    {
        lambda=sqrt(-lambda);
        if( false )
            printF("CgWave::getRayleighQuotient: component=%d: lambda(RQ)=%16.8e (Lap(phi) = -lambda^2 phi)\n",
                      component,lambda);
    }
    else
    {
        printF("Cgwave::getRayleighQuotient::ERROR: rayleigh quotient for Laplacian is postive!? lambda=%12.4e\n",lambda);
        OV_ABORT("error");
    }

    return lambda;
}



// ==========================================================================================
///  \brief Compute the relative residual in the eigenvalue equation:
///      || L_h v + lambda^2 v ||/lambda^2 
/// 
// ==========================================================================================
Real CgWave::getEigenPairResidual( Real lambda, realCompositeGridFunction & v,
                                                                      realCompositeGridFunction & res,  int component /* =0 */ )
{

    const Real & c  = dbase.get<real>("c");
    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
    

  // realCompositeGridFunction res(cg);   // ***** do this for now ... is there a work space we can use instead?

    Index I1,I2,I3;
    Index D1,D2,D3;
    int i1,i2,i3;

  // Compute   lap = L v 
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
        OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
        OV_GET_SERIAL_ARRAY(Real,res[grid],resLocal);


    // Compute at active interior points **CHECK ME ***  avoid interp points
        int extra=-1; // ** check me ** fix for periodic or interp boundaries
        getIndex(cg[grid].gridIndexRange(),I1,I2,I3,extra);

        getIndex(cg[grid].dimension(),D1,D2,D3);  
        bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,D1,D2,D3);
        RealArray lapLocal(D1,D2,D3);

        resLocal(D1,D2,D3,component)=0.;

        ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
        if( ok )
        {
            operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lapLocal,I1,I2,I3,component);
      // res = c^2 Lap (v) + lambda^2 * v
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3)>0 )
                {
                    resLocal(i1,i2,i3,component) = (c*c)*lapLocal(i1,i2,i3) + (lambda*lambda) * vLocal(i1,i2,i3,component);
                }
                else
                {
                    resLocal(i1,i2,i3,component)=0.;
                }
            }
        }
    }
    Real maxRes;

    maxRes = maxNorm( res,component )/(lambda*lambda);

  // bool plotResiduals = false; // true; // false; // Make an option 
  // if( plotResiduals )
  // {
  //   printF("CgWave::getEigenPairResidual: lambda=%10.3e,component=%d,  maxRes=%9.2e\n",lambda,component,maxRes);
  //   res.setName("res",0);

  //   GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("test"));

  //   ps.erase();
  //   PlotStuffParameters psp;
  //   psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  //   psp.set(GI_COMPONENT_FOR_CONTOURS,component);
  //   PlotIt::contour(ps,v,psp);
  //   ps.erase();
  //   PlotIt::contour(ps,res,psp);
  //   psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true); 
  // }
    return maxRes;


}





