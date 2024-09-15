// This file automatically generated from testAugmentedGmres.bC with bpp.
//===============================================================================
//
//  Test the Augmented GMRES algorithm
//  set testag = $HOME/Dropbox/research/cgwave/bin/testAugmentedGmres 
//==============================================================================


#include "AugmentedGmres.h"
#include "display.h"


// lapack routines
#ifdef OV_USE_DOUBLE
    #define GETRF EXTERN_C_NAME(dgetrf)
    #define GETRS EXTERN_C_NAME(dgetrs)
  // #define GECON EXTERN_C_NAME(dgecon)
  // #define LANGE EXTERN_C_NAME(dlange)
#else
    #define GETRF EXTERN_C_NAME(sgetrf)
    #define GETRS EXTERN_C_NAME(sgetrs)
  // #define GETRI EXTERN_C_NAME(sgetri)
  // #define GECON EXTERN_C_NAME(sgecon)
  // #define LANGE EXTERN_C_NAME(slange)
#endif

extern "C"
{
    void GETRF( int & m, int & n, real & a, const int & lda, int & ipvt, int & info );

    void GETRS( char *trans, const int & n, const int & nrhs, const Real & a, const int & lda, int & ipiv, Real & b, const int & ldb, const int & info );
  // void GETRI( int & n, real & a, const int & lda, const int & ipvt, real & work, const int & iwork, int & info );

  // void GECON( char *norm, int & n, real & a, const int & lda, real & anorm, real & rcond, real & work, int & iwork, int & info );
  // real LANGE( char *norm, int & m, int & n, real & a, const int & lda, real & work );

  //  void GEEV( char *jobvl, char* jobvr, int & n, real & a, const int & lda,
  //             real & wr, real & wi, real &vl, int & ldvl, real & vr, int & ldvr, real & work, int & lwork, int & info );

}


RealArray *Aptr=NULL; // set below to point to the matrix A 

void matVectFunction( const RealArray & x, RealArray & y )
{
    RealArray & A = *Aptr;

    const int m = A.getLength(0);
    const int n = A.getLength(1);
  // printF("matVec: m=%d n=%d transpose=%d\n",m,n,transpose);
    for( int i=0; i<m; i++ )
    {
        Real temp=0.;
        for( int j=0; j<n; j++ )
            temp += A(i,j)*x(j);

        y(i) = temp; 
    }

} 

// ================================== MAIN =================================================
int 
main(int argc, char *argv[])
{
    Overture::start(argc,argv);  // initialize Overture


    printF("Usage: testAugmentedGmres -Nx=<i> -nit=<i> -numToDeflate=<i> -tol=<f>\n");

    int debug=0; 
    int Nx=100; 
    int maxit=50;
    int numToDeflate=10;
    int includeBCs = 0;
    int initialGuess=0; // 1 = use a non-zero initial guess
    int useFunction=0;  // 1 = use matrix free mat-vec function
    Real tol = 1e-10;

    
    int len=0;
    if( argc >= 1 )
    { 
        for( int i=1; i<argc; i++ )
        {
            aString arg = argv[i];
            if( (len=arg.matches("-debug=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&debug);
                printF("Setting debug=%i\n",debug);
            }
            else if( (len=arg.matches("-Nx=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&Nx);
                printF("Setting Nx=%i\n",Nx );
            }
            else if( (len=arg.matches("-maxit=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&maxit);
                printF("Setting maxit=%i\n",maxit );
            }
            else if( (len=arg.matches("-numToDeflate=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&numToDeflate);
                printF("Setting numToDeflate=%i\n",numToDeflate );
            } 
            else if( (len=arg.matches("-initialGuess=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&initialGuess);
                printF("Setting initialGuess=%i\n",initialGuess );
            }
            else if( (len=arg.matches("-includeBCs=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&includeBCs);
                printF("Setting includeBCs=%i\n",includeBCs );
            } 
            else if( (len=arg.matches("-useFunction=")) )
            {
                sScanF(arg(len,arg.length()-1),"%i",&useFunction);
                printF("Setting useFunction=%i\n",useFunction );
            }                       
            else if( (len=arg.matches("-tol=")) )
            {
                sScanF(arg(len,arg.length()-1),"%e",&tol);
                printF("Setting tol=%i\n",tol );
            }                   
      // else if( (len=arg.matches("-orderOfAccuracy=")) )
      // {
      //   sScanF(arg(len,arg.length()-1),"%i",&orderOfAccuracy);
      //   printF("Setting orderOfAccuracy=%i\n",orderOfAccuracy);
      // }      
            else
            {
                printF("Unknown command line argument=[%s]\n",(const char*)arg);
            }
        }
    }

    printF("=================================================================================\n"
                  " ---  test the Augmented GMRES algorithm ---  \n"
                  "   Nx=%d, maxit=%d, numToDefalate=%d, includeBCs=%d, initialGuess=%d, useFunction=%d tol=%9.2e \n",
                  Nx,maxit,numToDeflate,includeBCs, initialGuess,useFunction, tol);



  // eigenfunction sin( k*pi*x )
  // eigenvalue is   -4*sin^2(k*pi*dx/2)/dx^2 
  //   |----+----+----+---- ... ---+----|
  //   0    1    2   3                  Nx
  //   ia    
  // 
  //   dx = (bx-ax)/N
  //   N = Nx+1
  // 
    Real ax=0, bx=1;
    Real dx = (bx-ax)/Nx; 

    int ia=0, ib=Nx;
    int Ng = ib-ia+1; 

    RealArray xg(Ng);
    for( int i=ia; i<=ib; i++ )
    {
        xg(i) = ax + dx*(i-ia); 
    }

    int i1,i2;       // matrix includes these point
    if( includeBCs )
    {
        i1=ia; i2=ib;
    }
    else
    {
        i1=ia+1; i2=ib-1; 
    }
    

    Real beta=5; 
    Real xe = .4; 
    
    int m=i2-i1+1; 
    RealArray A(m,m);
    RealArray b(m);
    A=0.; 
    Real dxSq = dx*dx; 
    for( int i=0; i<m; i++ )
    {
        if( i-1>=0 ){ A(i,i-1)= -1./dxSq; }
                                    A(i,i  ) = +2/dxSq;
        if( i+1<m  ){ A(i,i+1) =-1./dxSq; }

        int j=i+i1; // true index
        b(i) = (-((-2.*beta+SQR(2.*beta*(xg(j)-xe)))*(exp(-beta*(xg(j)-xe)*(xg(j)-xe)))));
    }

    if( includeBCs )
    {
        Real bcScale=1./dxSq; 
        A(ia,ia)=bcScale;  A(ia,ia+1)=0.; b(ia) = bcScale*(exp(-beta*(ax-xe)*(ax-xe)));
        A(ib,ib)=bcScale;  A(ib,ib-1)=0.; b(ib) = bcScale*(exp(-beta*(bx-xe)*(bx-xe)));
    }
    else
    {
    // adjust RHS with BC
        b(0  ) += (exp(-beta*(ax-xe)*(ax-xe)))/dxSq;
        b(m-1) += (exp(-beta*(bx-xe)*(bx-xe)))/dxSq; 
    }

    ::display(transpose(A),"A (transposed for output to look like a Matrix)","%7.2f ");

    Aptr = &A; // for matVectFunction above 

  // ----- Discrete eigenvectors ----
    RealArray vTrue(Nx+1,m);
    Range Rx = Nx+1;
    Range I(i1,i2); 
    for( int i=0; i<m; i++ )
    {
        Real freq= (i+1);   //  adjust the frequency 
        vTrue(Rx,i) = sin(freq*Pi*xg(Rx));
    }

    RealArray W;
    if( numToDeflate>0 )
    {
        W.redim(m,numToDeflate);
        Range M = m;
        for( int j=0; j<numToDeflate; j++ )
        {
            W(M,j) = vTrue(I,j);
        }
    }
    
    RealArray x0(m), r0(m); 
    if( initialGuess==0 )
    {
        x0 = 0;
    }
    else
    {
        x0 = (exp(-beta*(xg(I)-xe)*(xg(I)-xe)));
    }
    RealArray y(m); 
    AugmentedGmres::matVect( A,x0, y); 
    r0 = b - y; 

    Real resid = AugmentedGmres::norm(r0);
    printF("Init: resid=%9.2e\n",resid);

    AugmentedGmres auGmres;

    RealArray x(m); // put solution here

  // ========== CALL AUGMENTED GMRES =========
    if( useFunction==0 )
    {
    // pass matrix
        auGmres.solve( A, b, x0, W, maxit, tol, x );
    }
    else
    {
    // Use a function to compute the matrix vector product
        auGmres.solve( matVectFunction, b, x0, W, maxit, tol, x );
    }


    ::display(x,"Done: solution x","%7.3f ");

    int numberOfIterations = auGmres.getNumberOfIterations();
    printF("number of GMRES iterations=%d\n",numberOfIterations);

    const RealArray & resVect = auGmres.getResidualVector();
    ::display(resVect, "resVect","%9.2e ");


  // --- Check errors ----
    IntegerArray ipvt(m);
    int info;
    GETRF( m, m, A(0,0), m, ipvt(0), info );
    RealArray xTrue(m);
    xTrue=b;

    int nrhs=1; 
    GETRS( "N", m, nrhs, A(0,0), m, ipvt(0), xTrue(0), m, info );

    ::display(xTrue,"xTrue","%7.3f ");

    Real maxErr;

    maxErr = max(fabs(x-xTrue));
    printF("maxErr=%9.2e\n",maxErr);





    Overture::finish();          

    return(0);
}




