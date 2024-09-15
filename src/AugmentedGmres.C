// This file automatically generated from AugmentedGmres.bC with bpp.
#include "AugmentedGmres.h"

#include "display.h"

#define DGEQRF EXTERN_C_NAME(dgeqrf)
#define DORGQR EXTERN_C_NAME(dorgqr)

extern "C" 
{
  // QR factorization
    void DGEQRF( const int & M,
                              const int  &  N,
                              const double &  A,
                              const int  &  LDA,
                              double &  TAU,
                              double &  WORK,
                              int & LWORK,
                              int & INFO ); 

  // Recover Q from QR factorization
    void DORGQR( const int & M, const int & N, const int & K, Real & A, const int & LDA, Real & TAU, Real & WORK, int & LWORK, int & INFO );

}

// ===================================================================================================
/// \brief Constructor for the Augmented GMRES algorithm
// ===================================================================================================
AugmentedGmres::AugmentedGmres()
{

    dbase.put<int>("numberOfIterations")=0;
    dbase.put<RealArray>("resVect");
    dbase.put<int>("augmentedVectorsAreEigenvectors")=false;
    dbase.put<RealArray>("eig");  // holds eigenvalues for augmented vectors
    dbase.put<Real>("residual")=0.;  // holds residual from last solve

}

// ===================================================================================================
/// \brief Destructor for the Augmented GMRES algorithm
// ===================================================================================================
AugmentedGmres::~AugmentedGmres()
{

}

// ===================================================================================================
/// \brief  return the residual from the last solve
// ===================================================================================================
Real AugmentedGmres::getResidual() const
{
    return dbase.get<Real>("residual");
}


// ===================================================================================================
/// \brief Ssupply eigenvalues if augmented vectors are eigenvectors
/// \param augEigs (input) : eigenvalues for the augmented vectors
// ===================================================================================================
// 
int AugmentedGmres::setAugmentedEigenvalues( const RealArray & augEigs )
{
    int & augmentedVectorsAreEigenvectors = dbase.get<int>("augmentedVectorsAreEigenvectors");
    augmentedVectorsAreEigenvectors=true;

    RealArray & eig = dbase.get<RealArray>("eig"); 
    eig.redim(0);
    eig = augEigs;
}

//
// Compute y = A*x 
// \param transpose : 1 = compute y = A^Y x
// \param m0, n0 : if specified use this as the dimensions of A
//
void AugmentedGmres::matVect(const RealArray & A, const RealArray & x, RealArray & y, int transpose /* =0 */, int m0 /* =-1 */, int n0 /* =-1 */ )
{
    const int m = m0==-1 ? A.getLength(0) : m0;
    const int n = n0==-1 ? A.getLength(1) : n0;
  // printF("matVec: m=%d n=%d transpose=%d\n",m,n,transpose);
    if( transpose==0 )
    {
        for( int i=0; i<m; i++ )
        {
            Real temp=0.;
            for( int j=0; j<n; j++ )
                temp += A(i,j)*x(j);

            y(i) = temp; 
        }
    }
    else
    {
        for( int j=0; j<n; j++ )
        {
            Real temp=0;
            for( int i=0; i<m; i++ )
                temp += A(i,j)*x(i);
            y(j) = temp; 
        }    

    }
} 


// void matVectMatrixFree( const RealArray & x, RealArray & y )
// {
//   OV_ABORT("FINISH ME");

//   // const int m = x.getLength(0);
//   // for( int i=0; i<m; i++ )
//   // {
//   //   Real temp=0;
//   //   for( int j=0; j<m; j++ )
//   //   {
//   //     temp += A(i,j)*x(j);
//   //   }
//   //   y(i) = temp; 
//   // }
// }

// Compute the 2-norm of x 
Real AugmentedGmres::norm( const RealArray & x )
{
    const int m = x.getLength(0);
    Real temp=0; 
    for( int i=0; i<m; i++ )
        temp += SQR(x(i));

    return sqrt(temp);
} 

// Compute the inner product (x,y)
Real AugmentedGmres::innerProduct( const RealArray & x, const RealArray & y )
{
    const int m = x.getLength(0); 
    Real dot=0.;
    for( int i=0; i<m; i++ )
        dot += x(i)*y(i);

    return dot;
}

int AugmentedGmres::getNumberOfIterations() const
{
    return dbase.get<int>("numberOfIterations");
}

RealArray & AugmentedGmres::getResidualVector() const
{
    return dbase.get<RealArray>("resVect");
}



// ===================================================================================================
/// \brief Solve A x = b using the Augmented GMRES algorithm
/// \param x (output) 
/// \Return value: relative residual || b - A x ||/|| b ||
// ===================================================================================================
Real AugmentedGmres::solve( const RealArray & A, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, 
                                                        RealArray & x )
{
    return solve( A, NULL, b, x0, W, maxit, tol, x );
}


// use a matrix-vector multiply function
Real AugmentedGmres::solve( MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x )
{
    RealArray A; 
    return solve( A, matVectFunction, b, x0, W, maxit, tol, x );
}


// ===================================================================================================
/// \brief Solve A x = b using the Augmented GMRES algorithm
/// \param x (output) 
/// \Return value: relative residual || b - A x ||/|| b ||
// ===================================================================================================
Real AugmentedGmres::solve( const RealArray & A, MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, 
                                                        RealArray & x )
{


    const int matrixIsAFunction = matVectFunction==NULL ? 0 : 1; 

    const int m = b.getLength(0);
    const int numToDeflate = W.getLength(1);

    printf("AugmentedGmres::solve: m=%d, numToDeflate=%d, matrixIsAFunction=%d, maxit=%d, tol=%9.2e\n",m,numToDeflate,matrixIsAFunction,maxit,tol);

    Real resid =0.;

    

    Real bNorm = norm(b);

  // A*(x+x0) = b0 
  // A*x = b0 - A*x0
    RealArray b0(m);
    RealArray y(m);  

    if( matrixIsAFunction )
          matVectFunction( x0, y ); // y = A*x0
    else
        matVect( A, x0, y ); // y = A*x0

    b0 = b - y;
  // ::display(b0,"b0","%6.2f ");

    const int p = numToDeflate;

    int n = min(p+maxit+1,m); // max number of columns in V -- could do better here
    RealArray V(m,n+1);
    RealArray H(n+1,n+1);

    H=0; 
    V=0; 
  

    int & augmentedVectorsAreEigenvectors = dbase.get<int>("augmentedVectorsAreEigenvectors");
  // augmentedVectorsAreEigenvectors=false;

    const RealArray & eig = dbase.get<RealArray>("eig");
    if( augmentedVectorsAreEigenvectors )
    {
        if( eig.getLength(0) < p )
        {
            printF("AUG-GMRES: eigenvalues have been supplied for augmented eigenvectors but there are not enough!\n"
                          "  supplied eignvalues = %d, numberOfAugmented vectors =%d\n",eig.getLength(0),p); 
            OV_ABORT("ERROR");
        }
    }

    Range M=m;

    RealArray vp1(m);
    RealArray y1(m);  // temp space
    vp1=b0;
    if( p>0 )
    {
        RealArray AW(m,p);  // holds A*W 
        for( int j=0; j<p; j++ )
        {
            if( augmentedVectorsAreEigenvectors )
            {
                printF("AUG-GMRES: Compute A*W(:,%d) ASSUMING augmented vectors are eigenvectors, eig=%10.3e, lam=1-eig=%10.3e\n",j,eig(j),1.0-eig(j));
                
        // **NOTE** A*w = (I-S)*w 
                AW(M,j) = (1.0-eig(j))*W(M,j); // A*W(:,j) = (I-S)  W(:,j)  

        // -- now check:
                bool checkEigenvectorMatVect = false; // set to true to check  
                if( checkEigenvectorMatVect )        
                {
                    y1 = W(M,j);
                    if( matrixIsAFunction ) 
                        matVectFunction( y1, y ); // y = A*y1 
                    else
                        matVect( A,y1, y );

                    Real maxErr = max(fabs( AW(M,j) - y ))/ max(fabs(y1));
                    printF("++AUG-GMRES: j=%d: max | (1-eig(j))*W(:,j) - A*W(:,j) |/|w_j|_inf = %9.2e\n",j,maxErr);

                    Real lamEst1 =  sum( y/y1 )/m; 
                    Real lamEst2 =  y(5)/y1(5);

                    Real residEst1 = max(fabs( y - lamEst1*y1 ))/ max(fabs(y1));
                    Real residEst2 = max(fabs( y - lamEst2*y1 ))/ max(fabs(y1));
                        
                    printF("++AUG-GMRES: j=%d: eig=%12.5e, lam=1-eig=%12.5e, lamEst=%12.5e or %12.5e, residEst=[%8.2e,%8.2e]\n",j,eig(j),1.-eig(j),lamEst1,lamEst2,residEst1,residEst2); 
                    
                }      
            }
            else
            {      
        // printF("Compute A*W(:,%d)..\n",j);
                y1 = W(M,j);
        // ::display(y1,"W(:,j)","%6.2f ");

                printF("AUG-GMRES: Compute A*W(:,%d) (Augmented vectors)\n",j);

                if( matrixIsAFunction ) 
                    matVectFunction( y1, y ); // y = A*x0
                else
                    matVect( A,y1, y );

        // ::display(y,"A*W(:,j)","%6.2f ");
                AW(M,j) = y;
        // AW(M,j) = A(W(:,j));
            }

        }
    // if( augmentedVectorsAreEigenvectors )
    //   OV_ABORT("AUG-GMRES: Stop here for now");

    // printF("AugmentedGmres::solve:Compute W = QR ...\n");
    // ::display(W(M,Range(p)),"W","%6.2f ");
    // ::display(AW(M,Range(p)),"AW","%6.2f ");

    // [Vp,H] = qr(AW,0); % reduced QR factorization of A * (matrix of augmented vectors) 
        int lda=m, info=0;
        RealArray tau(p);
        int lwork= p*4; // n*nb nb = optimum block size
        RealArray work(lwork);
        DGEQRF( m,p,AW(0,0),lda,tau(0),work(0),lwork,info );
        if( info<0 )
        {
            printF("AugmentedGmres::solve: Error return from DGQRF = info=%d\n",info);
            OV_ABORT("error");
        }
        for( int i=0; i<p; i++ )
        {
            for( int j=i; j<p; j++ )
            {
                H(i,j)=AW(i,j);         // Upper triangular part of AW holds "R"
            }
        }
    // Eval "Q" 
        DORGQR( m, p, p, AW(0,0), lda,tau(0),work(0),lwork,info );  // on output AW holds Q 
        if( info<0 )
        {
            printF("AugmentedGmres::solve: Error return from DORGQR = info=%d\n",info);
            OV_ABORT("error");
        }  
        for( int j=0; j<p; j++ ) 
        {
            V(M,j) = AW(M,j); 
        } 

    // ::display(H(Range(p),Range(p)),"H after A*W=Vp*H","%6.2f ");
    // ::display(V(M,Range(p)),"Vp after A*W=Vp*H","%6.2f ");

        for( int i=1; i<=2; i++ )  // do twice for re-orthogonalization
        {
      //     vp1 = vp1 - Vp*(Vp'*vp1);

            int transpose=1; 
            matVect( V,vp1,y, transpose, m,p  );
      /// ::display(y(Range(p))," y = V^T vp1","%6.2f ");
            matVect( V,y, y1, 0, m,p  );
      // ::display(y1," y1 = V*(V^T vp1)","%6.2f ");
            vp1 -= y1; 
            
        }
    }
    vp1 = vp1/norm(vp1);

  // ::display(vp1,"vp1 at start","%6.2f ");

  // V = zeros(m,p+maxit+1);   % allocate space 

  // if( p>0 )
  //   V(:,1:p) = Vp;
  // end
    V(M,p) = vp1;  

  // % Q holds the products of Givens rotations
  // Q = eye(m,p+1);
  // Hg = H; % holds a copy of H for testing Givens

    RealArray Q(n+1,n+1); // holds the products of Givens rotations
    Q=0;
    for( int i=0; i<p+1; i++ )
        Q(i,i)=1; 


    RealArray g(n+1);
    g=0; 
    for( int k=0; k<p+1; k++ )
    {
        y = V(M,k); 
        g(k) = innerProduct( y,b0  );
  //   g(k,1) = V(:,k)' * b0; 
    }

  // ::display(g,"g at start","%6.2f ");

  // ---------- Arnoldi ---------
  // printF("AugmentedGmres::solve: start Arnoldi iterations bNorm=%9.2e...\n",bNorm);

    RealArray & resVect = dbase.get<RealArray>("resVect");
    resVect.redim(maxit+1);
    resVect=0; 

    RealArray resv(maxit+p);
    resv=0; 
  // iteration : counts actual number of mat-vects 
  //           : don't count augmented vectors if we do need to use a mat-vect
    int iteration= augmentedVectorsAreEigenvectors ? -1 : p-1;

    int it=0;

    for( int k=p; k<p+maxit; k++ )
    {

    // printF("AugmentedGmres::solve: Arnoldi iteration k=%d\n",k);
        iteration++;

        y1 = V(M,k);
        if( matrixIsAFunction ) 
            matVectFunction( y1, vp1 ); // vp1 = A*V(:,k)
        else
            matVect( A,y1, vp1 );

        for( int i=0; i<=k; i++ )
        {
       // H(i,k) = V(:,i)'*vkp1; 
       // vkp1 = vkp1 -H(i,k)*V(:,i);
              y = V(M,i); 
              H(i,k) = innerProduct( y,vp1 );

       // printF("H(%d,%d)=%9.3e\n",i,k,H(i,k));

              vp1 -= H(i,k)*V(M,i);
        }
    // re-orthogonalization -- optional -- do for now 
        for( int i=0; i<=k; i++ )
        {
      //     alpha= V(:,i)' * vkp1;
      //     vkp1 = vkp1 - alpha*V(:,i);
      //     H(i,k)=H(i,k)+alpha;    
            y = V(M,i); 
            Real alpha = innerProduct( y,vp1 );
            vp1 = vp1 - alpha*V(M,i);
            H(i,k)=H(i,k)+alpha;
          
        } 

        H(k+1,k) = norm(vp1);
        if( fabs(H(k+1,k)) < 1e-10*bNorm )
        {
            printF("INFO: in AuGmres: H(k+1,k)=%9.2e *happy break down*. Solution must have converged.\n",H(k+1,k));
      // break;
      // OV_ABORT("FIX ME"); 
        }
        else
        {
            V(M,k+1) = vp1/H(k+1,k);  
        }

    // printF("H(%d,%d)=%9.3e\n",k+1,k,H(k+1,k));
    // ::display(V(M,k+1),"New column in V","%6.2f ");

    //   H(1:k,k) = Q(1:k,1:k)*H(1:k,k); % Apply previous rotations to part of new column 
        Range K(0,k);
        y1(K) = H(K,k);
        matVect( Q,y1,y,  0,k+1,k+1 );
        H(K,k) = y(K); 
    // ::display(H(K,k),"New column in H","%6.2f ");

        Real rho = H(k,k);
        H(k,k) = sqrt( SQR(rho) + SQR(H(k+1,k)) ) ;
        Real c = rho/H(k,k);
        Real s = H(k+1,k)/H(k,k);
    // printF("c=%9.3e s=%9.3e\n",c,s);
        H(k+1,k) = 0;
        
    // Apply Givens rotation to Q 
        Q(k+1,K) = -s*Q(k,K);
        Q(k  ,K) =  c*Q(k,K);
        Q(k+1,k+1) = c;
        Q(k  ,k+1) = s;
        
    // Apply Givens rotation to g 
        g(k+1) = -s*g(k);  // this is the current 2-norm of the residual
        g(k  ) =  c*g(k); 

        it=k-p; // current Arnoldi iteration

    //  current L2 residual is g(k+1)

        Real resid = abs(g(k+1))/sqrt(m); // do this to match other solvers
    // Real resid = abs(g(k+1))/bNorm; 

        if( iteration==p && !augmentedVectorsAreEigenvectors )
            resv(Range(0,p))=resid; // just set residuals for augmented steps to the first residual
        else
            resv(iteration)=resid; 

        resVect(it)= resid;
        if( it>0 )
        {
            printF("AUG-GMRES: it=%3d, || r ||_2h = %8.2e, ratio=%8.2e, (matrixIsAFunction=%d)\n",it+1,resid,resVect(it)/resVect(it-1),matrixIsAFunction);
        }
    
        if( resid < tol )
        {
            break; 
        }

    }

    const int jp = p+it+1; // total number of iterations

  // printF("AGmres: jp=%d, iteration=%d\n",jp,iteration);
  // ::display(resv,"resv (*NEW WAY*)");

  // ::display(H(Range(jp),Range(jp)),"H for triangular solve","%6.2f ");

  // ---- y = H(1:jp,1:jp)\g(1:jp); :  upper triangular solve ----
    for( int i=jp-1; i>=0; i-- )
    {
        Real temp = g(i);
        for( int j=i+1; j<jp; j++ )
            temp -= H(i,j)*y(j); 
        y(i) = temp/H(i,i); 
    }
  // ::display(y(Range(jp)),"Solution y","%6.2f ");

  // x = x0 + [W , V(:,p+1:p+it) ]*y;

    for( int i=0; i<m; i++ )
    {  
        Real temp=0.;
        for( int j=0; j<p; j++ )
            temp += W(i,j)*y(j);     
        for( int j=p; j<jp; j++ )
            temp += V(i,j)*y(j);

        x(i) = x0(i) + temp;
    }
  // ::display(x,"Solution x","%6.2f ");

    if( true )
    {
     // *new* way
          int & numberOfIterations = dbase.get<int>("numberOfIterations");
          numberOfIterations = iteration+1;
          resVect.redim(numberOfIterations);
          resVect=resv(Range(0,numberOfIterations-1));
    }
    else
    {
        resVect.resize(it+1);
        ::display(resVect,"resv (*OLD WAY*)");

        dbase.get<int>("numberOfIterations")=it+1;

        dbase.get<Real>("residual")=resVect(it); // save the final residual
    }


    return resid;
}
