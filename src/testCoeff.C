// This file automatically generated from testCoeff.bC with bpp.
// ---------------------------------------------------------------------------
//
// Test routine to check IMPLICIT COEFFICIENTS
//
//    CHECK SQUARE OF LAPLACIAN COEFF FOR FOURTH ORDER IMPLICIT
//
// ---------------------------------------------------------------------------

#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "CgWave.h"
#include "SparseRep.h"


#include "OGTrigFunction.h"
#include "OGPolyFunction.h"
#include "display.h"

// #include "incompressible/SmCylinderExactSolution.h"
#include "ParallelUtility.h"
// #include "Oges.h"
// #include "CgSolverUtil.h"

#include "gridFunctionNorms.h"


// forward declaration
int coefficientsByDelta( CompositeGrid & cg, realArray & coeff, int grid, Index Rv[3], int coeffOption =4 );



#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )
int 
getLineFromFile( FILE *file, char s[], int lim);


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



// ============================================================================================
// Declare variables needed for Lap^2 in curvilinear coordinates
// Input: 
//   R = r,s,t
//   N = 0,1,2  
// ============================================================================================  



// // ============================================================================================
// // Macro to add foruth order correction terms to the implicit matrix
// //  
// // Add modified equation term 
// //        (L_2h)^2  
// // ============================================================================================
// #beginMacro addFourthOrderCorrectionTermsToMatrixOLD()
//   if( isRectangular )
//   {
//     // 2D: 
//     // (Delta_2h)^2  = (D+xD-x)^2 + (D+yD-y)^2 + 2 (D+xD-x)(D+yD-y)
//     // 3D 
//     // (Delta_2h)^2  = (D+xD-x)^2 + (D+yD-y)^2 + (D+zD-z)^2 + 2 (D+xD-x)(D+yD-y) + 2 (D+xD-x)(D+zD-z)+ 2 (D+yD-y)(D+zD-z)

//     Range R5 = Range(-2,2);
//     RealArray lapSq(R5,R5,R5);
//     lapSq = 0.;

//     const Real dx4 = dx[0]*dx[0]*dx[0]*dx[0];
//     const Real dy4 = dx[1]*dx[1]*dx[1]*dx[1];
//     lapSq(-2,0,0) +=  1./dx4;  lapSq(0,-2,0) +=  1./dy4;
//     lapSq(-1,0,0) += -4./dx4;  lapSq(0,-1,0) += -4./dy4;
//     lapSq( 0,0,0) +=  6./dx4;  lapSq(0, 0,0) +=  6./dy4;
//     lapSq( 1,0,0) += -4./dx4;  lapSq(0, 1,0) += -4./dy4;
//     lapSq( 2,0,0) +=  1./dx4;  lapSq(0, 2,0) +=  1./dy4;

//     const Real dxy2 = dx[0]*dx[0]*dx[1]*dx[1];
//     lapSq(-1, 1,0) += 2./dxy2; lapSq(0, 1,0) +=-4./dxy2; lapSq(1, 1,0) += 2./dxy2;
//     lapSq(-1, 0,0) +=-4./dxy2; lapSq(0, 0,0) += 8./dxy2; lapSq(1, 0,0) +=-4./dxy2;
//     lapSq(-1,-1,0) += 2./dxy2; lapSq(0,-1,0) +=-4./dxy2; lapSq(1,-1,0) += 2./dxy2;

//     if( numberOfDimensions==3 )
//     {
//       const Real dz4 = dx[2]*dx[2]*dx[2]*dx[2];
//       lapSq(0,0,-2) +=  1./dz4;  
//       lapSq(0,0,-1) += -4./dz4;  
//       lapSq(0,0, 0) +=  6./dz4;  
//       lapSq(0,0, 1) += -4./dz4;  
//       lapSq(0,0, 2) +=  1./dz4; 

//       const Real dxz2 = dx[0]*dx[0]*dx[2]*dx[2];
//       lapSq(-1,0, 1) += 2./dxz2; lapSq(0,0, 1) +=-4./dxz2; lapSq(1,0, 1) += 2./dxz2;
//       lapSq(-1,0, 0) +=-4./dxz2; lapSq(0,0, 0) += 8./dxz2; lapSq(1,0, 0) +=-4./dxz2;
//       lapSq(-1,0,-1) += 2./dxz2; lapSq(0,0,-1) +=-4./dxz2; lapSq(1,0,-1) += 2./dxz2;

//       Real dyz2 = dx[1]*dx[1]*dx[2]*dx[2];
//       lapSq(0,-1, 1) += 2./dyz2; lapSq(0,0, 1) +=-4./dyz2; lapSq(0,1, 1) += 2./dyz2;
//       lapSq(0,-1, 0) +=-4./dyz2; lapSq(0,0, 0) += 8./dyz2; lapSq(0,1, 0) +=-4./dyz2;
//       lapSq(0,-1,-1) += 2./dyz2; lapSq(0,0,-1) +=-4./dyz2; lapSq(0,1,-1) += 2./dyz2;
//     }

//     // coefficients in implicit time-stepping  
//     //  D+t D-t u =              c^2 Delta( cImp(1,0) *u^{n+1} + cImp(0,0) *u^n + cImp(-1,0)* u^{n-1} )   :  second-order coeff cImp(-1:1,0)
//     //              -(c^4*dt^2/12) Delta^2( cImp(1,1) *u^{n+1} + cImp(0,1) *u^n + cImp(-1,1)* u^{n-1}  )  :  fourth-order ceoff cImp(-1:1,1) 
//     // For accuracy the weights depend on one parameter beta2 for second-order,
//     // and a second parameter beta4 for fourth-order: (See notes in research/timeStepping/implicitTaylorSchemes.pdf)

//     const Real cLapSq = cImp(-1,1)*(c*dt)*(c*dt)*(c*dt)*(c*dt);  // CHECK ME 
//     ForStencil(m1,m2,m3)
//     {
//       int m  = M123(m1,m2,m3);   
//       coeffLocal(m,I1,I2,I3) += cLapSq*lapSq(m1,m2,m3);
//     }

//   }
//   else
//   {
//     // OV_ABORT("implicit matrix: finish me for order=4 CURVLINEAR");

//     printF("implicit matrix: order=4 CURVLINEAR : **CHECK ME**\n");

//     // 

//     //    RealArray lapSqCoeff(M0,I1,I2,I3);             // *************** TEMP FOR TESTING *************

//     OV_GET_SERIAL_ARRAY(Real,lapSqCoeff[grid],lapSqCoeffLocal); 

//     OV_GET_SERIAL_ARRAY(Real,mg.center(),xLocal);  // *************** TEMP FOR TESTING *************

//     OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rxLocal);
//     // macro to make the rxLocal array look 5-dimensional 
//     #define DD(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2))     

//     // ----- COMPUTE DERIVATIVES OF METRICS -----
//     mgop.setOrderOfAccuracy(4);  // ***************
//     Range Rd2=numberOfDimensions*numberOfDimensions;
//     RealArray ddx(I1,I2,I3,Rd2), ddy(I1,I2,I3,Rd2);
//     mgop.derivative( MappedGridOperators::xDerivative,rxLocal,ddx,I1,I2,I3,Rd2);            
//     mgop.derivative( MappedGridOperators::yDerivative,rxLocal,ddy,I1,I2,I3,Rd2);

//     const int extra=1; 
//     Index J1,J2,J3;

//     mgop.setOrderOfAccuracy(4); // ******************** 
//     getIndex(mg.dimension(),J1,J2,J3);
//     RealArray ddxx(J1,J2,J3,Rd2), ddxy(J1,J2,J3,Rd2), ddyy(J1,J2,J3,Rd2);          // COMPUTE THESE AT MORE POINTS FOR BELOW
//     getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
//     mgop.derivative( MappedGridOperators::xxDerivative,rxLocal,ddxx,J1,J2,J3,Rd2); // WE NEED AN EXTRA GHOST HERE !!!!          
//     mgop.derivative( MappedGridOperators::xyDerivative,rxLocal,ddxy,J1,J2,J3,Rd2);            
//     mgop.derivative( MappedGridOperators::yyDerivative,rxLocal,ddyy,J1,J2,J3,Rd2);             

// //     ddxx=0.;
// //     mgop.derivative( MappedGridOperators::yDerivative,rxLocal,ddxx,J1,J2,J3,0); // WE NEED AN EXTRA GHOST HERE !!!!          

// // ::display(rxLocal(J1,2,2,0),"rxLocal","%g ");
// // ::display(ddxx(J1,2,2,0),"(rx).xx","%g ");
// // ::display(ddxx(J1,3,3,0),"(rx).xx","%g ");

//     #define DDX(i1,i2,i3,m1,m2) ddx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
//     #define DDY(i1,i2,i3,m1,m2) ddy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 

//     #define DDXX(i1,i2,i3,m1,m2) ddxx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
//     #define DDXY(i1,i2,i3,m1,m2) ddxy(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
//     #define DDYY(i1,i2,i3,m1,m2) ddyy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 




//     // Define stencils for parametric derivatives 
//     Range R5 = Range(-2,2);
//     RealArray rrrrCoeff(R5), ssssCoeff(R5);
//     rrrrCoeff = 0.;  ssssCoeff=0.;
//     const Real dr4 = dr[0]*dr[0]*dr[0]*dr[0];
//     const Real ds4 = dr[1]*dr[1]*dr[1]*dr[1];
//     // (D+D-)^2 stencil 
//     rrrrCoeff(-2) =  1./dr4;  ssssCoeff(-2) =  1./ds4;
//     rrrrCoeff(-1) = -4./dr4;  ssssCoeff(-1) = -4./ds4;
//     rrrrCoeff( 0) =  6./dr4;  ssssCoeff( 0) =  6./ds4;
//     rrrrCoeff( 1) = -4./dr4;  ssssCoeff( 1) = -4./ds4;
//     rrrrCoeff( 2) =  1./dr4;  ssssCoeff( 2) =  1./ds4;            

//     RealArray rrrCoeff(R5), sssCoeff(R5);
//     rrrCoeff=0.; sssCoeff=0.;
//     // D0(D+D-) stencil 
//     rrrCoeff(-2) = -1./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff(-2) = -1./(2.*dr[1]*dr[1]*dr[1]);
//     rrrCoeff(-1) = +2./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff(-1) = +2./(2.*dr[1]*dr[1]*dr[1]);
//     rrrCoeff( 0) =  0.;                         sssCoeff( 0) =  0.;    
//     rrrCoeff( 1) = -2./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff( 1) = -2./(2.*dr[1]*dr[1]*dr[1]);
//     rrrCoeff( 2) =  1./(2.*dr[0]*dr[0]*dr[0]);  sssCoeff( 2) =  1./(2.*dr[1]*dr[1]*dr[1]);

//     RealArray rrCoeff(R5), ssCoeff(R5);
//     rrCoeff=0.; ssCoeff=0.;
//     // D+D-
//     rrCoeff(-1) = 1./(dr[0]*dr[0]); ssCoeff(-1) = 1./(dr[1]*dr[1]); 
//     rrCoeff( 0) =-2./(dr[0]*dr[0]); ssCoeff( 0) =-2./(dr[1]*dr[1]); 
//     rrCoeff( 1) = 1./(dr[0]*dr[0]); ssCoeff( 1) = 1./(dr[1]*dr[1]); 

//     RealArray rCoeff(R5), sCoeff(R5);
//     rCoeff=0.; sCoeff=0.;
//     // Dz stencil 
//     rCoeff(-1) = -1./(2.*dr[0]); sCoeff(-1) = -1./(2.*dr[1]); 
//     rCoeff( 0) =  0.;            sCoeff( 0) =  0.;
//     rCoeff( 1) =  1./(2.*dr[0]); sCoeff( 1) =  1./(2.*dr[1]); 

//     // Identity stencil 
//     RealArray iCoeff(R5);
//     iCoeff=0.;
//     iCoeff(0)=1.; 

//     const Real cLapSq = cImp(-1,1)*(c*dt)*(c*dt)*(c*dt)*(c*dt);
//     if( numberOfDimensions==2 )
//     {      
//       RealArray ddxxx(I1,I2,I3,Rd2), ddxyy(I1,I2,I3,Rd2), ddyyy(I1,I2,I3,Rd2);

//       // Operators have no third x derivative :(
//       // Take x derivative of xx derivative 
//       for( int dir=0; dir<=Rd2.getBound(); dir++ )
//       {
//         // rxxx.x
//         ddxxx(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddxx(I1+1,I2,I3,dir) - ddxx(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,0)*( ddxx(I1,I2+1,I3,dir) - ddxx(I1,I2-1,I3,dir) )/(2.*dr[1]);

//         // rxyy.x 
//         ddxyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,0)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1]);                                    

//         // rxyy.y
//         ddyyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,1)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1]);  
                                                                                                                                    
//       }  

//       #define DDXXX(i1,i2,i3,m1,m2) ddxxx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
//       #define DDXYY(i1,i2,i3,m1,m2) ddxyy(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
//       #define DDYYY(i1,i2,i3,m1,m2) ddyyy(i1,i2,i3,(m1)+numberOfDimensions*(m2))  

                                            
//       FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the domain
//       {
//         if( maskLocal(i1,i2,i3)>0 )
//         {
//           Real rx    =    DD(i1,i2,i3,0,0);
//           Real ry    =    DD(i1,i2,i3,0,1);
//           Real rxx   =   DDX(i1,i2,i3,0,0);
//           Real rxy   =   DDY(i1,i2,i3,0,0);
//           Real ryy   =   DDY(i1,i2,i3,0,1);  // ry.y
//           Real rxxx  =  DDXX(i1,i2,i3,0,0);
//           Real rxxy  =  DDXY(i1,i2,i3,0,0);  // rx.xy
//           Real rxyy  =  DDYY(i1,i2,i3,0,0);  // rx.yy
//           Real ryyy  =  DDYY(i1,i2,i3,0,1);  // ry.yy
//           Real rxxxx = DDXXX(i1,i2,i3,0,0);  // rx.xxx
//           Real rxxyy = DDXYY(i1,i2,i3,0,1);  // ry.xyy
//           Real ryyyy = DDYYY(i1,i2,i3,0,1);  // ry.yyy   

//           Real sx    =    DD(i1,i2,i3,1,0);
//           Real sy    =    DD(i1,i2,i3,1,1);
//           Real sxx   =   DDX(i1,i2,i3,1,0);
//           Real sxy   =   DDY(i1,i2,i3,1,0);
//           Real syy   =   DDY(i1,i2,i3,1,1);  
//           Real sxxx  =  DDXX(i1,i2,i3,1,0);
//           Real sxxy  =  DDXY(i1,i2,i3,1,0);
//           Real sxyy  =  DDYY(i1,i2,i3,1,0);  // sx.yy
//           Real syyy  =  DDYY(i1,i2,i3,1,1);  
//           Real sxxxx = DDXXX(i1,i2,i3,1,0);  // sx.xxx
//           Real sxxyy = DDXYY(i1,i2,i3,1,1);  // sy.xyy
//           Real syyyy = DDYYY(i1,i2,i3,1,1);  // sy.yyy                                                       

//           // ---- COEFFICIENTS OF LAPLACIAN SQUARED  : from laplacianCoefficients.mpl ----
//           Real urrrr = pow(rx, 0.4e1) + 0.2e1 * ry * ry * rx * rx + pow(ry, 0.4e1);
//           Real urrrs = 0.4e1 * pow(rx, 0.3e1) * sx + 0.4e1 * pow(ry, 0.3e1) * sy + 0.2e1 * ry * (sy * rx * rx + 0.2e1 * ry * sx * rx) + 0.2e1 * sy * ry * rx * rx;
//           Real urrss = 6 * rx * rx * sx * sx + 6 * ry * ry * sy * sy + 2 * ry * (2 * sy * sx * rx + ry * sx * sx) + 2 * sy * (sy * rx * rx + 2 * ry * sx * rx);
//           Real ursss = 0.4e1 * rx * pow(sx, 0.3e1) + 0.4e1 * ry * pow(sy, 0.3e1) + 0.2e1 * ry * sy * sx * sx + 0.2e1 * sy * (0.2e1 * sy * sx * rx + ry * sx * sx);
//           Real ussss = pow(sx, 0.4e1) + 0.2e1 * sy * sy * sx * sx + pow(sy, 0.4e1);
//           Real urrr = 6 * rx * rx * rxx + 6 * ry * ry * ryy + 2 * ry * (2 * rxy * rx + ry * rxx) + 2 * ryy * rx * rx + 4 * ry * rxy * rx;
//           Real urrs = ry * (3 * syy * ry + 3 * sy * ryy) + 7 * sy * ry * ryy + ry * (2 * syy * ry + 2 * sy * ryy) + syy * ry * ry + 2 * ry * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * sy * (2 * rxy * rx + ry * rxx) + 4 * ryy * sx * rx + 2 * ry * (2 * sxy * rx + 2 * sx * rxy) + 2 * syy * rx * rx + 4 * sy * rxy * rx + rx * (3 * sxx * rx + 3 * sx * rxx) + 7 * sx * rx * rxx + rx * (2 * sxx * rx + 2 * sx * rxx) + sxx * rx * rx;
//           Real urss = 7 * ry * sy * syy + sy * (3 * syy * ry + 3 * sy * ryy) + ryy * sy * sy + sy * (2 * syy * ry + 2 * sy * ryy) + 2 * ry * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * ryy * sx * sx + 4 * ry * sx * sxy + 4 * syy * sx * rx + 2 * sy * (2 * sxy * rx + 2 * sx * rxy) + 7 * rx * sx * sxx + sx * (3 * sxx * rx + 3 * sx * rxx) + rxx * sx * sx + sx * (2 * sxx * rx + 2 * sx * rxx);
//           Real usss = 6 * sx * sx * sxx + 6 * sy * sy * syy + 2 * sy * (2 * sx * sxy + sxx * sy) + 2 * syy * sx * sx + 4 * sy * sx * sxy;
//           Real urr = 4 * rx * rxxx + 4 * rxyy * rx + 3 * rxx * rxx + 2 * ryy * rxx + 4 * ry * rxxy + 4 * rxy * rxy + 4 * ry * ryyy + 3 * ryy * ryy;
//           Real urs = 4 * sxxx * rx + 4 * sxyy * rx + 6 * sxx * rxx + 2 * syy * rxx + 4 * sx * rxxx + 4 * sy * rxxy + 8 * sxy * rxy + 4 * sx * rxyy + 4 * ry * sxxy + 4 * syyy * ry + 2 * ryy * sxx + 6 * syy * ryy + 4 * sy * ryyy;
//           Real uss = 4 * sx * sxxx + 4 * sx * sxyy + 3 * sxx * sxx + 2 * sxx * syy + 4 * sxxy * sy + 4 * sxy * sxy + 4 * sy * syyy + 3 * syy * syy;
//           Real ur = rxxxx + 2 * rxxyy + ryyyy;
//           Real us = sxxxx + 2 * sxxyy + syyyy;

//           // printF(" (i1,i2)=(%d,%d) urrrr=%e, urrrs=%e, urrss=%e, ursss=%e, ussss=%e, urrr=%e urrs=%e urss=%e urr=%e urs=%e uss=%e ur=%e us=%e\n",
//           //     i1,i2,urrrr,urrrs,urrss,ursss,ussss,urrr,urrs,urss,urss,usss,urr,urs,uss,ur,us);

//           // printF(" (i1,i2)=(%d,%d) ur=%10.2e us=%10.2e rxxxx=%10.2e sxxxx=%10.2e ryyyy=%10.2e syyyy=%10.2e\n",i1,i2,ur,us,rxxxx,sxxxx,ryyyy,syyyy);

//           ForStencil(m1,m2,m3)
//           {
//             int m  = M123(m1,m2,m3);
//             // if( i1==5 && i2==2 )
//             // {
//             //    printF("Before O4 correction (i1,i2)=(%d,%d) coeff(%d)=%10.2e, ur=%10.2e\n",i1,i2,m,coeffLocal(m,i1,i2,i3),ur); 
//             // }                        
//             coeffLocal(m,i1,i2,i3) += 
//                             cLapSq*(  urrrr*rrrrCoeff(m1)*   iCoeff(m2) 
//                                     + urrrs* rrrCoeff(m1)*   sCoeff(m2)
//                                     + urrss*  rrCoeff(m1)*  ssCoeff(m2)
//                                     + ursss*   rCoeff(m1)* sssCoeff(m2)
//                                     + ussss*   iCoeff(m1)*ssssCoeff(m2)
//                                     + urrr * rrrCoeff(m1)*   iCoeff(m2)
//                                     + urrs *  rrCoeff(m1)*   sCoeff(m2)
//                                     + urss *   rCoeff(m1)*  ssCoeff(m2)
//                                     + usss *   iCoeff(m1)* sssCoeff(m2)
//                                     + urr  *  rrCoeff(m1)*   iCoeff(m2)
//                                     + urs  *   rCoeff(m1)*   sCoeff(m2)
//                                     + uss  *   iCoeff(m1)*  ssCoeff(m2)
//                                     + ur   *   rCoeff(m1)*   iCoeff(m2)
//                                     + us   *   iCoeff(m1)*   sCoeff(m2)
//                                     );

//             // save for testing : 
//             lapSqCoeffLocal(m,i1,i2,i3) = 
//                                   (  urrrr*rrrrCoeff(m1)*   iCoeff(m2) 
//                                     + urrrs* rrrCoeff(m1)*   sCoeff(m2)
//                                     + urrss*  rrCoeff(m1)*  ssCoeff(m2)
//                                     + ursss*   rCoeff(m1)* sssCoeff(m2)
//                                     + ussss*   iCoeff(m1)*ssssCoeff(m2)
//                                     + urrr * rrrCoeff(m1)*   iCoeff(m2)
//                                     + urrs *  rrCoeff(m1)*   sCoeff(m2)
//                                     + urss *   rCoeff(m1)*  ssCoeff(m2)
//                                     + usss *   iCoeff(m1)* sssCoeff(m2)
//                                     + urr  *  rrCoeff(m1)*   iCoeff(m2)
//                                     + urs  *   rCoeff(m1)*   sCoeff(m2)
//                                     + uss  *   iCoeff(m1)*  ssCoeff(m2)
//                                     + ur   *   rCoeff(m1)*   iCoeff(m2)
//                                     + us   *   iCoeff(m1)*   sCoeff(m2)
//                                     );            


//             // if( i1==5 && i2==2 )
//             // {
//             //    printF(" After O4 correction(i1,i2)=(%d,%d) coeff(%d)=%10.2e\n",i1,i2,m,coeffLocal(m,i1,i2,i3));
//             //    printF(" urrrr*rrrrCoeff(m1)*   iCoeff(m2) =%10.2e\n",urrrr*rrrrCoeff(m1)*   iCoeff(m2) );
//             //    printF(" urrrs* rrrCoeff(m1)*   sCoeff(m2) =%10.2e\n",urrrs* rrrCoeff(m1)*   sCoeff(m2));
//             //    printF(" urrss*  rrCoeff(m1)*  ssCoeff(m2) =%10.2e\n",urrss*  rrCoeff(m1)*  ssCoeff(m2));
//             //    printF(" ursss*   rCoeff(m1)* sssCoeff(m2) =%10.2e\n",ursss*   rCoeff(m1)* sssCoeff(m2));
//             //    printF(" ussss*   iCoeff(m1)*ssssCoeff(m2) =%10.2e\n",ussss*   iCoeff(m1)*ssssCoeff(m2));
//             //    printF(" urrr * rrrCoeff(m1)*   iCoeff(m2) =%10.2e\n",urrr * rrrCoeff(m1)*   iCoeff(m2));
//             //    printF(" urrs *  rrCoeff(m1)*   sCoeff(m2) =%10.2e\n",urrs *  rrCoeff(m1)*   sCoeff(m2));
//             //    printF(" urss *   rCoeff(m1)*  ssCoeff(m2) =%10.2e\n",urss *   rCoeff(m1)*  ssCoeff(m2));
//             //    printF(" usss *   iCoeff(m1)* sssCoeff(m2) =%10.2e\n",usss *   iCoeff(m1)* sssCoeff(m2));
//             //    printF(" urr  *  rrCoeff(m1)*   iCoeff(m2) =%10.2e\n",urr  *  rrCoeff(m1)*   iCoeff(m2));
//             //    printF(" urs  *   rCoeff(m1)*   sCoeff(m2) =%10.2e\n",urs  *   rCoeff(m1)*   sCoeff(m2));
//             //    printF(" uss  *   iCoeff(m1)*  ssCoeff(m2) =%10.2e\n",uss  *   iCoeff(m1)*  ssCoeff(m2));
//             //    printF(" ur   *   rCoeff(m1)*   iCoeff(m2) =%10.2e\n",ur   *   rCoeff(m1)*   iCoeff(m2));
//             //    printF(" us   *   iCoeff(m1)*   sCoeff(m2) =%10.2e\n",us   *   iCoeff(m1)*   sCoeff(m2));
//             //    printf(" ur=%10.2e rxxxx=%10.2e rxxyy=%10.2e ryyyy=%10.2e rCoeff=%10.2e iCoeff=%10.2e\n",ur,rxxxx,rxxyy,ryyyy,rCoeff(m1),iCoeff(m2));                        
//             // }                                   
//           } // end for stencil

//           if( false )
//           {
//             // check stencil for Delta^2 
//             // Real kx=1., ky=1., kz=1.;
//             // #define UE(j1,j2,j3) sin(kx*xLocal(j1,j2,j3,0))*sin(ky*xLocal(j1,j2,j3,1))*sin(kz*xLocal(j1,j2,j3,2))
//             // #define LAPSQE(j1,j2,j3) SQR(kx*kx+ky*ky+kz*kz) * UE(j1,j2,j3)

//             Real lapSq = 0.;
//             ForStencil(m1,m2,m3)
//             {
//               int m  = M123(m1,m2,m3); 

//               // lapSq += lapSqCoeff(m,i1,i2,i3)*UE(i1+m1,i2+m2,i3+m3);
//               if( true )
//               {
//                 lapSq += lapSqCoeffLocal(m,i1,i2,i3)*uLocal(i1+m1,i2+m2,i3+m3);  // chain rule coeff
//               }
//               else
//               {
//                 lapSq += coeffLapSqLocal(m,i1,i2,i3)*uLocal(i1+m1,i2+m2,i3+m3);  // delta coeff
//               }

//             }
//             // Real lapSqe = LAPSQE(i1,i2,i3); 

//             Real lapSqe = lapSqu(i1,i2,i3);

//             Real err = fabs( lapSq - lapSqe )/(fabs(lapSqe) + 1.e-10 );
//             // printf("(i1,i2,i3)=(%3d,%3d,%3d) : lapSq(u)[coeff] = %10.2e, lapSq(u)[stencil] = %10.2e, relErr= %9.2e\n",i1,i2,i3,lapSq,lapSqe, err);


//             lapSqErrLocal(i1,i2,i3)=err;
//             maxErr = max(maxErr,err);

//             if( i1==2 && i2==2 )
//             {
//               Real maxErrCoeff=0.;
//               Real coeffNorm = max(fabs(lapSqCoeffLocal(M0,i1,i2,i3)));

//               printf("lapSqCoeff(m1,m2,m3) and relative errors (i1,i2,i3)=(%3d,%3d,%3d):\n",i1,i2,i3);
//               for( int m3=-halfWidth3; m3<=halfWidth3; m3++ )
//               {
//                 printf("m3=%2d\n",m3);
//                 for( int m2=-halfWidth2; m2<=halfWidth2; m2++ ) 
//                 {
//                   printf("m2=%2d, m1=[%2d,%2d] : ",m2,-halfWidth1,halfWidth1);
//                   for( int m1=-halfWidth1; m1<=halfWidth1; m1++ )
//                   {
//                     int m  = M123(m1,m2,m3); 
//                     printF(" %10.2e ",lapSqCoeffLocal(m,i1,i2,i3));
//                   }
//                   printF("   ");
//                   for( int m1=-halfWidth1; m1<=halfWidth1; m1++ )
//                   {
//                     int m  = M123(m1,m2,m3); 
//                     Real err = fabs(lapSqCoeffLocal(m,i1,i2,i3)-coeffLapSqLocal(m,i1,i2,i3))/(coeffNorm);
//                     maxErrCoeff=max(maxErrCoeff,err);
//                     // Real err = coeffLapSqLocal(m,i1,i2,i3);
//                     printF(" %10.2e ",err);
//                   }                  
                                    
//                   printF("\n");
//                 }
//               }
//               printF(">>>(i1,i2,i3)=(%3d,%3d,%3d) : maxRelErr in coefficients = %9.2e\n",i1,i2,i3,maxErrCoeff);              
//             }


//           // Real lapSqe = lapSqu(i1,i2,i3);

//           // Real err = fabs( lapSq - lapSqe )/(fabs(lapSqe) + 1.e-10 );
//           // printf("(i1,i2,i3)=(%3d,%3d,%3d) : lapSq = %10.2e, lapSqe = %10.2e, relErr= %9.2e\n",i1,i2,i3,lapSq,lapSqe, err);

//           // maxErr = max(maxErr,err);          
//           }
//         }
//       }

//     }
//     else
//     {
//       // --------------- THREE DIMENSIONS ----------------

//       RealArray rxza(I1,I2,I3,Rd2);
//       mgop.derivative( MappedGridOperators::zDerivative,rxLocal,rxza,I1,I2,I3,Rd2);            

//       #define DDZ(i1,i2,i3,m1,m2) rxza(i1,i2,i3,(m1)+numberOfDimensions*(m2))

//       const int extra=1; 
//       Index J1,J2,J3;
//       getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
//       RealArray ddxz(J1,J2,J3,Rd2), ddyz(J1,J2,J3,Rd2), ddzz(J1,J2,J3,Rd2);          // COMPUTE THESE AT MORE POINTS FOR BELOW
//       mgop.setOrderOfAccuracy(4); // ******************** 
//       mgop.derivative( MappedGridOperators::xzDerivative,rxLocal,ddxz,J1,J2,J3,Rd2);            
//       mgop.derivative( MappedGridOperators::yzDerivative,rxLocal,ddyz,J1,J2,J3,Rd2);            
//       mgop.derivative( MappedGridOperators::zzDerivative,rxLocal,ddzz,J1,J2,J3,Rd2);             

//       #define DDXZ(i1,i2,i3,m1,m2) ddxz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
//       #define DDYZ(i1,i2,i3,m1,m2) ddyz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
//       #define DDZZ(i1,i2,i3,m1,m2) ddzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  

//       // metric third derivatives 
//       RealArray ddxxx(I1,I2,I3,Rd2), ddxyy(I1,I2,I3,Rd2), ddyyy(I1,I2,I3,Rd2);
//       RealArray ddxzz(I1,I2,I3,Rd2), ddzzz(I1,I2,I3,Rd2), ddyzz(I1,I2,I3,Rd2);

//       // Operators have no third x derivative :(
//       // Take x derivative of xx derivative 
//       for( int dir=0; dir<=Rd2.getBound(); dir++ )
//       {
//         // rxxx.x
//         ddxxx(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddxx(I1+1,I2,I3,dir) - ddxx(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,0)*( ddxx(I1,I2+1,I3,dir) - ddxx(I1,I2-1,I3,dir) )/(2.*dr[1])
//                              +DD(I1,I2,I3,2,0)*( ddxx(I1,I2,I3+1,dir) - ddxx(I1,I2,I3-1,dir) )/(2.*dr[2]);

//         // rxyy = ryy.x 
//         ddxyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,0)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1])
//                              +DD(I1,I2,I3,2,0)*( ddyy(I1,I2,I3+1,dir) - ddyy(I1,I2,I3-1,dir) )/(2.*dr[2]);


//         // ryyy = ryy.y
//         ddyyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,1)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1])
//                              +DD(I1,I2,I3,2,1)*( ddyy(I1,I2,I3+1,dir) - ddyy(I1,I2,I3-1,dir) )/(2.*dr[2]);

//         // rxzz = rzz.x 
//         ddxzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,0)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
//                              +DD(I1,I2,I3,2,0)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);

//         // ryzz = rzz.y 
//         ddyzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,1)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
//                              +DD(I1,I2,I3,2,1)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);

//         // rzzz = rzz.z
//         ddzzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,2)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
//                              +DD(I1,I2,I3,1,2)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
//                              +DD(I1,I2,I3,2,2)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);

//       } 

//       #define DDXXX(i1,i2,i3,m1,m2) ddxxx(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
//       #define DDXYY(i1,i2,i3,m1,m2) ddxyy(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
//       #define DDYYY(i1,i2,i3,m1,m2) ddyyy(i1,i2,i3,(m1)+numberOfDimensions*(m2)) 
//       #define DDXZZ(i1,i2,i3,m1,m2) ddxzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
//       #define DDZZZ(i1,i2,i3,m1,m2) ddzzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))   
//       #define DDYZZ(i1,i2,i3,m1,m2) ddyzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))   

//       RealArray ttttCoeff(R5);
//       ttttCoeff = 0.; 
//       const Real dt4 = dr[2]*dr[2]*dr[2]*dr[2];
//       // (D+D-)^2 stencil 
//       ttttCoeff(-2) =  1./dt4; 
//       ttttCoeff(-1) = -4./dt4; 
//       ttttCoeff( 0) =  6./dt4; 
//       ttttCoeff( 1) = -4./dt4; 
//       ttttCoeff( 2) =  1./dt4;             

//       RealArray tttCoeff(R5);
//       tttCoeff=0.; 
//       // D0(D+D-) stencil 
//       tttCoeff(-2) = -1./(2.*dr[2]*dr[2]*dr[2]); 
//       tttCoeff(-1) = +2./(2.*dr[2]*dr[2]*dr[2]); 
//       tttCoeff( 0) =  0.;                        
//       tttCoeff( 1) = -2./(2.*dr[2]*dr[2]*dr[2]); 
//       tttCoeff( 2) =  1./(2.*dr[2]*dr[2]*dr[2]); 

//       RealArray ttCoeff(R5);
//       ttCoeff=0.;
//       // D+D-
//       ttCoeff(-1) = 1./(dr[2]*dr[2]);
//       ttCoeff( 0) =-2./(dr[2]*dr[2]);
//       ttCoeff( 1) = 1./(dr[2]*dr[2]);

//       RealArray tCoeff(R5);
//       tCoeff=0.; 
//       // Dz stencil 
//       tCoeff(-1) = -1./(2.*dr[2]);
//       tCoeff( 0) =  0.;           
//       tCoeff( 1) =  1./(2.*dr[2]); 

//       FOR_3D(i1,i2,i3,I1,I2,I3) // loop over points on the domain
//       {
//         if( maskLocal(i1,i2,i3)>0 )
//         {
//           declareMetricDerivatives3d(r,0)
//           declareMetricDerivatives3d(s,1)
//           declareMetricDerivatives3d(t,2)

//           // ---- *** COEFFICIENTS OF 3D LAPLACIAN SQUARED : from laplacianCoefficients.mpl ----

//           Real urrrr = pow(rx, 0.4e1) + 0.2e1 * ry * ry * rx * rx + 0.2e1 * rz * rz * rx * rx + pow(ry, 0.4e1) + 0.2e1 * rz * rz * ry * ry + pow(rz, 0.4e1);
//           Real urrrs = 0.4e1 * pow(rx, 0.3e1) * sx + 0.4e1 * pow(rz, 0.3e1) * sz + 0.2e1 * rz * (sz * ry * ry + 0.2e1 * rz * sy * ry) + 0.2e1 * sz * rz * ry * ry + 0.2e1 * rz * (sz * rx * rx + 0.2e1 * rz * sx * rx) + 0.2e1 * sz * rz * rx * rx + 0.4e1 * pow(ry, 0.3e1) * sy + 0.2e1 * ry * (sy * rx * rx + 0.2e1 * ry * sx * rx) + 0.2e1 * sy * ry * rx * rx;
//           Real urrss = 6 * rx * rx * sx * sx + 6 * rz * rz * sz * sz + 2 * rz * (2 * sz * sy * ry + rz * sy * sy) + 2 * sz * (sz * ry * ry + 2 * rz * sy * ry) + 2 * rz * (2 * sz * sx * rx + rz * sx * sx) + 2 * sz * (sz * rx * rx + 2 * rz * sx * rx) + 6 * ry * ry * sy * sy + 2 * ry * (2 * sy * sx * rx + ry * sx * sx) + 2 * sy * (sy * rx * rx + 2 * ry * sx * rx);
//           Real ursss = 0.4e1 * rx * pow(sx, 0.3e1) + 0.4e1 * rz * pow(sz, 0.3e1) + 0.2e1 * rz * sz * sy * sy + 0.2e1 * sz * (0.2e1 * sz * sy * ry + rz * sy * sy) + 0.2e1 * rz * sz * sx * sx + 0.2e1 * sz * (0.2e1 * sz * sx * rx + rz * sx * sx) + 0.4e1 * ry * pow(sy, 0.3e1) + 0.2e1 * ry * sy * sx * sx + 0.2e1 * sy * (0.2e1 * sy * sx * rx + ry * sx * sx);
//           Real ussss = pow(sx, 0.4e1) + 0.2e1 * sy * sy * sx * sx + 0.2e1 * sz * sz * sx * sx + pow(sy, 0.4e1) + 0.2e1 * sz * sz * sy * sy + pow(sz, 0.4e1);
//           Real urrrt = 0.4e1 * pow(rz, 0.3e1) * tz + 0.2e1 * rz * (tz * ry * ry + 0.2e1 * rz * ty * ry) + 0.2e1 * tz * rz * ry * ry + 0.2e1 * rz * (tz * rx * rx + 0.2e1 * rz * tx * rx) + 0.2e1 * tz * rz * rx * rx + 0.4e1 * pow(ry, 0.3e1) * ty + 0.2e1 * ry * (ty * rx * rx + 0.2e1 * ry * tx * rx) + 0.2e1 * ty * ry * rx * rx + 0.4e1 * pow(rx, 0.3e1) * tx;
//           Real urrtt = 6 * rz * rz * tz * tz + 2 * rz * (2 * tz * ty * ry + rz * ty * ty) + 2 * tz * (tz * ry * ry + 2 * rz * ty * ry) + 2 * rz * (2 * tz * tx * rx + rz * tx * tx) + 2 * tz * (tz * rx * rx + 2 * rz * tx * rx) + 2 * ry * (2 * ty * tx * rx + ry * tx * tx) + 2 * ty * (ty * rx * rx + 2 * ry * tx * rx) + 6 * ry * ry * ty * ty + 6 * rx * rx * tx * tx;
//           Real urttt = 0.4e1 * rz * pow(tz, 0.3e1) + 0.2e1 * rz * tz * ty * ty + 0.2e1 * tz * (0.2e1 * tz * ty * ry + rz * ty * ty) + 0.2e1 * rz * tz * tx * tx + 0.2e1 * tz * (0.2e1 * tz * tx * rx + rz * tx * tx) + 0.2e1 * ry * ty * tx * tx + 0.2e1 * ty * (0.2e1 * ty * tx * rx + ry * tx * tx) + 0.4e1 * ry * pow(ty, 0.3e1) + 0.4e1 * rx * pow(tx, 0.3e1);
//           Real utttt = pow(tx, 0.4e1) + 0.2e1 * ty * ty * tx * tx + 0.2e1 * tz * tz * tx * tx + pow(ty, 0.4e1) + 0.2e1 * tz * tz * ty * ty + pow(tz, 0.4e1);
//           Real ussst = 0.4e1 * pow(sz, 0.3e1) * tz + 0.2e1 * sz * (tz * sy * sy + 0.2e1 * sz * ty * sy) + 0.2e1 * tz * sz * sy * sy + 0.2e1 * sz * (tz * sx * sx + 0.2e1 * sz * tx * sx) + 0.2e1 * tz * sz * sx * sx + 0.2e1 * sy * (ty * sx * sx + 0.2e1 * sy * tx * sx) + 0.2e1 * ty * sy * sx * sx + 0.4e1 * pow(sy, 0.3e1) * ty + 0.4e1 * pow(sx, 0.3e1) * tx;
//           Real usstt = 6 * sz * sz * tz * tz + 2 * sz * (2 * tz * ty * sy + sz * ty * ty) + 2 * tz * (tz * sy * sy + 2 * sz * ty * sy) + 2 * sz * (2 * tz * tx * sx + sz * tx * tx) + 2 * tz * (tz * sx * sx + 2 * sz * tx * sx) + 2 * sy * (2 * ty * tx * sx + sy * tx * tx) + 2 * ty * (ty * sx * sx + 2 * sy * tx * sx) + 6 * sy * sy * ty * ty + 6 * sx * sx * tx * tx;
//           Real usttt = 0.4e1 * sz * pow(tz, 0.3e1) + 0.2e1 * sz * tz * ty * ty + 0.2e1 * tz * (0.2e1 * tz * ty * sy + sz * ty * ty) + 0.2e1 * sz * tz * tx * tx + 0.2e1 * tz * (0.2e1 * tz * tx * sx + sz * tx * tx) + 0.2e1 * sy * ty * tx * tx + 0.2e1 * ty * (0.2e1 * ty * tx * sx + sy * tx * tx) + 0.4e1 * sy * pow(ty, 0.3e1) + 0.4e1 * sx * pow(tx, 0.3e1);
                    
//           Real urrst = 2 * rz * (2 * tz * sy * ry + 2 * sz * ty * ry + 2 * rz * ty * sy) + 2 * sz * (tz * ry * ry + 2 * rz * ty * ry) + 2 * tz * (sz * ry * ry + 2 * rz * sy * ry) + 12 * rz * rz * tz * sz + 2 * rz * (2 * tz * sx * rx + 2 * sz * tx * rx + 2 * rz * tx * sx) + 2 * sz * (tz * rx * rx + 2 * rz * tx * rx) + 2 * tz * (sz * rx * rx + 2 * rz * sx * rx) + 12 * ry * ry * ty * sy + 2 * ry * (2 * ty * sx * rx + 2 * sy * tx * rx + 2 * ry * tx * sx) + 2 * sy * (ty * rx * rx + 2 * ry * tx * rx) + 2 * ty * (sy * rx * rx + 2 * ry * sx * rx) + 12 * rx * rx * tx * sx;
//           Real ursst = 2 * rz * (tz * sy * sy + 2 * sz * ty * sy) + 2 * sz * (2 * tz * sy * ry + 2 * sz * ty * ry + 2 * rz * ty * sy) + 2 * tz * (2 * sz * sy * ry + rz * sy * sy) + 12 * rz * sz * sz * tz + 2 * rz * (tz * sx * sx + 2 * sz * tx * sx) + 2 * sz * (2 * tz * sx * rx + 2 * sz * tx * rx + 2 * rz * tx * sx) + 2 * tz * (2 * sz * sx * rx + rz * sx * sx) + 12 * ry * sy * sy * ty + 2 * ry * (ty * sx * sx + 2 * sy * tx * sx) + 2 * sy * (2 * ty * sx * rx + 2 * sy * tx * rx + 2 * ry * tx * sx) + 2 * ty * (2 * sy * sx * rx + ry * sx * sx) + 12 * rx * sx * sx * tx;
//           Real urstt = 2 * rz * (2 * tz * ty * sy + sz * ty * ty) + 2 * sz * (2 * tz * ty * ry + rz * ty * ty) + 2 * tz * (2 * tz * sy * ry + 2 * sz * ty * ry + 2 * rz * ty * sy) + 12 * rz * sz * tz * tz + 2 * rz * (2 * tz * tx * sx + sz * tx * tx) + 2 * sz * (2 * tz * tx * rx + rz * tx * tx) + 2 * tz * (2 * tz * sx * rx + 2 * sz * tx * rx + 2 * rz * tx * sx) + 12 * ry * sy * ty * ty + 2 * ry * (2 * ty * tx * sx + sy * tx * tx) + 2 * sy * (2 * ty * tx * rx + ry * tx * tx) + 2 * ty * (2 * ty * sx * rx + 2 * sy * tx * rx + 2 * ry * tx * sx) + 12 * rx * sx * tx * tx;

//           Real urrr = 6. * rx * rx * rxx + 6. * rz * rz * rzz + 2. * rz * (2. * ryz * ry + rz * ryy) + 2. * rzz * ry * ry + 4. * rz * ryz * ry + 2. * rz * (2. * rxz * rx + rz * rxx) + 2. * rzz * rx * rx + 4. * rz * rxz * rx + 6. * ry * ry * ryy + 2. * ry * (2. * rxy * rx + ry * rxx) + 2. * ryy * rx * rx + 4. * ry * rxy * rx;
//           Real urrs = 7 * sx * rx * rxx + rz * (3 * rz * szz + 3 * sz * rzz) + rz * (2 * rz * szz + 2 * sz * rzz) + szz * rz * rz + 7 * sz * rz * rzz + 2 * rz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * sz * (2 * ryz * ry + rz * ryy) + 2 * rz * (2 * ry * syz + 2 * ryz * sy) + 2 * szz * ry * ry + 4 * rzz * sy * ry + 4 * sz * ryz * ry + 2 * rz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * sz * (2 * rxz * rx + rz * rxx) + 2 * rz * (2 * rx * sxz + 2 * rxz * sx) + 2 * szz * rx * rx + 4 * rzz * sx * rx + 4 * sz * rxz * rx + ry * (3 * syy * ry + 3 * sy * ryy) + ry * (2 * syy * ry + 2 * sy * ryy) + syy * ry * ry + 7 * sy * ry * ryy + 4 * ryy * sx * rx + 4 * sy * rxy * rx + 2 * ry * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * sy * (2 * rxy * rx + ry * rxx) + 2 * ry * (2 * sxy * rx + 2 * sx * rxy) + 2 * syy * rx * rx + rx * (2 * sxx * rx + 2 * sx * rxx) + sxx * rx * rx + rx * (3 * sxx * rx + 3 * sx * rxx);
//           Real urss = 7 * rx * sx * sxx + sz * (3 * rz * szz + 3 * sz * rzz) + rzz * sz * sz + sz * (2 * rz * szz + 2 * sz * rzz) + 7 * rz * sz * szz + 2 * sz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * rzz * sy * sy + 2 * sz * (2 * ry * syz + 2 * ryz * sy) + 2 * rz * (2 * syz * sy + sz * syy) + 4 * rz * syz * sy + 4 * szz * sy * ry + 2 * rz * (2 * sxz * sx + sz * sxx) + 2 * sz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * rzz * sx * sx + 2 * sz * (2 * rx * sxz + 2 * rxz * sx) + 4 * rz * sxz * sx + 4 * szz * sx * rx + sy * (3 * syy * ry + 3 * sy * ryy) + ryy * sy * sy + sy * (2 * syy * ry + 2 * sy * ryy) + 7 * ry * sy * syy + 4 * ry * sx * sxy + 4 * syy * sx * rx + 2 * ry * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * ryy * sx * sx + 2 * sy * (2 * sxy * rx + 2 * sx * rxy) + sx * (3 * sxx * rx + 3 * sx * rxx) + rxx * sx * sx + sx * (2 * sxx * rx + 2 * sx * rxx);
//           Real usss = 6 * sx * sx * sxx + 6 * sz * sz * szz + 2 * sz * (2 * syz * sy + sz * syy) + 2 * szz * sy * sy + 4 * sz * syz * sy + 2 * sz * (2 * sxz * sx + sz * sxx) + 2 * szz * sx * sx + 4 * sz * sxz * sx + 6 * sy * sy * syy + 2 * sy * (2 * sx * sxy + sxx * sy) + 2 * syy * sx * sx + 4 * sy * sx * sxy;
//           Real urrt = rz * (3 * tzz * rz + 3 * tz * rzz) + rz * (2 * tzz * rz + 2 * tz * rzz) + tzz * rz * rz + 7 * tz * rz * rzz + 2 * rz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * tz * (2 * ryz * ry + rz * ryy) + 2 * rz * (2 * tyz * ry + 2 * ty * ryz) + 2 * tzz * ry * ry + 4 * rzz * ty * ry + 4 * tz * ryz * ry + 2 * rz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * tz * (2 * rxz * rx + rz * rxx) + 2 * rz * (2 * txz * rx + 2 * tx * rxz) + 2 * tzz * rx * rx + 4 * rzz * tx * rx + 4 * tz * rxz * rx + ry * (2 * tyy * ry + 2 * ty * ryy) + tyy * ry * ry + ry * (3 * tyy * ry + 3 * ty * ryy) + 2 * ry * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ty * (2 * rxy * rx + ry * rxx) + 2 * ry * (2 * txy * rx + 2 * tx * rxy) + 2 * tyy * rx * rx + 7 * ty * ry * ryy + 4 * ryy * tx * rx + 4 * ty * rxy * rx + rx * (3 * txx * rx + 3 * tx * rxx) + rx * (2 * txx * rx + 2 * tx * rxx) + txx * rx * rx + 7 * tx * rx * rxx;
//           Real urtt = tz * (3 * tzz * rz + 3 * tz * rzz) + rzz * tz * tz + tz * (2 * tzz * rz + 2 * tz * rzz) + 7 * rz * tz * tzz + 2 * rz * (2 * ty * tyz + tyy * tz) + 2 * tz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * rzz * ty * ty + 2 * tz * (2 * tyz * ry + 2 * ty * ryz) + 4 * rz * ty * tyz + 4 * tzz * ty * ry + 2 * rz * (2 * tx * txz + txx * tz) + 2 * tz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * rzz * tx * tx + 2 * tz * (2 * txz * rx + 2 * tx * rxz) + 4 * rz * tx * txz + 4 * tzz * tx * rx + ty * (3 * tyy * ry + 3 * ty * ryy) + ryy * ty * ty + ty * (2 * tyy * ry + 2 * ty * ryy) + 2 * ry * (2 * tx * txy + txx * ty) + 2 * ty * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ryy * tx * tx + 2 * ty * (2 * txy * rx + 2 * tx * rxy) + 7 * ry * ty * tyy + 4 * ry * tx * txy + 4 * tyy * tx * rx + tx * (3 * txx * rx + 3 * tx * rxx) + rxx * tx * tx + tx * (2 * txx * rx + 2 * tx * rxx) + 7 * rx * tx * txx;
//           Real uttt = 6 * tz * tz * tzz + 2 * tz * (2 * ty * tyz + tyy * tz) + 2 * tzz * ty * ty + 4 * tz * ty * tyz + 2 * tz * (2 * tx * txz + txx * tz) + 2 * tzz * tx * tx + 4 * tz * tx * txz + 2 * ty * (2 * tx * txy + txx * ty) + 2 * tyy * tx * tx + 4 * ty * tx * txy + 6 * ty * ty * tyy + 6 * tx * tx * txx;
//           Real usst = sz * (3 * tzz * sz + 3 * tz * szz) + sz * (2 * tzz * sz + 2 * tz * szz) + tzz * sz * sz + 7 * tz * sz * szz + 2 * sz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * tz * (2 * syz * sy + sz * syy) + 2 * sz * (2 * tyz * sy + 2 * ty * syz) + 2 * tzz * sy * sy + 4 * szz * ty * sy + 4 * tz * syz * sy + 2 * sz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * tz * (2 * sxz * sx + sz * sxx) + 2 * sz * (2 * txz * sx + 2 * tx * sxz) + 2 * tzz * sx * sx + 4 * szz * tx * sx + 4 * tz * sxz * sx + sy * (2 * tyy * sy + 2 * ty * syy) + tyy * sy * sy + sy * (3 * tyy * sy + 3 * ty * syy) + 2 * sy * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * ty * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * txy * sx + 2 * tx * sxy) + 2 * tyy * sx * sx + 7 * ty * sy * syy + 4 * syy * tx * sx + 4 * ty * sx * sxy + sx * (3 * txx * sx + 3 * tx * sxx) + sx * (2 * txx * sx + 2 * tx * sxx) + txx * sx * sx + 7 * tx * sx * sxx;
//           Real ustt = tz * (2 * tzz * sz + 2 * tz * szz) + tz * (3 * tzz * sz + 3 * tz * szz) + szz * tz * tz + 7 * sz * tz * tzz + 2 * sz * (2 * ty * tyz + tyy * tz) + 2 * tz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * szz * ty * ty + 2 * tz * (2 * tyz * sy + 2 * ty * syz) + 4 * sz * ty * tyz + 4 * tzz * ty * sy + 2 * sz * (2 * tx * txz + txx * tz) + 2 * tz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * szz * tx * tx + 2 * tz * (2 * txz * sx + 2 * tx * sxz) + 4 * sz * tx * txz + 4 * tzz * tx * sx + ty * (3 * tyy * sy + 3 * ty * syy) + syy * ty * ty + ty * (2 * tyy * sy + 2 * ty * syy) + 2 * sy * (2 * tx * txy + txx * ty) + 2 * ty * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * syy * tx * tx + 2 * ty * (2 * txy * sx + 2 * tx * sxy) + 7 * sy * ty * tyy + 4 * sy * tx * txy + 4 * tyy * tx * sx + tx * (3 * txx * sx + 3 * tx * sxx) + sxx * tx * tx + tx * (2 * txx * sx + 2 * tx * sxx) + 7 * sx * tx * txx;
                    
//           Real urst = sz * (3 * tzz * rz + 3 * tz * rzz) + tz * (3 * rz * szz + 3 * sz * rzz) + rz * (2 * tzz * sz + 2 * tz * szz) + sz * (2 * tzz * rz + 2 * tz * rzz) + tz * (2 * rz * szz + 2 * sz * rzz) + rz * (3 * tzz * sz + 3 * tz * szz) + 2 * rzz * tz * sz + 2 * szz * tz * rz + 2 * tzz * rz * sz + 2 * rz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * sz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * tz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * rz * (2 * tyz * sy + 2 * ty * syz) + 2 * sz * (2 * tyz * ry + 2 * ty * ryz) + 2 * tz * (2 * ry * syz + 2 * ryz * sy) + 4 * rzz * ty * sy + 4 * szz * ty * ry + 4 * tzz * sy * ry + 2 * rz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * sz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * tz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * rz * (2 * txz * sx + 2 * tx * sxz) + 2 * sz * (2 * txz * rx + 2 * tx * rxz) + 2 * tz * (2 * rx * sxz + 2 * rxz * sx) + 4 * rzz * tx * sx + 4 * szz * tx * rx + 4 * tzz * sx * rx + ty * (2 * syy * ry + 2 * sy * ryy) + ry * (3 * tyy * sy + 3 * ty * syy) + sy * (3 * tyy * ry + 3 * ty * ryy) + ty * (3 * syy * ry + 3 * sy * ryy) + ry * (2 * tyy * sy + 2 * ty * syy) + sy * (2 * tyy * ry + 2 * ty * ryy) + 2 * ry * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * sy * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ty * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * ry * (2 * txy * sx + 2 * tx * sxy) + 2 * sy * (2 * txy * rx + 2 * tx * rxy) + 2 * ty * (2 * sxy * rx + 2 * sx * rxy) + 2 * ryy * ty * sy + 2 * syy * ty * ry + 2 * tyy * sy * ry + 4 * ryy * tx * sx + 4 * syy * tx * rx + 4 * tyy * sx * rx + rx * (3 * txx * sx + 3 * tx * sxx) + sx * (3 * txx * rx + 3 * tx * rxx) + tx * (3 * sxx * rx + 3 * sx * rxx) + rx * (2 * txx * sx + 2 * tx * sxx) + sx * (2 * txx * rx + 2 * tx * rxx) + tx * (2 * sxx * rx + 2 * sx * rxx) + 2 * rxx * tx * sx + 2 * sxx * tx * rx + 2 * txx * sx * rx;

//           Real urr = 4. * rx * rxxx + 4 * rxyy * rx + 4 * rxzz * rx + 3 * rxx * rxx + 2 * ryy * rxx + 2 * rzz * rxx + 4 * ry * rxxy + 4 * rz * rxxz + 4 * rxy * rxy + 4 * rxz * rxz + 4 * ry * ryyy + 4 * ryzz * ry + 3 * ryy * ryy + 2 * rzz * ryy + 4 * rz * ryyz + 4 * ryz * ryz + 4 * rz * rzzz + 3 * rzz * rzz;
//           Real urs = 4 * rz * szzz + 4 * sz * rzzz + 6 * rzz * szz + 4 * rz * syyz + 4 * sz * ryyz + 2 * rzz * syy + 2 * szz * ryy + 4 * ryzz * sy + 8 * ryz * syz + 4 * ry * syzz + 4 * rz * sxxz + 4 * sz * rxxz + 2 * rzz * sxx + 2 * szz * rxx + 4 * rxzz * sx + 8 * rxz * sxz + 4 * rx * sxzz + 4 * sy * ryyy + 6 * syy * ryy + 4 * syyy * ry + 4 * ry * sxxy + 4 * sy * rxxy + 2 * ryy * sxx + 2 * syy * rxx + 4 * sxyy * rx + 8 * sxy * rxy + 4 * sx * rxyy + 4 * sx * rxxx + 6 * sxx * rxx + 4 * sxxx * rx;
//           Real uss = 4 * sx * sxxx + 4 * sx * sxyy + 4 * sxzz * sx + 3 * sxx * sxx + 2 * sxx * syy + 2 * szz * sxx + 4 * sxxy * sy + 4 * sz * sxxz + 4 * sxy * sxy + 4 * sxz * sxz + 4 * sy * syyy + 4 * syzz * sy + 3 * syy * syy + 2 * szz * syy + 4 * sz * syyz + 4 * syz * syz + 4 * sz * szzz + 3 * szz * szz;
//           Real urt = 4 * tz * rzzz + 6 * tzz * rzz + 4 * tzzz * rz + 4 * rz * tyyz + 4 * tz * ryyz + 2 * rzz * tyy + 2 * tzz * ryy + 4 * tyzz * ry + 8 * tyz * ryz + 4 * ty * ryzz + 4 * rz * txxz + 4 * tz * rxxz + 2 * rzz * txx + 2 * tzz * rxx + 4 * txzz * rx + 8 * txz * rxz + 4 * tx * rxzz + 4 * ty * ryyy + 6 * tyy * ryy + 4 * tyyy * ry + 2 * tyy * rxx + 4 * txyy * rx + 8 * txy * rxy + 4 * tx * rxyy + 4 * ry * txxy + 4 * ty * rxxy + 2 * ryy * txx + 4 * tx * rxxx + 6 * txx * rxx + 4 * txxx * rx;
//           Real utt = 4 * tx * txxx + 4 * tx * txyy + 4 * tx * txzz + 3 * txx * txx + 2 * txx * tyy + 2 * txx * tzz + 4 * txxy * ty + 4 * txxz * tz + 4 * txy * txy + 4 * txz * txz + 4 * ty * tyyy + 4 * ty * tyzz + 3 * tyy * tyy + 2 * tyy * tzz + 4 * tyyz * tz + 4 * tyz * tyz + 4 * tz * tzzz + 3 * tzz * tzz;
//           Real ust = 4 * tz * szzz + 6 * tzz * szz + 4 * tzzz * sz + 4 * sz * tyyz + 4 * tz * syyz + 2 * szz * tyy + 2 * tzz * syy + 4 * tyzz * sy + 8 * tyz * syz + 4 * ty * syzz + 4 * sz * txxz + 4 * tz * sxxz + 2 * szz * txx + 2 * tzz * sxx + 4 * txzz * sx + 8 * txz * sxz + 4 * tx * sxzz + 4 * ty * syyy + 6 * tyy * syy + 4 * tyyy * sy + 4 * sy * txxy + 4 * ty * sxxy + 2 * syy * txx + 2 * tyy * sxx + 4 * txyy * sx + 8 * txy * sxy + 4 * tx * sxyy + 4 * tx * sxxx + 6 * txx * sxx + 4 * txxx * sx;
                    
//           Real ur = rxxxx + 2 * rxxyy + 2 * rxxzz + ryyyy + 2 * ryyzz + rzzzz;
//           Real us = sxxxx + 2 * sxxyy + 2 * sxxzz + syyyy + 2 * syyzz + szzzz;
//           Real ut = txxxx + 2 * txxyy + 2 * txxzz + tyyyy + 2 * tyyzz + tzzzz;



//           if( true && i1==2 && i2==2 && i3==2  )
//           {
//             printF(" (i1,i2,i3)=(%3d,%3d,%3d) urrrr=%g urrrs=%g, urrss=%g ursss=%g urrrt=%g urrtt=%g urttt=%g ussss=%g ussst=%g usstt=%g usttt=%g utttt=%g\n",
//                 i1,i2,i3,
//                 urrrr, urrrs, urrss, ursss, urrrt, urrtt, urttt, ussss, ussst, usstt, usttt, utttt );
//             printF(" urrr,urrs,urss,usss,urrt,uttt,usst,ustt,urst= %g %g %g %g %g %g %g %g %g\n",urrr,urrs,urss,usss,urrt,uttt,usst,ustt,urst);
//             printF(" urr,urs,uss,urt,utt,ust,ur,us,ut= %g %g %g %g %g %g %g %g %g\n",urr,urs,uss,urt,utt,ust,ur,us,ut);
//             printF(" rx,ry,rz,rxx,rxy,rxz,ryy,ryz,rzz,rxxx,rxxy,rxxz,rxyy,rxzz= %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",rx,ry,rz,rxx,rxy,rxz,ryy,ryz,rzz,rxxx,rxxy,rxxz,rxyy,rxzz);
//             printF(" sx,sy,sz,sxx,sxy,sxz,syy,syz,szz,sxxx,sxxy,sxxz,sxyy,sxzz= %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",sx,sy,sz,sxx,sxy,sxz,syy,syz,szz,sxxx,sxxy,sxxz,sxyy,sxzz);
//             printF(" tx,ty,tz,txx,txy,txz,tyy,tyz,tzz,txxx,txxy,txxz,txyy,txzz= %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",tx,ty,tz,txx,txy,txz,tyy,tyz,tzz,txxx,txxy,txxz,txyy,txzz);
//             printF(" rxxxx,rxxyy,rxxzz,ryyyy,ryyzz,rzzzz= %g %g %g %g %g %g\n",rxxxx,rxxyy,rxxzz,ryyyy,ryyzz,rzzzz);
//             printF(" sxxxx,sxxyy,sxxzz,syyyy,syyzz,szzzz= %g %g %g %g %g %g\n",sxxxx,sxxyy,sxxzz,syyyy,syyzz,szzzz);
//             printF(" txxxx,txxyy,txxzz,tyyyy,tyyzz,tzzzz= %g %g %g %g %g %g\n",txxxx,txxyy,txxzz,tyyyy,tyyzz,tzzzz);
//             printF(" rx*sx+ry*sy+rz*sz=%g, rx*tx+ry*ty+rz*tz=%g, sx*tx+sy*ty+sz*tz=%g\n",rx*sx+ry*sy+rz*sz,rx*tx+ry*ty+rz*tz,sx*tx+sy*ty+sz*tz);
//           }

//           ForStencil(m1,m2,m3)
//           {
//             int m  = M123(m1,m2,m3);   
//             coeffLocal(m,i1,i2,i3) += 
//                             cLapSq*(  urrrr*rrrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + urrrs* rrrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + urrss*  rrCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
//                                     + ursss*   rCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
//                                     + ussss*   iCoeff(m1)*ssssCoeff(m2)*   iCoeff(m3) 
//                                     + urrrt* rrrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     + urrtt*  rrCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
//                                     + urttt*   rCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
//                                     + utttt*   iCoeff(m1)*   iCoeff(m2)*ttttCoeff(m3) 
//                                     + ussst*   iCoeff(m1)* sssCoeff(m2)*   tCoeff(m3) 
//                                     + usstt*   iCoeff(m1)*  ssCoeff(m2)*  ttCoeff(m3) 
//                                     + usttt*   iCoeff(m1)*   sCoeff(m2)* tttCoeff(m3) 

//                                     + urrst*  rrCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
//                                     + ursst*   rCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
//                                     + urstt*   rCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 

//                                     + urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + urrs *  rrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + urss *   rCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
//                                     + usss *   iCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
//                                     + urrt *  rrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     + urtt *   rCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
//                                     + uttt *   iCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
//                                     + usst *   iCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
//                                     + ustt *   iCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 

//                                     + urst *   rCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 


//                                     + urr  *  rrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + urs  *   rCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + uss  *   iCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
//                                     + urt  *   rCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     + ust  *   iCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
//                                     + utt  *   iCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3)                                             

//                                     + ur   *   rCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + us   *   iCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + ut   *   iCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     );

//             // save for testing : 
//             lapSqCoeffLocal(m,i1,i2,i3) = 
//                                    (  urrrr*rrrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + urrrs* rrrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + urrss*  rrCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
//                                     + ursss*   rCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
//                                     + ussss*   iCoeff(m1)*ssssCoeff(m2)*   iCoeff(m3) 
//                                     + urrrt* rrrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     + urrtt*  rrCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
//                                     + urttt*   rCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
//                                     + utttt*   iCoeff(m1)*   iCoeff(m2)*ttttCoeff(m3) 
//                                     + ussst*   iCoeff(m1)* sssCoeff(m2)*   tCoeff(m3) 
//                                     + usstt*   iCoeff(m1)*  ssCoeff(m2)*  ttCoeff(m3) 
//                                     + usttt*   iCoeff(m1)*   sCoeff(m2)* tttCoeff(m3) 

//                                     + urrst*  rrCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
//                                     + ursst*   rCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
//                                     + urstt*   rCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 

//                                     + urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + urrs *  rrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + urss *   rCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
//                                     + usss *   iCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
//                                     + urrt *  rrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     + urtt *   rCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
//                                     + uttt *   iCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
//                                     + usst *   iCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
//                                     + ustt *   iCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 

//             // TROUBLE HERE: 
//                                     + urst *   rCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 

//                                     + urr  *  rrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + urs  *   rCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + uss  *   iCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
//                                     + urt  *   rCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     + ust  *   iCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
//                                     + utt  *   iCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3)                                             

//                                     + ur   *   rCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//                                     + us   *   iCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
//                                     + ut   *   iCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
//                                     );  
//  // *****TEST
//             // lapSqCoeff(m,i1,i2,i3) = 
//             //                        (  urrrr*rrrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//             //                         + urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//             //                         + urr  *  rrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//             //                         + ur   *   rCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
//             //                       );

//           } // end for stencil


//           if( true && i1==2 && i2==2 && i3==2  )
//           {
//             // check stencil for Delta^2 
//             // Real kx=1., ky=1., kz=1.;
//             // #define UE(j1,j2,j3) sin(kx*xLocal(j1,j2,j3,0))*sin(ky*xLocal(j1,j2,j3,1))*sin(kz*xLocal(j1,j2,j3,2))

//             // #define LAPSQE(j1,j2,j3) SQR(kx*kx+ky*ky+kz*kz) * UE(j1,j2,j3)

//             Real lapSq = 0.;
//             ForStencil(m1,m2,m3)
//             {
//               int m  = M123(m1,m2,m3); 

//               // lapSq += lapSqCoeff(m,i1,i2,i3)*UE(i1+m1,i2+m2,i3+m3);
//               lapSq += lapSqCoeffLocal(m,i1,i2,i3)*uLocal(i1+m1,i2+m2,i3+m3);

//             }
//             // Real lapSqe = LAPSQE(i1,i2,i3); 

//             Real lapSqe = lapSqu(i1,i2,i3);

//             Real err = fabs( lapSq - lapSqe )/(fabs(lapSqe) + 1.e-10 );
//             printf("(i1,i2,i3)=(%3d,%3d,%3d) : lapSq = %10.2e, lapSqe = %10.2e, relErr= %9.2e\n",i1,i2,i3,lapSq,lapSqe, err);

//             maxErr = max(maxErr,err);

//             if( i1==2 )
//             {
//               Real maxErrCoeff=0.;
//               printf("lapSqCoeff(m1,m2,m3) and relative errors (i1,i2,i3)=(%3d,%3d,%3d):\n",i1,i2,i3);
//               for( int m3=-halfWidth3; m3<=halfWidth3; m3++ )
//               {
//                 printf("m3=%2d\n",m3);
//                 for( int m2=-halfWidth2; m2<=halfWidth2; m2++ ) 
//                 {
//                   printf("m2=%2d, m1=[%2d,%2d] : ",m2,-halfWidth1,halfWidth1);
//                   for( int m1=-halfWidth1; m1<=halfWidth1; m1++ )
//                   {
//                     int m  = M123(m1,m2,m3); 
//                     printF(" %10.2e ",lapSqCoeffLocal(m,i1,i2,i3));
//                   }
//                   printF("   ");
//                   for( int m1=-halfWidth1; m1<=halfWidth1; m1++ )
//                   {
//                     int m  = M123(m1,m2,m3); 
//                     Real err = fabs(lapSqCoeffLocal(m,i1,i2,i3)-coeffLapSqLocal(m,i1,i2,i3))/(fabs(coeffLapSqLocal(m,i1,i2,i3))+1.);
//                     maxErrCoeff=max(maxErrCoeff,err);
//                     // Real err = coeffLapSqLocal(m,i1,i2,i3);
//                     printF(" %10.2e ",err);
//                   }                  
                                    
//                   printF("\n");
//                 }
//               }
//               printF(">>>(i1,i2,i3)=(%3d,%3d,%3d) : maxRelErr in coefficients = %9.2e\n",i1,i2,i3,maxErrCoeff);              
//             }
                        
//           }


//         } // end if maskLocal > 0 

//       } // end for3d

//       // OV_ABORT("Stop here for now");


//     }



//   } // end curvilinear
// #endMacro  



// ============================================================================================
// Macro to add foruth order correction terms to the implicit matrix
//  
// Add modified equation term 
//        (L_2h)^2  
// ============================================================================================

// void display(realArray & u )
// {
//   printF("u.getlength(0)=%i\n",u.getLength(0));
    
//   ::display(u,"u");
// }

// #define ogf EXTERN_C_NAME(ogf)
// #define ogderiv EXTERN_C_NAME(ogderiv)
// extern "C"
// {

//    /* Here are functions for TZ flow that can be called from fortran */

//   real
//   ogf(OGFunction *&ep, const real &x, const real &y,const real &z, const int & c, const real & t )
//   {
//     return (*ep)(x,y,z,c,t);
//   }
    
    
//   /* return a general derivative */
//   void
//   ogderiv(OGFunction *&ep, const int & ntd, const int & nxd, const int & nyd, const int & nzd, 
//            const real &x, const real &y, const real &z, const real & t, const int & n, real & ud )
//   {
//     ud=(*ep).gd(ntd,nxd,nyd,nzd,x,y,z,n,t);
//   }


// }

// #define getWaveDerivatives EXTERN_C_NAME(getwavederivatives)
// #define hierDeriv EXTERN_C_NAME(hierderiv)

// extern "C"
// {
//   void getWaveDerivatives( const int&nd, const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
//                            const int&gridIndexRange, const int&dimRange, const int&isPeriodic, 
//                            const real & u, const int&mask, const real&rsxy, const real& xy, const int&boundaryCondition, 
//                            const int&ipar, const real &rpar, real & maxErr, const int&ierr );

//   void hierDeriv( const int&nd, const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
//                   const int&gridIndexRange, const int&dimRange, const int&isPeriodic, 
//                   const real & u, const int&mask, const real&rsxy, const real& xy, 
//                   real & ud, real & ude,
//                   const int&boundaryCondition, 
//                   const int&ipar, const real &rpar, real & maxErr, real & l2Err, const int&ierr );
// }



// =================================================================================================
//   Create a twilightzone exact solution object
// =================================================================================================
int createTwilightZoneSolution( CgWave::TwilightZoneEnum & twilightZone, OGFunction *& tz,
                          int numberOfDimensions, int degreeInSpace, int degreeInTime, RealArray & trigFreq )
{

    if( twilightZone==CgWave::polynomial )
    {
        if( degreeInSpace>6 )
        {
            printF("createTZ:ERROR: degreeInSpace=%d is not supported by OGPolyFunction, Reducing to 6\n",degreeInSpace);
            degreeInSpace=6;
        }
        int numberOfComponentsForTZ=1;
        tz = new OGPolyFunction(degreeInSpace,numberOfDimensions,numberOfComponentsForTZ,degreeInTime);

        const int ndp=max(max(5,degreeInSpace+1),degreeInTime+1);
    
        printF("\n $$$$$$$ setup TZ: build OGPolyFunction: numCompTz=%i degreeSpace=%i, degreeTime=%i ndp=%i $$$$\n",
                      numberOfComponentsForTZ,degreeInSpace,degreeInTime,ndp);

        RealArray spatialCoefficientsForTZ(ndp,ndp,ndp,numberOfComponentsForTZ);  
        spatialCoefficientsForTZ=0.;
        RealArray timeCoefficientsForTZ(ndp,numberOfComponentsForTZ);      
        timeCoefficientsForTZ=0.;

        const int degreeInSpaceZ = numberOfDimensions==2 ? 0 : degreeInSpace;
        for( int iz=0; iz<=degreeInSpaceZ; iz++ )
        {
            for( int iy=0; iy<=degreeInSpace; iy++ )
            {
                for( int ix=0; ix<=degreeInSpace; ix++ )
                {
                    for( int n=0; n<numberOfComponentsForTZ; n++ )
                    {
            // coeff of x^ix * y^iy * z^iz 
                        if( ix+iy+iz <= degreeInSpace )
                        {
                            spatialCoefficientsForTZ(ix,iy,iz,n) = 1./( 1. + ix + 1.5*iy + 1.25*iz + n );
                        }
                        
                    }
                    
                }
            }
        }
        
        for( int n=0; n<numberOfComponentsForTZ; n++ )
        {
            for( int i=0; i<ndp; i++ )
                timeCoefficientsForTZ(i,n)= i<=degreeInTime ? 1./(i+1) : 0. ;
        }
        ::display(timeCoefficientsForTZ,"timeCoefficientsForTZ","%6.3f ");

        ((OGPolyFunction*)tz)->setCoefficients( spatialCoefficientsForTZ,timeCoefficientsForTZ );       

    }
    else if( twilightZone==CgWave::trigonometric )
    {

        const int numberOfComponents=1; 
        RealArray fx( numberOfComponents),fy( numberOfComponents),fz( numberOfComponents),ft( numberOfComponents);
        RealArray gx( numberOfComponents),gy( numberOfComponents),gz( numberOfComponents),gt( numberOfComponents);
        gx=0.;
        gy=0.;
        gz=0.;
        gt=0.;
        RealArray amplitude( numberOfComponents), cc( numberOfComponents);
        amplitude=1.;
        cc=0.;

    // RealArray & trigFreq = dbase.get<RealArray>("trigFreq");

    // fx= dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[0];
    // fy =  numberOfDimensions>1 ?  dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[1] : 0.;
    // fz =  numberOfDimensions>2 ?  dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[2] : 0.;
    // ft =  dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[3];

        fx=trigFreq(0); 
        fy=trigFreq(1); 
        fz=trigFreq(2); 
        ft=trigFreq(3); 

        tz = new OGTrigFunction(fx,fy,fz,ft);

        OGTrigFunction & trig = *((OGTrigFunction*)tz);  // cast tz to be an OGTrigFunction
        
        trig.setShifts(gx,gy,gz,gt);
        trig.setAmplitudes(amplitude);
        trig.setConstants(cc);

    }
    else
    {
        OV_ABORT("createTZ:ERROR: unknown twilightZone");
    }
    
    

    return 0;

}



int 
main(int argc, char *argv[])
{
    Overture::start(argc,argv);
  // Use this to avoid un-necessary communication: 
    Optimization_Manager::setForceVSG_Update(Off);
    const int myid=Communication_Manager::My_Process_Number;

  // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  // INIT_PETSC_SOLVER();

    int maxDeriv=2; // max derivative to compute
    int numResolutions=1;

    int orderOfAccuracyInSpace = 4;
    int numGhost = orderOfAccuracyInSpace/2; 

    CgWave::TwilightZoneEnum twilightZone = CgWave::polynomial;
    int degreeInSpace=4, degreeInTime=4;
    Real fx=1.; 
    RealArray trigFreq(4);
    trigFreq=fx; 

    int plotOption=true;
    bool smartRelease=false;
    bool reportMemory=false;
    bool loadBalance=false;
    int numberOfParallelGhost=2;

    aString caseName = "nonBox"; 

    aString nameOfOGFile= "nonBox1.order4.ng3"; 

    aString commandFileName="";
    if( argc > 1 )
    { // look at arguments for "noplot" or some other name
        int len=0;
        aString line;
        for( int i=1; i<argc; i++ )
        {
            line=argv[i];
            if( line=="-noplot" || line=="noplot" )
                plotOption=false;
            else if( len=line.matches("-g=") )
            {
                nameOfOGFile=line(len,line.length()-1);
        // printf("\n$$$$ node %i : use grid=[%s]\n",myid,(const char*)nameOfOGFile);
            }      
            else if( line=="-nopause" || line=="-abortOnEnd" || line=="-nodirect" ||
                              line=="-readCollective" || line=="-writeCollective" ||
                              line=="nopause" || line=="abortOnEnd" || line=="nodirect" )
                continue; // these commands are processed by getGraphicsInterface below 

            else if( len=line.matches("-tz=") )
            {
                aString tzType = line(len,line.length()-1);
                if( tzType == "poly" )
                {
                    twilightZone = CgWave::polynomial;
                }
                else if( tzType == "trig" )
                {
                    twilightZone = CgWave::trigonometric;
                }
                else
                {
                    printF("ERROR: unknown tzType=[%s]\n",(const char*)tzType);
                }
            } 
            else if( len=line.matches("-numRes=") )
            {
                sScanF(line(len,line.length()-1),"%i",&numResolutions); 
                printF("Setting numResolutions=%d\n",numResolutions);
            }

            else if( len=line.matches("-degreeInSpace=") )
            {
                sScanF(line(len,line.length()-1),"%i",&degreeInSpace); 
                printF("Setting degreeInSpace=%d\n",degreeInSpace);
            }
            else if( len=line.matches("-order=") )
            {
                sScanF(line(len,line.length()-1),"%i",&orderOfAccuracyInSpace); 
                printF("Setting orderOfAccuracyInSpace=%d\n",orderOfAccuracyInSpace);
            }  
            else if( len=line.matches("-deriv=") )
            {
                sScanF(line(len,line.length()-1),"%i",&maxDeriv); 
                printF("Setting maxDeriv=%d\n",maxDeriv);
            }              
            else if( len=line.matches("-fx=") )
            {
                sScanF(line(len,line.length()-1),"%e",&fx); 
                printF("Setting trig frequency: fx=%g\n",fx);
                trigFreq=fx; 
            }           

            else if( line=="memory" )
            {
                reportMemory=true;
                Diagnostic_Manager::setTrackArrayData(TRUE);
            }
            else if( line=="loadBalance" || line=="-loadBalance" )
            {
                loadBalance=true;
            }
            else if( len=line.matches("-numberOfParallelGhost=") )
            {
                sScanF(line(len,line.length()-1),"%i",&numberOfParallelGhost);
                if( numberOfParallelGhost<0 || numberOfParallelGhost>10 )
                {
                    printF("ERROR: numberOfParallelGhost=%i is no valid!\n",numberOfParallelGhost);
                    OV_ABORT("error");
                }
                printF("Setting numberOfParallelGhost=%i\n",numberOfParallelGhost);
            }
            else if( len=line.matches("-caseName=") )
            {
                caseName = line(len,line.length()-1);
                printF("Setting caseName=[%s]\n",(const char*)caseName);
            }      
            else if( commandFileName=="" )
            {
                commandFileName=line;    
                printF("testCoeff: reading commands from file [%s]\n",(const char*)commandFileName);
            }
            
        }
    }
    else
        printF("Usage: `testCoeff [options][file.cmd]' \n"
                        "     options:                            \n" 
                        "          -noplot:   run without graphics \n" 
                        "          -caseName: defines the exact solution \n" 
                        "          -abortOnEnd: abort if command file ends \n" 
                        "          -numberOfParallelGhost=<num> : number of parallel ghost lines \n" 
                        "          memory:   run with A++ memory tracking\n" 
                        "          release:  run with A++ smart release of memory\n"
                        "     file.cmd: read this command file \n");

    GenericGraphicsInterface & ps = *Overture::getGraphicsInterface("testCoeff",false,argc,argv);
    PlotStuffParameters psp;


  // By default start saving the command file called "testCoeff.cmd"
    aString logFile="testCoeff.cmd";
    ps.saveCommandFile(logFile);
    printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

    ps.appendToTheDefaultPrompt("testCoeff>");

  // read from a command file if given
    if( commandFileName!="" )
    {
        printF("read command file =%s\n",(const char*)commandFileName);
        ps.readCommandFile(commandFileName);
    }

    const int maxResolutions=6; 
  // const int maxDerivDim=20; 
  // RealArray maxErr(Range(1,maxDerivDim),maxResolutions);
  // maxErr=0.;
  // RealArray maxErr2(Range(1,maxDerivDim),maxResolutions);
  // maxErr2=0.;

  // RealArray l2Err(Range(1,maxDerivDim),maxResolutions);
  // l2Err=0.;
  // RealArray l2Err2(Range(1,maxDerivDim),maxResolutions);
  // l2Err2=0.;  

    aString gridNames[maxResolutions];
    if( caseName=="annulus" )
    {
          gridNames[0]= "annulus1.order10.ng6";
          gridNames[1]= "annulus2.order10.ng6";
          gridNames[2]= "annulus4.order10.ng6";
          gridNames[3]= "annulus8.order10.ng6";
          gridNames[4]="annulus16.order10.ng6";
    }
    else if( caseName=="square" )
    {
          gridNames[0]= "square16.order8";
          gridNames[1]= "square32.order8";
          gridNames[2]= "square64.order8";
          gridNames[3]="square128.order8";    
    }
    else if( caseName=="nonSquare" )
    {
          gridNames[0]= "nonSquare16.order8";
          gridNames[1]= "nonSquare32.order8";
          gridNames[2]= "nonSquare64.order8";
          gridNames[3]="nonSquare128.order8";    
    }  
    else if( caseName=="rotatedSquare" )
    {
          gridNames[0]= "rotatedSquare16.order8";
          gridNames[1]= "rotatedSquare32.order8";
          gridNames[2]= "rotatedSquare64.order8";
          gridNames[3]= "rotatedSquare128.order8";    
    }  
    else if( caseName=="washer" )
    {  // smoothed polygon "washer"
          gridNames[0]="washere1.order8";
          gridNames[1]="washere2.order8";
          gridNames[2]="washere4.order8";
          gridNames[3]="washere8.order8";    
    }     
    else if( caseName=="wiggley" )
    {  // TFI for a 2d wiggley domain 
          gridNames[0]="wiggley2.order8";
          gridNames[1]="wiggley4.order8";
          gridNames[2]="wiggley8.order8";
          gridNames[3]="wiggley16.order8";    
          gridNames[4]="wiggley32.order8";    
    }
  else if( caseName=="wiggley3d" )
    {  // TFI for a 3d wiggley domain 
          gridNames[0]="wiggley3de1.order4.ng3";
          gridNames[1]="wiggley3de2.order4.ng3";
          gridNames[2]="wiggley3de4.order4.ng3";
          gridNames[3]="wiggley3de8.order4.ng3";
          gridNames[4]="wiggley3de16.order4.ng3";    
          gridNames[5]="wiggley3de32.order4.ng3";    
    }   
    else if( caseName=="rhombus" )
    {  // TFI for a wiggley domain 
          gridNames[0]="rhombus2.order8";
          gridNames[1]="rhombus4.order8";
          gridNames[2]="rhombus8.order8";
          gridNames[3]="rhombus16.order8";    
          gridNames[4]="rhombus32.order8";    
    } 
    else if( caseName=="nonBox" )
    {  
          gridNames[0]="nonBox1.order4.ng3";
          gridNames[1]="nonBox2.order4.ng3";
          gridNames[2]="nonBox4.order4.ng3";
          gridNames[3]="nonBox8.order4.ng3";
          gridNames[4]="nonBox16.order4.ng3";    
          gridNames[5]="nonBox32.order4.ng3";    
    } 
    else if( caseName=="quarterCyl" )
    {  
          gridNames[0]="quarterCyl1.order4.ng3";
          gridNames[1]="quarterCyl2.order4.ng3";
          gridNames[2]="quarterCyl4.order4.ng3";
          gridNames[3]="quarterCyl8.order4.ng3";
          gridNames[4]="quarterCyl16.order4.ng3";    
          gridNames[5]="quarterCyl32.order4.ng3";    
    }  
    else if( caseName=="orthoSphere" )
    {  
          gridNames[0]="orthoSpheree1.order4.ng3";
          gridNames[1]="orthoSpheree2.order4.ng3";
          gridNames[2]="orthoSpheree4.order4.ng3";
          gridNames[3]="orthoSpheree8.order4.ng3";
          gridNames[4]="orthoSpheree16.order4.ng3";    
          gridNames[5]="orthoSpheree32.order4.ng3";    
    }                      
    else
    {
        OV_ABORT("unknown caseName");
    }
    


    for( int res=0; res<numResolutions; res++ )
    {
        nameOfOGFile = gridNames[res];

    // CompositeGrid cg;
    // const int maxWidthExtrapInterpNeighbours=4;  // This means we support 3rd-order extrap, (1,-3,3,-1)
    // aString nameOfGridFile="";
    // nameOfGridFile = readOrBuildTheGrid(ps, cg, loadBalance, numberOfParallelGhost,maxWidthExtrapInterpNeighbours );

   // create and read in a CompositeGrid
        #ifdef USE_PPP
      // On Parallel machines always add at least this many ghost lines on local arrays
            const int numGhost=2;
            MappedGrid::setMinimumNumberOfDistributedGhostLines(numGhost);
        #endif
        CompositeGrid cg;
        getFromADataBase(cg,nameOfOGFile,loadBalance);  

        Range all;
        realCompositeGridFunction u(cg,all,all,all);
        u.setName("u",0);
        u=0.;


    // --- Assign the solution from the twilightzone ----
        OGFunction *tz=NULL;
        int numberOfDimensions = cg.numberOfDimensions(); 

        createTwilightZoneSolution( twilightZone, tz, numberOfDimensions, degreeInSpace, degreeInTime, trigFreq );

        assert( tz!=NULL );
        OGFunction & exact = *tz;
        Index I1,I2,I3;
        Real t=0.;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            mg.update(MappedGrid::THEmask |MappedGrid::THEcenter | MappedGrid::THEvertex );

            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
            OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

            getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points    

            const int includeGhost=1;
            bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
            if( ok )
            {
                int numberOfComponents=1;
                Range C=numberOfComponents;
                int isRectangular=0;
                exact.gd( uLocal ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,C,t);
            }
        } // end for grid 

        const int orderOfAccuracy=orderOfAccuracyInSpace;


        CompositeGridOperators op(cg);
        op.setOrderOfAccuracy(orderOfAccuracy);   

        realCompositeGridFunction lapSq(cg,all,all,all);     // holds Lap_h^2( u ) 
        realCompositeGridFunction lapSqTrue(cg,all,all,all); // holds Lap^2( ue ) 
        realCompositeGridFunction lapSqErr(cg,all,all,all);  // holds relative err in Lap_^2( u ) 
        lapSq.setName("lapSq",0);
        lapSq=0.;

        lapSqTrue=0.;
        lapSq.setName("lapSqTrue",0);

        lapSqErr=0.;
        lapSq.setName("lapSqErr",0);

    // -- COMPUTE Lap^2 (u)
                
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            MappedGridOperators & mgop = op[grid];

            getIndex(cg[grid].dimension(),I1,I2,I3); 

            OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
            OV_GET_SERIAL_ARRAY(real,lapSq[grid],lapSqLocal);
            OV_GET_SERIAL_ARRAY(real,lapSqTrue[grid],lapSqTrueLocal);
      // OV_GET_SERIAL_ARRAY(real,lapSqErr[grid],lapSqErrLocal);

            getIndex(mg.dimension(),I1,I2,I3);
            RealArray lap(I1,I2,I3);

            mgop.setOrderOfAccuracy(2); 
            mgop.derivative(MappedGridOperators::laplacianOperator,uLocal,lap       ,I1,I2,I3); 

            getIndex(mg.gridIndexRange(),I1,I2,I3);

            mgop.derivative(MappedGridOperators::laplacianOperator,lap,   lapSqLocal,I1,I2,I3);     
            mgop.setOrderOfAccuracy(orderOfAccuracy);   

      // ----- compute lapSq(exact) -------
            bool isRectangular=false;
            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
            RealArray uxxxx(I1,I2,I3), uxxyy(I1,I2,I3), uyyyy(I1,I2,I3);
            exact.gd( uxxxx ,xLocal,numberOfDimensions,isRectangular,0,4,0,0,I1,I2,I3,0,t);
            exact.gd( uxxyy ,xLocal,numberOfDimensions,isRectangular,0,2,2,0,I1,I2,I3,0,t);
            exact.gd( uyyyy ,xLocal,numberOfDimensions,isRectangular,0,0,4,0,I1,I2,I3,0,t);

            lapSqTrueLocal(I1,I2,I3) = uxxxx(I1,I2,I3) + 2.*uxxyy(I1,I2,I3) + uyyyy(I1,I2,I3);
            if( numberOfDimensions==3 )
            {
                RealArray uxxzz(I1,I2,I3), uyyzz(I1,I2,I3), uzzzz(I1,I2,I3);
                exact.gd( uxxzz ,xLocal,numberOfDimensions,isRectangular,0,2,0,2,I1,I2,I3,0,t);
                exact.gd( uyyzz ,xLocal,numberOfDimensions,isRectangular,0,0,2,2,I1,I2,I3,0,t);
                exact.gd( uzzzz ,xLocal,numberOfDimensions,isRectangular,0,0,0,4,I1,I2,I3,0,t); 
                lapSqTrueLocal(I1,I2,I3) += uzzzz(I1,I2,I3) + 2.*uxxzz(I1,I2,I3) + 2.*uyyzz(I1,I2,I3);

            }

        } // end for grid
        
        const Real lapSqTrueNorm = maxNorm( lapSqTrue );
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];

            OV_GET_SERIAL_ARRAY(real,lapSq[grid],lapSqLocal);
            OV_GET_SERIAL_ARRAY(real,lapSqTrue[grid],lapSqTrueLocal);
            OV_GET_SERIAL_ARRAY(real,lapSqErr[grid],lapSqErrLocal);      

            lapSqErrLocal(I1,I2,I3) = fabs(lapSqLocal(I1,I2,I3)-lapSqTrueLocal(I1,I2,I3)) / lapSqTrueNorm;
        
        }
        const Real maxErrLapSq = maxNorm( lapSqErr );

        printF(">>> res=%d: max relative error in LapSq_{2h}(uExact) = %9.2e (computed with difference operators).\n",res,maxErrLapSq);


    // --- plot the solution and Lap^2(u) ----
        if( plotOption )
        {
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            psp.set(GI_TOP_LABEL,"Solution");
            ps.erase();
            PlotIt::contour( ps,u,psp );

            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            psp.set(GI_TOP_LABEL,"LapSq(u)");
            ps.erase();
            PlotIt::contour( ps,lapSq,psp );  

            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            psp.set(GI_TOP_LABEL,"LapSqErr - err in LapSq_h(ue)");
            ps.erase();
            PlotIt::contour( ps,lapSqErr,psp );            
        }



    // --------------------------------------
    // ------- EVALUATE COEFFICIENTS --------
    // --------------------------------------


        Real c=1.; 
        Real dt=.1;
        RealArray cImp(Range(-1,1),orderOfAccuracy); 
        cImp=0.; 
    // For accuracy the weights depend on one parameter beta2 for second-order,
    // and a second parameter beta4 for fourth-order: (See notes in research/timeStepping/implicitTaylorSchemes.pdf)
    // Full-weighting for second-order part by default: 
        Real beta2=.5, beta4=0.; 
        Real alpha2 = (1.-beta2)/2.;
        Real alpha4 = (alpha2-beta4-1./12.)/2.; 
        cImp(-1,0)=alpha2;
        cImp( 0,0)= beta2;
        cImp( 1,0)=alpha2;
        cImp(-1,1)=alpha4;
        cImp( 0,1)= beta4;
        cImp( 1,1)=alpha4;      

        const int e=0, cc=0; // equation number and component number 

        int stencilWidth = orderOfAccuracy + 1;
        int numberOfGhostLines= orderOfAccuracy/2;  
    
        int extraEntries = 1;  // we add 1 extra entry for interpolation equations

        const int baseStencilSize = pow(stencilWidth,cg.numberOfDimensions());   // number of entries in default stencil 
        const int stencilSize=int( baseStencilSize + extraEntries );             // add extra for interpolation and upwind equations

        const int numberOfComponentsForCoefficients=1;
        const int stencilDimension=stencilSize*SQR(numberOfComponentsForCoefficients);

        const int baseStencilDimension=baseStencilSize*SQR(numberOfComponentsForCoefficients);


        printF(">>>> stenciWidth=%d, stencilSize=%d, numberOfGhostLines=%d\n",stencilWidth,stencilSize,numberOfGhostLines);

      
        realCompositeGridFunction impCoeff;

        impCoeff.updateToMatchGrid(cg,stencilDimension,all,all,all); 
    // impCoeff.setIsACoefficientMatrix(true,baseStencilSize,numberOfGhostLines,numberOfComponentsForCoefficients);
        impCoeff.setIsACoefficientMatrix(true,stencilSize,numberOfGhostLines,numberOfComponentsForCoefficients);


    // Use these for indexing into coefficient matrices representing systems of equations
    // #define CE(c,e) (baseStencilSize*((c)+numberOfComponentsForCoefficients*(e)))
        #define M123(m1,m2,m3) (m1+halfWidth1+width*(m2+halfWidth2+width*(m3+halfWidth3)))
    // #define M123CE(m1,m2,m3,c,e) (M123(m1,m2,m3)+CE(c,e))    

    // Index I1,I2,I3;
        Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
        Index Jbv[3], &Jb1=Jbv[0], &Jb2=Jbv[1], &Jb3=Jbv[2];
        Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];
        int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];
        int jv[3], &j1=jv[0], &j2=jv[1], &j3=jv[2];
        int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];
        int m1,m2,m3; 

    // ===============================================
    // ========= GET COEFFICIENTS BY DELTA ===========
    // ===============================================

        Range M0=pow(stencilWidth,cg.numberOfDimensions());
        realCompositeGridFunction coeffLapSq(cg,M0,all,all,all);    

        
    // coefficientsByDelta( cg, coeffLapSq );
        Index Rv[3];

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            getIndex(cg[grid].gridIndexRange(),I1,I2,I3);
            if( true )
            {
                Rv[0]=I1; Rv[1]=I2; Rv[2]=I3;
            }
            else
            {
                Rv[0]=Range(2,2); 
                Rv[1]=Range(2,2); 
                Rv[2]= numberOfDimensions==2 ? Range(0,0) : Range(2,2);        
            }
            coefficientsByDelta( cg, coeffLapSq[grid],grid,Rv );
        }

    // --- compute errors in LapSq_h( u ) (FROM DELTA) -----
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            MappedGridOperators & mgop = op[grid];      

            OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
            OV_GET_SERIAL_ARRAY(Real,coeffLapSq[grid],coeffLapSqLocal);
            OV_GET_SERIAL_ARRAY(Real,lapSqTrue[grid],lapSqTrueLocal);
            OV_GET_SERIAL_ARRAY(Real,lapSqErr[grid],lapSqErrLocal);

      // ----- compute error in lapSq(u) -------
            realMappedGridFunction & coeff = impCoeff[grid];
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

            getIndex(mg.gridIndexRange(),I1,I2,I3);
      // int m1,m2,m3;
            
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                Real lapSq = 0.;
                ForStencil(m1,m2,m3)
                {
                    int m  = M123(m1,m2,m3); 
                    lapSq += coeffLapSqLocal(m,i1,i2,i3)*uLocal(i1+m1,i2+m2,i3+m3);  // delta coeff
                }

                Real err = fabs( lapSq - lapSqTrueLocal(i1,i2,i3) )/lapSqTrueNorm;
                lapSqErrLocal(i1,i2,i3)=err;
            }
        }
        const Real maxErrLapSqDelta = maxNorm(lapSqErr);

        printF(">>> res=%d: max relative error in LapSq_{2h}(uExact) = %9.2e (Coeff from DELTA).\n",res,maxErrLapSqDelta );

    // ::display(coeffLapSq[0](all,2,2,2),"coeffLapSq[0](all,2,2,2)");




    // realCompositeGridFunction lapSqErr(cg,all,all,all);
    // lapSqErr.setName("lapSqErr",0);
        lapSqErr=0.; // reset 

        realCompositeGridFunction lapSqCoeff(cg,M0,all,all,all); // holds coefficients ofr LapSq from chain rule formula

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

      // **** Here are coefficients by delta ****
            OV_GET_SERIAL_ARRAY(Real,coeffLapSq[grid],coeffLapSqLocal);


            OV_GET_SERIAL_ARRAY(Real,lapSqErr[grid],lapSqErrLocal);

            const int mDiag = M123(0,0,0);              // index of diagonal entry

            if( true )
            {
        // ----- this grid is advanced with IMPLICIT time-stepping ----
                OV_GET_SERIAL_ARRAY(Real,mg.center(),xLocal);  // *************** TEMP FOR TESTING *************
  
                OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rxLocal); 
                #define RXLOCAL(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2))   


        // Add fake metric derivatives
                if( false )
                {
                    getIndex(mg.dimension(),I1,I2,I3);
                    rxLocal=0.;
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        Real r1 = i1*dr[0]; 

            // RXLOCAL(i1,i2,i3,0,0) = 1. + r1;
            // RXLOCAL(i1,i2,i3,0,0) = 1. + xLocal(i1,i2,i3,0);
            // RXLOCAL(i1,i2,i3,1,1) = 1. + xLocal(i1,i2,i3,1);
            // RXLOCAL(i1,i2,i3,2,2) = 1. + xLocal(i1,i2,i3,2);
                        for( int m2=0; m2<numberOfDimensions; m2++)
                        for( int m1=0; m1<numberOfDimensions; m1++)
                        {
              // RXLOCAL(i1,i2,i3,m1,m2) = m1+m2 + .1*xLocal(i1,i2,i3,0) + .2*xLocal(i1,i2,i3,1)+ .1*xLocal(i1,i2,i3,2);
                            RXLOCAL(i1,i2,i3,m1,m2) = m1+m2 + .1*SQR(xLocal(i1,i2,i3,0)) + .2*xLocal(i1,i2,i3,1)+ .1*xLocal(i1,i2,i3,2);
                        }

                        // RXLOCAL(i1,i2,i3,0,0) = 1. + xLocal(i1,i2,i3,0);
            // RXLOCAL(i1,i2,i3,1,1) = 0.; // 1. + xLocal(i1,i2,i3,0);
            // RXLOCAL(i1,i2,i3,2,2) = 0.; // xLocal(i1,i2,i3,0);

                    }
                }

                getIndex(mg.gridIndexRange(),I1,I2,I3);
                RealArray lapCoeff(M0,I1,I2,I3);   

                mgop.assignCoefficients(MappedGridOperators::laplacianOperator,lapCoeff,I1,I2,I3,0,0); // 

        // eval Lap^2 u 
                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

                getIndex(mg.dimension(),I1,I2,I3);
                RealArray lap(I1,I2,I3), lapSqu(I1,I2,I3);  
                mgop.setOrderOfAccuracy(2); 
                mgop.derivative(MappedGridOperators::laplacianOperator,uLocal,lap   ,I1,I2,I3); 
                mgop.derivative(MappedGridOperators::laplacianOperator,lap,   lapSqu,I1,I2,I3);     
                mgop.setOrderOfAccuracy(orderOfAccuracy);     

                getIndex(mg.gridIndexRange(),I1,I2,I3);

                printF("\n ***ADD Fourth-order in time coefficients to the Implicit Matrix ****\n\n");

                Real maxErr = 0.;
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
                        printF("implicit matrix: order=4 CURVLINEAR ADD (Lap_2h)^2\n");
                        OV_GET_SERIAL_ARRAY(Real,lapSqCoeff[grid],lapSqCoeffLocal); 
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
                  // save for testing : 
                                    ForStencil(m1,m2,m3)
                                    {
                                        int m  = M123(m1,m2,m3);
                                        lapSqCoeffLocal(m,i1,i2,i3) = cp(m1,m2);
                                    }          
                                } // end if mask
                            } // end for 3d 
                        }
                        else
                        {
              // --------------- THREE DIMENSIONS ----------------
                            Index J1,J2,J3;
                            int extra=1;
                            getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
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
                  // save for testing : 
                                    ForStencil(m1,m2,m3)
                                    {
                                        int m  = M123(m1,m2,m3);
                                        lapSqCoeffLocal(m,i1,i2,i3) = cp(m1,m2,m3);
                                    }          
                                } // end if mask
                            } // end for 3d 
                /* ---
              // --------------- THREE DIMENSIONS ----------------
                            RealArray rxza(I1,I2,I3,Rd2);
                            mgop.derivative( MappedGridOperators::zDerivative,rxLocal,rxza,I1,I2,I3,Rd2);            
                            #define DDZ(i1,i2,i3,m1,m2) rxza(i1,i2,i3,(m1)+numberOfDimensions*(m2))
                            const int extra=1; 
                            Index J1,J2,J3;
                            getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
                            RealArray ddxz(J1,J2,J3,Rd2), ddyz(J1,J2,J3,Rd2), ddzz(J1,J2,J3,Rd2);          // COMPUTE THESE AT MORE POINTS FOR BELOW
                            mgop.setOrderOfAccuracy(4); // ******************** 
                            mgop.derivative( MappedGridOperators::xzDerivative,rxLocal,ddxz,J1,J2,J3,Rd2);            
                            mgop.derivative( MappedGridOperators::yzDerivative,rxLocal,ddyz,J1,J2,J3,Rd2);            
                            mgop.derivative( MappedGridOperators::zzDerivative,rxLocal,ddzz,J1,J2,J3,Rd2);             
                            #define DDXZ(i1,i2,i3,m1,m2) ddxz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                            #define DDYZ(i1,i2,i3,m1,m2) ddyz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
                            #define DDZZ(i1,i2,i3,m1,m2) ddzz(i1,i2,i3,(m1)+numberOfDimensions*(m2))  
              // metric third derivatives 
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
                // rxyy = ryy.x 
                                ddxyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,0)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,0)*( ddyy(I1,I2,I3+1,dir) - ddyy(I1,I2,I3-1,dir) )/(2.*dr[2]);
                // ryyy = ryy.y
                                ddyyy(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddyy(I1+1,I2,I3,dir) - ddyy(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,1)*( ddyy(I1,I2+1,I3,dir) - ddyy(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,1)*( ddyy(I1,I2,I3+1,dir) - ddyy(I1,I2,I3-1,dir) )/(2.*dr[2]);
                // rxzz = rzz.x 
                                ddxzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,0)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,0)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,0)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);
                // ryzz = rzz.y 
                                ddyzz(I1,I2,I3,dir) = DD(I1,I2,I3,0,1)*( ddzz(I1+1,I2,I3,dir) - ddzz(I1-1,I2,I3,dir) )/(2.*dr[0])
                                                                          +DD(I1,I2,I3,1,1)*( ddzz(I1,I2+1,I3,dir) - ddzz(I1,I2-1,I3,dir) )/(2.*dr[1])
                                                                          +DD(I1,I2,I3,2,1)*( ddzz(I1,I2,I3+1,dir) - ddzz(I1,I2,I3-1,dir) )/(2.*dr[2]);
                // rzzz = rzz.z
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
                                        Real rx     =     DD(i1,i2,i3,0,0);  // rx, sx, tx
                                        Real ry     =     DD(i1,i2,i3,0,1);  // ry, sy, ty
                                        Real rz     =     DD(i1,i2,i3,0,2);  // rz, sz, tz
                                        Real rxx    =    DDX(i1,i2,i3,0,0);  // rxx, sxx, txx
                                        Real rxy    =    DDX(i1,i2,i3,0,1);  // rxy, sxy, txy  
                                        Real rxz    =    DDX(i1,i2,i3,0,2);  // rxz, sxz, txz
                                        Real ryy    =    DDY(i1,i2,i3,0,1);  // ryy, syy, tyy
                                        Real ryz    =    DDY(i1,i2,i3,0,2);  // ryz, syz, tyz
                                        Real rzz    =    DDZ(i1,i2,i3,0,2);  // rzz, szz, tzz
                                        Real rxxx   =   DDXX(i1,i2,i3,0,0);  // rxxx, sxxx, txxx
                                        Real rxxy   =   DDXX(i1,i2,i3,0,1);  // rxxy, sxxy, txxy
                                        Real rxyy   =   DDYY(i1,i2,i3,0,0);  // rxyy, sxyy, txyy
                                        Real ryyy   =   DDYY(i1,i2,i3,0,1);  // ryyy, syyy, tyyy
                                        Real rxxz   =   DDXX(i1,i2,i3,0,2);  // rxxz, sxxz, txxz
                                        Real rxzz   =   DDZZ(i1,i2,i3,0,1);  // rxzz, sxzz, txzz
                                        Real rzzz   =   DDZZ(i1,i2,i3,0,2);  // rzzz, szzz, tzzz
                                        Real ryyz   =   DDYY(i1,i2,i3,0,2);  // ryyz, syyz, tyyz
                                        Real ryzz   =   DDZZ(i1,i2,i3,0,1);  // ryzz, syzz, tyzz
                                        Real rxxxx  =  DDXXX(i1,i2,i3,0,0);  // rxxx, sxxx, txxx
                                        Real rxxyy  =  DDXYY(i1,i2,i3,0,0);  // rxxyy, sxxyy, txxyy
                                        Real ryyyy  =  DDYYY(i1,i2,i3,0,1);  // ryyyy, syyyy, tyyyy
                                        Real rxxzz  =  DDXZZ(i1,i2,i3,0,0);  // rxxzz, sxxzz, txxzz
                                        Real rzzzz  =  DDZZZ(i1,i2,i3,0,2);  // rzzzz, szzzz, tzzzz
                                        Real ryyzz  =  DDYZZ(i1,i2,i3,0,1);  // ryyzz, syyzz, tyyzz
                                        Real sx     =     DD(i1,i2,i3,1,0);  // rx, sx, tx
                                        Real sy     =     DD(i1,i2,i3,1,1);  // ry, sy, ty
                                        Real sz     =     DD(i1,i2,i3,1,2);  // rz, sz, tz
                                        Real sxx    =    DDX(i1,i2,i3,1,0);  // rxx, sxx, txx
                                        Real sxy    =    DDX(i1,i2,i3,1,1);  // rxy, sxy, txy  
                                        Real sxz    =    DDX(i1,i2,i3,1,2);  // rxz, sxz, txz
                                        Real syy    =    DDY(i1,i2,i3,1,1);  // ryy, syy, tyy
                                        Real syz    =    DDY(i1,i2,i3,1,2);  // ryz, syz, tyz
                                        Real szz    =    DDZ(i1,i2,i3,1,2);  // rzz, szz, tzz
                                        Real sxxx   =   DDXX(i1,i2,i3,1,0);  // rxxx, sxxx, txxx
                                        Real sxxy   =   DDXX(i1,i2,i3,1,1);  // rxxy, sxxy, txxy
                                        Real sxyy   =   DDYY(i1,i2,i3,1,0);  // rxyy, sxyy, txyy
                                        Real syyy   =   DDYY(i1,i2,i3,1,1);  // ryyy, syyy, tyyy
                                        Real sxxz   =   DDXX(i1,i2,i3,1,2);  // rxxz, sxxz, txxz
                                        Real sxzz   =   DDZZ(i1,i2,i3,1,1);  // rxzz, sxzz, txzz
                                        Real szzz   =   DDZZ(i1,i2,i3,1,2);  // rzzz, szzz, tzzz
                                        Real syyz   =   DDYY(i1,i2,i3,1,2);  // ryyz, syyz, tyyz
                                        Real syzz   =   DDZZ(i1,i2,i3,1,1);  // ryzz, syzz, tyzz
                                        Real sxxxx  =  DDXXX(i1,i2,i3,1,0);  // rxxx, sxxx, txxx
                                        Real sxxyy  =  DDXYY(i1,i2,i3,1,0);  // rxxyy, sxxyy, txxyy
                                        Real syyyy  =  DDYYY(i1,i2,i3,1,1);  // ryyyy, syyyy, tyyyy
                                        Real sxxzz  =  DDXZZ(i1,i2,i3,1,0);  // rxxzz, sxxzz, txxzz
                                        Real szzzz  =  DDZZZ(i1,i2,i3,1,2);  // rzzzz, szzzz, tzzzz
                                        Real syyzz  =  DDYZZ(i1,i2,i3,1,1);  // ryyzz, syyzz, tyyzz
                                        Real tx     =     DD(i1,i2,i3,2,0);  // rx, sx, tx
                                        Real ty     =     DD(i1,i2,i3,2,1);  // ry, sy, ty
                                        Real tz     =     DD(i1,i2,i3,2,2);  // rz, sz, tz
                                        Real txx    =    DDX(i1,i2,i3,2,0);  // rxx, sxx, txx
                                        Real txy    =    DDX(i1,i2,i3,2,1);  // rxy, sxy, txy  
                                        Real txz    =    DDX(i1,i2,i3,2,2);  // rxz, sxz, txz
                                        Real tyy    =    DDY(i1,i2,i3,2,1);  // ryy, syy, tyy
                                        Real tyz    =    DDY(i1,i2,i3,2,2);  // ryz, syz, tyz
                                        Real tzz    =    DDZ(i1,i2,i3,2,2);  // rzz, szz, tzz
                                        Real txxx   =   DDXX(i1,i2,i3,2,0);  // rxxx, sxxx, txxx
                                        Real txxy   =   DDXX(i1,i2,i3,2,1);  // rxxy, sxxy, txxy
                                        Real txyy   =   DDYY(i1,i2,i3,2,0);  // rxyy, sxyy, txyy
                                        Real tyyy   =   DDYY(i1,i2,i3,2,1);  // ryyy, syyy, tyyy
                                        Real txxz   =   DDXX(i1,i2,i3,2,2);  // rxxz, sxxz, txxz
                                        Real txzz   =   DDZZ(i1,i2,i3,2,1);  // rxzz, sxzz, txzz
                                        Real tzzz   =   DDZZ(i1,i2,i3,2,2);  // rzzz, szzz, tzzz
                                        Real tyyz   =   DDYY(i1,i2,i3,2,2);  // ryyz, syyz, tyyz
                                        Real tyzz   =   DDZZ(i1,i2,i3,2,1);  // ryzz, syzz, tyzz
                                        Real txxxx  =  DDXXX(i1,i2,i3,2,0);  // rxxx, sxxx, txxx
                                        Real txxyy  =  DDXYY(i1,i2,i3,2,0);  // rxxyy, sxxyy, txxyy
                                        Real tyyyy  =  DDYYY(i1,i2,i3,2,1);  // ryyyy, syyyy, tyyyy
                                        Real txxzz  =  DDXZZ(i1,i2,i3,2,0);  // rxxzz, sxxzz, txxzz
                                        Real tzzzz  =  DDZZZ(i1,i2,i3,2,2);  // rzzzz, szzzz, tzzzz
                                        Real tyyzz  =  DDYZZ(i1,i2,i3,2,1);  // ryyzz, syyzz, tyyzz
                  // ---- *** COEFFICIENTS OF 3D LAPLACIAN SQUARED : from laplacianCoefficients.mpl ----
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
                                    Real urrst = 2 * rz * (2 * tz * sy * ry + 2 * sz * ty * ry + 2 * rz * ty * sy) + 2 * sz * (tz * ry * ry + 2 * rz * ty * ry) + 2 * tz * (sz * ry * ry + 2 * rz * sy * ry) + 12 * rz * rz * tz * sz + 2 * rz * (2 * tz * sx * rx + 2 * sz * tx * rx + 2 * rz * tx * sx) + 2 * sz * (tz * rx * rx + 2 * rz * tx * rx) + 2 * tz * (sz * rx * rx + 2 * rz * sx * rx) + 12 * ry * ry * ty * sy + 2 * ry * (2 * ty * sx * rx + 2 * sy * tx * rx + 2 * ry * tx * sx) + 2 * sy * (ty * rx * rx + 2 * ry * tx * rx) + 2 * ty * (sy * rx * rx + 2 * ry * sx * rx) + 12 * rx * rx * tx * sx;
                                    Real ursst = 2 * rz * (tz * sy * sy + 2 * sz * ty * sy) + 2 * sz * (2 * tz * sy * ry + 2 * sz * ty * ry + 2 * rz * ty * sy) + 2 * tz * (2 * sz * sy * ry + rz * sy * sy) + 12 * rz * sz * sz * tz + 2 * rz * (tz * sx * sx + 2 * sz * tx * sx) + 2 * sz * (2 * tz * sx * rx + 2 * sz * tx * rx + 2 * rz * tx * sx) + 2 * tz * (2 * sz * sx * rx + rz * sx * sx) + 12 * ry * sy * sy * ty + 2 * ry * (ty * sx * sx + 2 * sy * tx * sx) + 2 * sy * (2 * ty * sx * rx + 2 * sy * tx * rx + 2 * ry * tx * sx) + 2 * ty * (2 * sy * sx * rx + ry * sx * sx) + 12 * rx * sx * sx * tx;
                                    Real urstt = 2 * rz * (2 * tz * ty * sy + sz * ty * ty) + 2 * sz * (2 * tz * ty * ry + rz * ty * ty) + 2 * tz * (2 * tz * sy * ry + 2 * sz * ty * ry + 2 * rz * ty * sy) + 12 * rz * sz * tz * tz + 2 * rz * (2 * tz * tx * sx + sz * tx * tx) + 2 * sz * (2 * tz * tx * rx + rz * tx * tx) + 2 * tz * (2 * tz * sx * rx + 2 * sz * tx * rx + 2 * rz * tx * sx) + 12 * ry * sy * ty * ty + 2 * ry * (2 * ty * tx * sx + sy * tx * tx) + 2 * sy * (2 * ty * tx * rx + ry * tx * tx) + 2 * ty * (2 * ty * sx * rx + 2 * sy * tx * rx + 2 * ry * tx * sx) + 12 * rx * sx * tx * tx;
                                    Real urrr = 6. * rx * rx * rxx + 6. * rz * rz * rzz + 2. * rz * (2. * ryz * ry + rz * ryy) + 2. * rzz * ry * ry + 4. * rz * ryz * ry + 2. * rz * (2. * rxz * rx + rz * rxx) + 2. * rzz * rx * rx + 4. * rz * rxz * rx + 6. * ry * ry * ryy + 2. * ry * (2. * rxy * rx + ry * rxx) + 2. * ryy * rx * rx + 4. * ry * rxy * rx;
                                    Real urrs = 7 * sx * rx * rxx + rz * (3 * rz * szz + 3 * sz * rzz) + rz * (2 * rz * szz + 2 * sz * rzz) + szz * rz * rz + 7 * sz * rz * rzz + 2 * rz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * sz * (2 * ryz * ry + rz * ryy) + 2 * rz * (2 * ry * syz + 2 * ryz * sy) + 2 * szz * ry * ry + 4 * rzz * sy * ry + 4 * sz * ryz * ry + 2 * rz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * sz * (2 * rxz * rx + rz * rxx) + 2 * rz * (2 * rx * sxz + 2 * rxz * sx) + 2 * szz * rx * rx + 4 * rzz * sx * rx + 4 * sz * rxz * rx + ry * (3 * syy * ry + 3 * sy * ryy) + ry * (2 * syy * ry + 2 * sy * ryy) + syy * ry * ry + 7 * sy * ry * ryy + 4 * ryy * sx * rx + 4 * sy * rxy * rx + 2 * ry * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * sy * (2 * rxy * rx + ry * rxx) + 2 * ry * (2 * sxy * rx + 2 * sx * rxy) + 2 * syy * rx * rx + rx * (2 * sxx * rx + 2 * sx * rxx) + sxx * rx * rx + rx * (3 * sxx * rx + 3 * sx * rxx);
                                    Real urss = 7 * rx * sx * sxx + sz * (3 * rz * szz + 3 * sz * rzz) + rzz * sz * sz + sz * (2 * rz * szz + 2 * sz * rzz) + 7 * rz * sz * szz + 2 * sz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * rzz * sy * sy + 2 * sz * (2 * ry * syz + 2 * ryz * sy) + 2 * rz * (2 * syz * sy + sz * syy) + 4 * rz * syz * sy + 4 * szz * sy * ry + 2 * rz * (2 * sxz * sx + sz * sxx) + 2 * sz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * rzz * sx * sx + 2 * sz * (2 * rx * sxz + 2 * rxz * sx) + 4 * rz * sxz * sx + 4 * szz * sx * rx + sy * (3 * syy * ry + 3 * sy * ryy) + ryy * sy * sy + sy * (2 * syy * ry + 2 * sy * ryy) + 7 * ry * sy * syy + 4 * ry * sx * sxy + 4 * syy * sx * rx + 2 * ry * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * ryy * sx * sx + 2 * sy * (2 * sxy * rx + 2 * sx * rxy) + sx * (3 * sxx * rx + 3 * sx * rxx) + rxx * sx * sx + sx * (2 * sxx * rx + 2 * sx * rxx);
                                    Real usss = 6 * sx * sx * sxx + 6 * sz * sz * szz + 2 * sz * (2 * syz * sy + sz * syy) + 2 * szz * sy * sy + 4 * sz * syz * sy + 2 * sz * (2 * sxz * sx + sz * sxx) + 2 * szz * sx * sx + 4 * sz * sxz * sx + 6 * sy * sy * syy + 2 * sy * (2 * sx * sxy + sxx * sy) + 2 * syy * sx * sx + 4 * sy * sx * sxy;
                                    Real urrt = rz * (3 * tzz * rz + 3 * tz * rzz) + rz * (2 * tzz * rz + 2 * tz * rzz) + tzz * rz * rz + 7 * tz * rz * rzz + 2 * rz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * tz * (2 * ryz * ry + rz * ryy) + 2 * rz * (2 * tyz * ry + 2 * ty * ryz) + 2 * tzz * ry * ry + 4 * rzz * ty * ry + 4 * tz * ryz * ry + 2 * rz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * tz * (2 * rxz * rx + rz * rxx) + 2 * rz * (2 * txz * rx + 2 * tx * rxz) + 2 * tzz * rx * rx + 4 * rzz * tx * rx + 4 * tz * rxz * rx + ry * (2 * tyy * ry + 2 * ty * ryy) + tyy * ry * ry + ry * (3 * tyy * ry + 3 * ty * ryy) + 2 * ry * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ty * (2 * rxy * rx + ry * rxx) + 2 * ry * (2 * txy * rx + 2 * tx * rxy) + 2 * tyy * rx * rx + 7 * ty * ry * ryy + 4 * ryy * tx * rx + 4 * ty * rxy * rx + rx * (3 * txx * rx + 3 * tx * rxx) + rx * (2 * txx * rx + 2 * tx * rxx) + txx * rx * rx + 7 * tx * rx * rxx;
                                    Real urtt = tz * (3 * tzz * rz + 3 * tz * rzz) + rzz * tz * tz + tz * (2 * tzz * rz + 2 * tz * rzz) + 7 * rz * tz * tzz + 2 * rz * (2 * ty * tyz + tyy * tz) + 2 * tz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * rzz * ty * ty + 2 * tz * (2 * tyz * ry + 2 * ty * ryz) + 4 * rz * ty * tyz + 4 * tzz * ty * ry + 2 * rz * (2 * tx * txz + txx * tz) + 2 * tz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * rzz * tx * tx + 2 * tz * (2 * txz * rx + 2 * tx * rxz) + 4 * rz * tx * txz + 4 * tzz * tx * rx + ty * (3 * tyy * ry + 3 * ty * ryy) + ryy * ty * ty + ty * (2 * tyy * ry + 2 * ty * ryy) + 2 * ry * (2 * tx * txy + txx * ty) + 2 * ty * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ryy * tx * tx + 2 * ty * (2 * txy * rx + 2 * tx * rxy) + 7 * ry * ty * tyy + 4 * ry * tx * txy + 4 * tyy * tx * rx + tx * (3 * txx * rx + 3 * tx * rxx) + rxx * tx * tx + tx * (2 * txx * rx + 2 * tx * rxx) + 7 * rx * tx * txx;
                                    Real uttt = 6 * tz * tz * tzz + 2 * tz * (2 * ty * tyz + tyy * tz) + 2 * tzz * ty * ty + 4 * tz * ty * tyz + 2 * tz * (2 * tx * txz + txx * tz) + 2 * tzz * tx * tx + 4 * tz * tx * txz + 2 * ty * (2 * tx * txy + txx * ty) + 2 * tyy * tx * tx + 4 * ty * tx * txy + 6 * ty * ty * tyy + 6 * tx * tx * txx;
                                    Real usst = sz * (3 * tzz * sz + 3 * tz * szz) + sz * (2 * tzz * sz + 2 * tz * szz) + tzz * sz * sz + 7 * tz * sz * szz + 2 * sz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * tz * (2 * syz * sy + sz * syy) + 2 * sz * (2 * tyz * sy + 2 * ty * syz) + 2 * tzz * sy * sy + 4 * szz * ty * sy + 4 * tz * syz * sy + 2 * sz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * tz * (2 * sxz * sx + sz * sxx) + 2 * sz * (2 * txz * sx + 2 * tx * sxz) + 2 * tzz * sx * sx + 4 * szz * tx * sx + 4 * tz * sxz * sx + sy * (2 * tyy * sy + 2 * ty * syy) + tyy * sy * sy + sy * (3 * tyy * sy + 3 * ty * syy) + 2 * sy * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * ty * (2 * sx * sxy + sxx * sy) + 2 * sy * (2 * txy * sx + 2 * tx * sxy) + 2 * tyy * sx * sx + 7 * ty * sy * syy + 4 * syy * tx * sx + 4 * ty * sx * sxy + sx * (3 * txx * sx + 3 * tx * sxx) + sx * (2 * txx * sx + 2 * tx * sxx) + txx * sx * sx + 7 * tx * sx * sxx;
                                    Real ustt = tz * (2 * tzz * sz + 2 * tz * szz) + tz * (3 * tzz * sz + 3 * tz * szz) + szz * tz * tz + 7 * sz * tz * tzz + 2 * sz * (2 * ty * tyz + tyy * tz) + 2 * tz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * szz * ty * ty + 2 * tz * (2 * tyz * sy + 2 * ty * syz) + 4 * sz * ty * tyz + 4 * tzz * ty * sy + 2 * sz * (2 * tx * txz + txx * tz) + 2 * tz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * szz * tx * tx + 2 * tz * (2 * txz * sx + 2 * tx * sxz) + 4 * sz * tx * txz + 4 * tzz * tx * sx + ty * (3 * tyy * sy + 3 * ty * syy) + syy * ty * ty + ty * (2 * tyy * sy + 2 * ty * syy) + 2 * sy * (2 * tx * txy + txx * ty) + 2 * ty * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * syy * tx * tx + 2 * ty * (2 * txy * sx + 2 * tx * sxy) + 7 * sy * ty * tyy + 4 * sy * tx * txy + 4 * tyy * tx * sx + tx * (3 * txx * sx + 3 * tx * sxx) + sxx * tx * tx + tx * (2 * txx * sx + 2 * tx * sxx) + 7 * sx * tx * txx;
                                    Real urst = sz * (3 * tzz * rz + 3 * tz * rzz) + tz * (3 * rz * szz + 3 * sz * rzz) + rz * (2 * tzz * sz + 2 * tz * szz) + sz * (2 * tzz * rz + 2 * tz * rzz) + tz * (2 * rz * szz + 2 * sz * rzz) + rz * (3 * tzz * sz + 3 * tz * szz) + 2 * rzz * tz * sz + 2 * szz * tz * rz + 2 * tzz * rz * sz + 2 * rz * (2 * tyz * sy + tz * syy + 2 * ty * syz + sz * tyy) + 2 * sz * (2 * tyz * ry + tz * ryy + 2 * ty * ryz + rz * tyy) + 2 * tz * (2 * ry * syz + sz * ryy + 2 * ryz * sy + rz * syy) + 2 * rz * (2 * tyz * sy + 2 * ty * syz) + 2 * sz * (2 * tyz * ry + 2 * ty * ryz) + 2 * tz * (2 * ry * syz + 2 * ryz * sy) + 4 * rzz * ty * sy + 4 * szz * ty * ry + 4 * tzz * sy * ry + 2 * rz * (2 * txz * sx + tz * sxx + 2 * tx * sxz + sz * txx) + 2 * sz * (2 * txz * rx + tz * rxx + 2 * tx * rxz + rz * txx) + 2 * tz * (2 * rx * sxz + sz * rxx + 2 * rxz * sx + rz * sxx) + 2 * rz * (2 * txz * sx + 2 * tx * sxz) + 2 * sz * (2 * txz * rx + 2 * tx * rxz) + 2 * tz * (2 * rx * sxz + 2 * rxz * sx) + 4 * rzz * tx * sx + 4 * szz * tx * rx + 4 * tzz * sx * rx + ty * (2 * syy * ry + 2 * sy * ryy) + ry * (3 * tyy * sy + 3 * ty * syy) + sy * (3 * tyy * ry + 3 * ty * ryy) + ty * (3 * syy * ry + 3 * sy * ryy) + ry * (2 * tyy * sy + 2 * ty * syy) + sy * (2 * tyy * ry + 2 * ty * ryy) + 2 * ry * (2 * txy * sx + ty * sxx + 2 * tx * sxy + sy * txx) + 2 * sy * (2 * txy * rx + ty * rxx + 2 * tx * rxy + ry * txx) + 2 * ty * (2 * sxy * rx + sy * rxx + 2 * sx * rxy + ry * sxx) + 2 * ry * (2 * txy * sx + 2 * tx * sxy) + 2 * sy * (2 * txy * rx + 2 * tx * rxy) + 2 * ty * (2 * sxy * rx + 2 * sx * rxy) + 2 * ryy * ty * sy + 2 * syy * ty * ry + 2 * tyy * sy * ry + 4 * ryy * tx * sx + 4 * syy * tx * rx + 4 * tyy * sx * rx + rx * (3 * txx * sx + 3 * tx * sxx) + sx * (3 * txx * rx + 3 * tx * rxx) + tx * (3 * sxx * rx + 3 * sx * rxx) + rx * (2 * txx * sx + 2 * tx * sxx) + sx * (2 * txx * rx + 2 * tx * rxx) + tx * (2 * sxx * rx + 2 * sx * rxx) + 2 * rxx * tx * sx + 2 * sxx * tx * rx + 2 * txx * sx * rx;
                                    Real urr = 4. * rx * rxxx + 4 * rxyy * rx + 4 * rxzz * rx + 3 * rxx * rxx + 2 * ryy * rxx + 2 * rzz * rxx + 4 * ry * rxxy + 4 * rz * rxxz + 4 * rxy * rxy + 4 * rxz * rxz + 4 * ry * ryyy + 4 * ryzz * ry + 3 * ryy * ryy + 2 * rzz * ryy + 4 * rz * ryyz + 4 * ryz * ryz + 4 * rz * rzzz + 3 * rzz * rzz;
                                    Real urs = 4 * rz * szzz + 4 * sz * rzzz + 6 * rzz * szz + 4 * rz * syyz + 4 * sz * ryyz + 2 * rzz * syy + 2 * szz * ryy + 4 * ryzz * sy + 8 * ryz * syz + 4 * ry * syzz + 4 * rz * sxxz + 4 * sz * rxxz + 2 * rzz * sxx + 2 * szz * rxx + 4 * rxzz * sx + 8 * rxz * sxz + 4 * rx * sxzz + 4 * sy * ryyy + 6 * syy * ryy + 4 * syyy * ry + 4 * ry * sxxy + 4 * sy * rxxy + 2 * ryy * sxx + 2 * syy * rxx + 4 * sxyy * rx + 8 * sxy * rxy + 4 * sx * rxyy + 4 * sx * rxxx + 6 * sxx * rxx + 4 * sxxx * rx;
                                    Real uss = 4 * sx * sxxx + 4 * sx * sxyy + 4 * sxzz * sx + 3 * sxx * sxx + 2 * sxx * syy + 2 * szz * sxx + 4 * sxxy * sy + 4 * sz * sxxz + 4 * sxy * sxy + 4 * sxz * sxz + 4 * sy * syyy + 4 * syzz * sy + 3 * syy * syy + 2 * szz * syy + 4 * sz * syyz + 4 * syz * syz + 4 * sz * szzz + 3 * szz * szz;
                                    Real urt = 4 * tz * rzzz + 6 * tzz * rzz + 4 * tzzz * rz + 4 * rz * tyyz + 4 * tz * ryyz + 2 * rzz * tyy + 2 * tzz * ryy + 4 * tyzz * ry + 8 * tyz * ryz + 4 * ty * ryzz + 4 * rz * txxz + 4 * tz * rxxz + 2 * rzz * txx + 2 * tzz * rxx + 4 * txzz * rx + 8 * txz * rxz + 4 * tx * rxzz + 4 * ty * ryyy + 6 * tyy * ryy + 4 * tyyy * ry + 2 * tyy * rxx + 4 * txyy * rx + 8 * txy * rxy + 4 * tx * rxyy + 4 * ry * txxy + 4 * ty * rxxy + 2 * ryy * txx + 4 * tx * rxxx + 6 * txx * rxx + 4 * txxx * rx;
                                    Real utt = 4 * tx * txxx + 4 * tx * txyy + 4 * tx * txzz + 3 * txx * txx + 2 * txx * tyy + 2 * txx * tzz + 4 * txxy * ty + 4 * txxz * tz + 4 * txy * txy + 4 * txz * txz + 4 * ty * tyyy + 4 * ty * tyzz + 3 * tyy * tyy + 2 * tyy * tzz + 4 * tyyz * tz + 4 * tyz * tyz + 4 * tz * tzzz + 3 * tzz * tzz;
                                    Real ust = 4 * tz * szzz + 6 * tzz * szz + 4 * tzzz * sz + 4 * sz * tyyz + 4 * tz * syyz + 2 * szz * tyy + 2 * tzz * syy + 4 * tyzz * sy + 8 * tyz * syz + 4 * ty * syzz + 4 * sz * txxz + 4 * tz * sxxz + 2 * szz * txx + 2 * tzz * sxx + 4 * txzz * sx + 8 * txz * sxz + 4 * tx * sxzz + 4 * ty * syyy + 6 * tyy * syy + 4 * tyyy * sy + 4 * sy * txxy + 4 * ty * sxxy + 2 * syy * txx + 2 * tyy * sxx + 4 * txyy * sx + 8 * txy * sxy + 4 * tx * sxyy + 4 * tx * sxxx + 6 * txx * sxx + 4 * txxx * sx;
                                    Real ur = rxxxx + 2 * rxxyy + 2 * rxxzz + ryyyy + 2 * ryyzz + rzzzz;
                                    Real us = sxxxx + 2 * sxxyy + 2 * sxxzz + syyyy + 2 * syyzz + szzzz;
                                    Real ut = txxxx + 2 * txxyy + 2 * txxzz + tyyyy + 2 * tyyzz + tzzzz;
                                    if( true && i1==2 && i2==2 && i3==2  )
                                    {
                                        printF(" (i1,i2,i3)=(%3d,%3d,%3d) urrrr=%g urrrs=%g, urrss=%g ursss=%g urrrt=%g urrtt=%g urttt=%g ussss=%g ussst=%g usstt=%g usttt=%g utttt=%g\n",
                                                i1,i2,i3,
                                                urrrr, urrrs, urrss, ursss, urrrt, urrtt, urttt, ussss, ussst, usstt, usttt, utttt );
                                        printF(" urrr,urrs,urss,usss,urrt,uttt,usst,ustt,urst= %g %g %g %g %g %g %g %g %g\n",urrr,urrs,urss,usss,urrt,uttt,usst,ustt,urst);
                                        printF(" urr,urs,uss,urt,utt,ust,ur,us,ut= %g %g %g %g %g %g %g %g %g\n",urr,urs,uss,urt,utt,ust,ur,us,ut);
                                        printF(" rx,ry,rz,rxx,rxy,rxz,ryy,ryz,rzz,rxxx,rxxy,rxxz,rxyy,rxzz= %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",rx,ry,rz,rxx,rxy,rxz,ryy,ryz,rzz,rxxx,rxxy,rxxz,rxyy,rxzz);
                                        printF(" sx,sy,sz,sxx,sxy,sxz,syy,syz,szz,sxxx,sxxy,sxxz,sxyy,sxzz= %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",sx,sy,sz,sxx,sxy,sxz,syy,syz,szz,sxxx,sxxy,sxxz,sxyy,sxzz);
                                        printF(" tx,ty,tz,txx,txy,txz,tyy,tyz,tzz,txxx,txxy,txxz,txyy,txzz= %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",tx,ty,tz,txx,txy,txz,tyy,tyz,tzz,txxx,txxy,txxz,txyy,txzz);
                                        printF(" rxxxx,rxxyy,rxxzz,ryyyy,ryyzz,rzzzz= %g %g %g %g %g %g\n",rxxxx,rxxyy,rxxzz,ryyyy,ryyzz,rzzzz);
                                        printF(" sxxxx,sxxyy,sxxzz,syyyy,syyzz,szzzz= %g %g %g %g %g %g\n",sxxxx,sxxyy,sxxzz,syyyy,syyzz,szzzz);
                                        printF(" txxxx,txxyy,txxzz,tyyyy,tyyzz,tzzzz= %g %g %g %g %g %g\n",txxxx,txxyy,txxzz,tyyyy,tyyzz,tzzzz);
                                        printF(" rx*sx+ry*sy+rz*sz=%g, rx*tx+ry*ty+rz*tz=%g, sx*tx+sy*ty+sz*tz=%g\n",rx*sx+ry*sy+rz*sz,rx*tx+ry*ty+rz*tz,sx*tx+sy*ty+sz*tz);
                                    }
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
                                                                                        + urrst*  rrCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
                                                                                        + ursst*   rCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
                                                                                        + urstt*   rCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 
                                                                                        + urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrs *  rrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
                                                                                        + urss *   rCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
                                                                                        + usss *   iCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrt *  rrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
                                                                                        + urtt *   rCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
                                                                                        + uttt *   iCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
                                                                                        + usst *   iCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
                                                                                        + ustt *   iCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 
                                                                                        + urst *   rCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
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
                    // save for testing : 
                                        lapSqCoeffLocal(m,i1,i2,i3) = 
                                                                                      (  urrrr*rrrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
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
                                                                                        + urrst*  rrCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
                                                                                        + ursst*   rCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
                                                                                        + urstt*   rCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 
                                                                                        + urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrs *  rrCoeff(m1)*   sCoeff(m2)*   iCoeff(m3) 
                                                                                        + urss *   rCoeff(m1)*  ssCoeff(m2)*   iCoeff(m3) 
                                                                                        + usss *   iCoeff(m1)* sssCoeff(m2)*   iCoeff(m3) 
                                                                                        + urrt *  rrCoeff(m1)*   iCoeff(m2)*   tCoeff(m3) 
                                                                                        + urtt *   rCoeff(m1)*   iCoeff(m2)*  ttCoeff(m3) 
                                                                                        + uttt *   iCoeff(m1)*   iCoeff(m2)* tttCoeff(m3) 
                                                                                        + usst *   iCoeff(m1)*  ssCoeff(m2)*   tCoeff(m3) 
                                                                                        + ustt *   iCoeff(m1)*   sCoeff(m2)*  ttCoeff(m3) 
                    // TROUBLE HERE: 
                                                                                        + urst *   rCoeff(m1)*   sCoeff(m2)*   tCoeff(m3) 
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
         // *****TEST
                    // lapSqCoeff(m,i1,i2,i3) = 
                    //                        (  urrrr*rrrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                    //                         + urrr * rrrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                    //                         + urr  *  rrCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                    //                         + ur   *   rCoeff(m1)*   iCoeff(m2)*   iCoeff(m3) 
                    //                       );
                                    } // end for stencil
                                    if( true && i1==2 && i2==2 && i3==2  )
                                    {
                    // check stencil for Delta^2 
                    // Real kx=1., ky=1., kz=1.;
                    // #define UE(j1,j2,j3) sin(kx*xLocal(j1,j2,j3,0))*sin(ky*xLocal(j1,j2,j3,1))*sin(kz*xLocal(j1,j2,j3,2))
                    // #define LAPSQE(j1,j2,j3) SQR(kx*kx+ky*ky+kz*kz) * UE(j1,j2,j3)
                                        Real lapSq = 0.;
                                        ForStencil(m1,m2,m3)
                                        {
                                            int m  = M123(m1,m2,m3); 
                      // lapSq += lapSqCoeff(m,i1,i2,i3)*UE(i1+m1,i2+m2,i3+m3);
                                            lapSq += lapSqCoeffLocal(m,i1,i2,i3)*uLocal(i1+m1,i2+m2,i3+m3);
                                        }
                    // Real lapSqe = LAPSQE(i1,i2,i3); 
                                        Real lapSqe = lapSqu(i1,i2,i3);
                                        Real err = fabs( lapSq - lapSqe )/(fabs(lapSqe) + 1.e-10 );
                                        printf("(i1,i2,i3)=(%3d,%3d,%3d) : lapSq = %10.2e, lapSqe = %10.2e, relErr= %9.2e\n",i1,i2,i3,lapSq,lapSqe, err);
                                        maxErr = max(maxErr,err);
                                        if( i1==2 )
                                        {
                                            Real maxErrCoeff=0.;
                                            printf("lapSqCoeff(m1,m2,m3) and relative errors (i1,i2,i3)=(%3d,%3d,%3d):\n",i1,i2,i3);
                                            for( int m3=-halfWidth3; m3<=halfWidth3; m3++ )
                                            {
                                                printf("m3=%2d\n",m3);
                                                for( int m2=-halfWidth2; m2<=halfWidth2; m2++ ) 
                                                {
                                                    printf("m2=%2d, m1=[%2d,%2d] : ",m2,-halfWidth1,halfWidth1);
                                                    for( int m1=-halfWidth1; m1<=halfWidth1; m1++ )
                                                    {
                                                        int m  = M123(m1,m2,m3); 
                                                        printF(" %10.2e ",lapSqCoeffLocal(m,i1,i2,i3));
                                                    }
                                                    printF("   ");
                                                    for( int m1=-halfWidth1; m1<=halfWidth1; m1++ )
                                                    {
                                                        int m  = M123(m1,m2,m3); 
                                                        Real err = fabs(lapSqCoeffLocal(m,i1,i2,i3)-coeffLapSqLocal(m,i1,i2,i3))/(fabs(coeffLapSqLocal(m,i1,i2,i3))+1.);
                                                        maxErrCoeff=max(maxErrCoeff,err);
                            // Real err = coeffLapSqLocal(m,i1,i2,i3);
                                                        printF(" %10.2e ",err);
                                                    }                  
                                                    printF("\n");
                                                }
                                            }
                                            printF(">>>(i1,i2,i3)=(%3d,%3d,%3d) : maxRelErr in coefficients = %9.2e\n",i1,i2,i3,maxErrCoeff);              
                                        }
                                    }
                                } // end if maskLocal > 0 
                            } // end for3d
              // OV_ABORT("Stop here for now");
                -- */
                        }
                    } // end curvilinear

                printF("\n DONE: Compute LapSq(u) : Max-rel-err=%9.2e\n",maxErr);

            }


        } // end for grid

   // --- compute errors in LapSq_h( u ) (FROM CHAIN RULE) -----
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            MappedGridOperators & mgop = op[grid];      

            OV_GET_SERIAL_ARRAY(Real,lapSqCoeff[grid],lapSqCoeffLocal); // coeff from chain rule

            OV_GET_SERIAL_ARRAY(Real,u[grid],uLocal);
      // OV_GET_SERIAL_ARRAY(Real,coeffLapSq[grid],coeffLapSqLocal);
            OV_GET_SERIAL_ARRAY(Real,lapSqTrue[grid],lapSqTrueLocal);
            OV_GET_SERIAL_ARRAY(Real,lapSqErr[grid],lapSqErrLocal);

      // ----- compute error in lapSq(u) -------
            realMappedGridFunction & coeff = impCoeff[grid];
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

            getIndex(mg.gridIndexRange(),I1,I2,I3);
      // int m1,m2,m3;
            
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                Real lapSq = 0.;
                ForStencil(m1,m2,m3)
                {
                    int m  = M123(m1,m2,m3); 
                    lapSq += lapSqCoeffLocal(m,i1,i2,i3)*uLocal(i1+m1,i2+m2,i3+m3);  // delta coeff
                }

                Real err = fabs( lapSq - lapSqTrueLocal(i1,i2,i3) )/lapSqTrueNorm;
                lapSqErrLocal(i1,i2,i3)=err;
            }
        }
        const Real maxErrLapSqChainRule = maxNorm(lapSqErr);

        printF(">>> res=%d: max relative error in LapSq_{2h}(uExact) = %9.2e (Coeff from Hierarchical).\n",res,maxErrLapSqChainRule );



        if( plotOption )
        {
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            psp.set(GI_TOP_LABEL,"Solution");
            ps.erase();
            PlotIt::contour( ps,lapSqErr,psp );
        }

    } // end resolutions




    Overture::finish();          
    return 0;
}
