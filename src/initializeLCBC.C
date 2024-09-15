#include "CgWave.h"
// #include "CompositeGridOperators.h";
// #include "PlotStuff.h"
#include "display.h"
// #include "ParallelOverlappingGridInterpolator.h"
// #include "ParallelUtility.h"
// #include "LoadBalancer.h"
// #include "gridFunctionNorms.h"
#include "OGPolyFunction.h"
#include "OGTrigFunction.h"
// #include "DialogData.h"
// #include "Ogshow.h"
#include "LCBC.h"
#include "ParallelUtility.h"

#define COEF_COMPUTE_HIGHER_DIM (0)
#define DEBUG_SET_EXACT (0)

#define ForBoundary(side,axis) for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) \
for( int side=0; side<=1; side++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )
#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )


/* Create Global Variable for LCBC function pointers */

#define dimBasedValue(dim,n1,n2) ((dim == 2) ? (n1) : (n2))
#undef ind2
#define ind2(i,j,N1) (i + j*(N1))

OGFunction *ogFun;
typedef double (*LCBCFn)(double *);

bool NGH_DEBUG = 0;

void printArray(double *arr, int L1, int L2);


// ==============================================================================================
/// \brief Compute the coefficients of the Laplacian in parameter space
///      Delta u = c11 u.rr + c12 u.rs + c22 u.ss + c1 u.r + c2 u.s
/// \details
/// We evaluate the coefficients at points on and near the boundary.
// ==============================================================================================
void getLcbcCoef(MappedGrid & mg, LcbcData *&LcbcCoef, RealArray coefArray[],  int p, bool evaluateAllFaces )
// void getLcbcCoef(MappedGrid & mg, double **&coef, RealArray coefArray[], RealArray & rxLocal, aString caseName, int p, int numberOfDimensions)
{
    
    const int debug = 0;
    
    const bool isRectangular     = mg.isRectangular();
    const int numberOfDimensions = mg.numberOfDimensions();
    
    int maxCoefNum = ((numberOfDimensions+1)*(numberOfDimensions + 2))/2;
    int faceNum = 2*numberOfDimensions;
    Index ib[3];
    Range ibf[3];
    Range ibf_ext[3];
    int includeGhost = p;
    int extraTangentialGhostPoints = (p-1); // Note that we are creating grid functions with extra ghost points but not evaluating the data at these ghost points. Note also that starting from July 24 update, it is enough for the data grids to have orderInSpace/2 ghost points in each direction only (The LCBC class extends the grid functions to the needed size internally)
    int lengthAxis[3], ghostAxis[3];
    
    ForBoundary(side,axis)
    {
        const int face = side + 2*axis;
        
        if( mg.boundaryCondition(side,axis)>0 || evaluateAllFaces )
        {
            
            getBoundaryIndex(mg.gridIndexRange(),side,axis,ib[0],ib[1],ib[2],includeGhost);
            const real dr[3]={mg.gridSpacing(0),mg.gridSpacing(1),mg.gridSpacing(2)};
            
            // We evaluate the coefficients at points on and near the boundary
            //
            // --- add p extra points in each normal direction
            // Boundary at r=0:
            //          +-------------+
            //          |      |      |
            //          |      |      |
            //          |      |      |
            //          |      |      |
            //          |      |      |
            //          +-------------+
            //        i1=-p   i1=0   i1=p
            
            for(int ax = 0; ax<3; ax++)
            {
                if(ax<numberOfDimensions)
                {
                    
                    if(ax == axis)
                    {
                        ibf[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                        ibf_ext[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                        
                        ghostAxis[ax] = p;
                    }
                    else
                    {
                        ibf[ax] = Range(ib[ax].getBase(),ib[ax].getBound());
                        ibf_ext[ax] = Range((ib[ax].getBase() - extraTangentialGhostPoints),(ib[ax].getBound() + extraTangentialGhostPoints));
                        ghostAxis[ax] = p + extraTangentialGhostPoints;
                    }
                }
                else
                {
                    ibf_ext[ax] = Range(0,0);
                    ibf[ax] = Range(0,0);
                    ghostAxis[ax] = 0;
                }
                lengthAxis[ax] =  ibf_ext[ax].getBound() - ibf_ext[ax].getBase() + 1;
                
            }
            for(int coefNum = 0; coefNum<maxCoefNum; coefNum++)
            {
                coefArray[ind2(face,coefNum,faceNum)].redim(ibf_ext[0],ibf_ext[1],ibf_ext[2]);
                // coefArray[ind2(face,coefNum,faceNum)] = RealArray(ibf_ext[0],ibf_ext[1],ibf_ext[2]);
            }
            
            LcbcCoef[face].initialize(lengthAxis, ghostAxis, maxCoefNum, NULL, true);

            if( isRectangular )
            {
                // ---- Rectangular grid case ----
                real dx[3]={1.,1.,1.};
                mg.getDeltaX(dx);

                if( numberOfDimensions==2 )
                {
                    // wdh: Order: c11 =rx^2 + ry^2
                    //             c22 =sx^2 + sy^2
                    //             c1  = rxx + ryy
                    //             c2  = sxx + syy
                    //             c12 = rx*sx + ry*sy
                    //             c0 ?
                    // for( int m=0; m<6; m++ )
                    // {
                    //   coefArray[ind2(face,m,faceNum)].redim(ibf[0],ibf[1],ibf[2])
                    // }
                    coefArray[ind2(face,0,faceNum)](ibf[0],ibf[1],ibf[2]) = SQR( dr[0]/dx[0] );
                    coefArray[ind2(face,1,faceNum)](ibf[0],ibf[1],ibf[2]) = SQR( dr[1]/dx[1] );
                    coefArray[ind2(face,2,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0; // c1 ux rxx + ryy
                    coefArray[ind2(face,3,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,4,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,5,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;

                }
                else
                {
                    // wdh: Order: c11 =rx^2 + ry^2 + rz^2
                    //             c22 =sx^2 + sy^2 + sz^2
                    //             c33 =tx^2 + ty^2 + tz^2
                    //             c1  = rxx + ryy
                    //             c2  = sxx + syy
                    //             c3  = sxx + syy
                    //             c12 = rx*sx + ry*sy + rz*sz
                    //             c13 = ?
                    //             c23 = rx*sx + ry*sy + rz*sz
                    //             c0  = ?
                    coefArray[ind2(face,0,faceNum)](ibf[0],ibf[1],ibf[2]) = SQR( dr[0]/dx[0] );
                    coefArray[ind2(face,1,faceNum)](ibf[0],ibf[1],ibf[2]) = SQR( dr[1]/dx[1] );
                    coefArray[ind2(face,2,faceNum)](ibf[0],ibf[1],ibf[2]) = SQR( dr[2]/dx[2] );
                    coefArray[ind2(face,3,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,4,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,5,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,6,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,7,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,8,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                    coefArray[ind2(face,9,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
                }


//                printf("Rectangular case\n") ;
//                getchar();
//                Index I1=ibf[0], I2=ibf[1], I3=ibf[2];
//                FOR_3D(i1,i2,i3,I1,I2,I3)
//                {
//                    for(int kap = 0; kap<9; kap++){
//                        printf("coef %d = %f\n", kap, coefArray[ind2(face,kap,faceNum)](i1,i2,i3)) ;
//                    }
//                }


            }
            else
            {
                // --- curvilinear grid case ---
                const IntegerArray & gid = mg.gridIndexRange();
                Mapping & map = mg.mapping().getMapping();
                
                Index I1=ibf[0], I2=ibf[1], I3=ibf[2];
                printF("getLcbcCoef (side,axis)=(%d,%d) I1=[%d,%d] I2=[%d,%d] I3=[%d,%d]\n",side,axis,
                       I1.getBase(),I1.getBound(),
                       I2.getBase(),I2.getBound(),
                       I3.getBase(),I3.getBound());
                const int numPts = I1.getLength()*I2.getLength()*I3.getLength();
                RealArray r(numPts,numberOfDimensions), xr(numPts,numberOfDimensions,numberOfDimensions);
                
                // Compute derivatives of rx by differencing
                //   First-order difference : use delta= eps^(1/3) :
                //   Set truncation error =  round-off error
                //   truncation error= delta^2    (using D0)
                //   round-off error = eps/delta
                const Real delta = pow(REAL_EPSILON,1./3.);
                
#if(COEF_COMPUTE_HIGHER_DIM==0)
                const int numCases=1 + 2*numberOfDimensions; // NGH: this is the number of pertubations of rx -> We may need more for higher order...
                RealArray *rxa = new RealArray [numCases];
                // Here are the pertubations we make for each case:
                Real delta1[] = {0., delta, -delta,     0.,     0.,     0.,     0.}; // perturb r1
                Real delta2[] = {0.,    0.,     0.,  delta, -delta,     0.,     0.}; // perturb r2
                Real delta3[] = {0.,    0.,     0.,     0.,     0.,  delta, -delta}; // perturb r3
#endif
                
#if(COEF_COMPUTE_HIGHER_DIM==1)
                const int numCases=1 + 6*numberOfDimensions; // NGH: this is the number of pertubations of rx -> We may need more for higher order...
                RealArray *rxa = new RealArray [numCases];
                // Here are the pertubations we make for each case:
                Real delta1[] = {0., -3*dr[0],-2*dr[0], -dr[0], dr[0], 2*dr[0], 3*dr[0],     0.,     0.,     0.,     0.,     0.,     0., 0., 0., 0., 0.,0.,0.}; // perturb r1
                Real delta2[] = {0.,    0.,     0., 0., 0., 0., 0., -3*dr[1],-2*dr[1], -dr[1], dr[1], 2*dr[1], 3*dr[1],      0.,     0., 0., 0.,0.,0.}; // perturb r2
                Real delta3[] = {0.,    0.,     0.,     0.,     0., 0., 0., 0., 0.,0.,0.,0.,0., -3*dr[2], -2*dr[2], -dr[2], dr[2], 2*dr[2],3*dr[2]}; // perturb r3
#endif
                
                for( int k=0; k<numCases; k++ )
                {
                    int i=0;
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        // parameter space coordinates
                        r(i,0) = (i1-gid(0,0))*dr[0] + delta1[k];
                        r(i,1) = (i2-gid(0,1))*dr[1] + delta2[k];
                        if( numberOfDimensions==3)
                            r(i,2) = (i3-gid(0,2))*dr[2] + delta3[k];
                        i++;
                    }
                    
                    // evaluate the Jacobian matrix
                    map.mapS( r,Overture::nullRealArray(),xr );
                    
                    if( debug & 4 )
                    {
                        ::display(r,"r");
                        ::display(xr,"xr");
                    }
                    
                    // --- Invert the Jacobian matrix to get rx array ----
                    
                    // fill in this array:
                    RealArray & rx = rxa[k];
                    rx.redim(numPts,numberOfDimensions,numberOfDimensions);
                    
                    Real det;
                    for( int i=0; i<numPts; i++ )
                    {
                        if( numberOfDimensions==2 )
                        {
                            det = xr(i,0,0) * xr(i,1,1) - xr(i,0,1) * xr(i,1,0);
                            det=1./det;
                            rx(i,0,0) = xr(i,1,1) * det;
                            rx(i,1,0) =-xr(i,1,0) * det;
                            rx(i,0,1) =-xr(i,0,1) * det;
                            rx(i,1,1) = xr(i,0,0) * det;
                        }
                        else
                        {
                            det =( (xr(i,0,0)*xr(i,1,1)-xr(i,0,1)*xr(i,1,0))*xr(i,2,2) +
                                  (xr(i,0,1)*xr(i,1,2)-xr(i,0,2)*xr(i,1,1))*xr(i,2,0) +
                                  (xr(i,0,2)*xr(i,1,0)-xr(i,0,0)*xr(i,1,2))*xr(i,2,1) );
                            det=1./det; // NGH: this was missing before. This should fix the bug in the 3D curvilinear case.
                            rx(i,0,0)=(xr(i,1,1)*xr(i,2,2)-xr(i,1,2)*xr(i,2,1))*det;
                            rx(i,1,0)=(xr(i,1,2)*xr(i,2,0)-xr(i,1,0)*xr(i,2,2))*det;
                            rx(i,2,0)=(xr(i,1,0)*xr(i,2,1)-xr(i,1,1)*xr(i,2,0))*det;
                            rx(i,0,1)=(xr(i,2,1)*xr(i,0,2)-xr(i,2,2)*xr(i,0,1))*det;
                            rx(i,1,1)=(xr(i,2,2)*xr(i,0,0)-xr(i,2,0)*xr(i,0,2))*det;
                            rx(i,2,1)=(xr(i,2,0)*xr(i,0,1)-xr(i,2,1)*xr(i,0,0))*det;
                            rx(i,0,2)=(xr(i,0,1)*xr(i,1,2)-xr(i,0,2)*xr(i,1,1))*det;
                            rx(i,1,2)=(xr(i,0,2)*xr(i,1,0)-xr(i,0,0)*xr(i,1,2))*det;
                            rx(i,2,2)=(xr(i,0,0)*xr(i,1,1)-xr(i,0,1)*xr(i,1,0))*det;
                        }
                    }
                    
                } // end for k
                
                RealArray & rx = rxa[0]; // rx (un-perturbed)
                if( numberOfDimensions==2 )
                { // ---- two dimensions -----
                    RealArray rxx(2,2), rxy(2,2);
                    int i=0;
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        // -- compute rx.x, rx.y, etc. ---
                        // we actually compute more values than we use below
                        for( int m2=0; m2<numberOfDimensions; m2++ )
                        {
                            for( int m1=0; m1<numberOfDimensions; m1++ )
                            {
                                
#if(COEF_COMPUTE_HIGHER_DIM==1)
                                Real rxr = (-rxa[1](i,m1,m2)+9.*rxa[2](i,m1,m2)-45.*rxa[3](i,m1,m2)+45.*rxa[4](i,m1,m2)-9.*rxa[5](i,m1,m2)+rxa[6](i,m1,m2))/(60.*dr[0]);  // finite difference approx.
                                Real rxs = (-rxa[7](i,m1,m2)+9.*rxa[8](i,m1,m2)-45.*rxa[9](i,m1,m2)+45.*rxa[10](i,m1,m2)-9.*rxa[11](i,m1,m2)+rxa[12](i,m1,m2))/(60.*dr[1]);
#endif
#if(COEF_COMPUTE_HIGHER_DIM==0)
                                Real rxr = (rxa[1](i,m1,m2)-rxa[2](i,m1,m2))/(2.*delta);  // finite difference approx.
                                Real rxs = (rxa[3](i,m1,m2)-rxa[4](i,m1,m2))/(2.*delta);
#endif
                                rxx(m1,m2) = rx(i,0,0)*rxr + rx(i,1,0)*rxs;
                                rxy(m1,m2) = rx(i,0,1)*rxr + rx(i,1,1)*rxs;
                                
                            }
                        }
                        
                        coefArray[ind2(face,0,faceNum)](i1,i2,i3) = SQR(rx(i,0,0)) + SQR(rx(i,0,1));
                        coefArray[ind2(face,1,faceNum)](i1,i2,i3) = SQR(rx(i,1,0)) + SQR(rx(i,1,1));
                        coefArray[ind2(face,2,faceNum)](i1,i2,i3) = rxx(0,0) + rxy(0,1);                        // c1 = rxx + ryy
                        coefArray[ind2(face,3,faceNum)](i1,i2,i3) = rxx(1,0) + rxy(1,1);                        // c2 = sxx + syy
                        coefArray[ind2(face,4,faceNum)](i1,i2,i3) = 2.0*(rx(i,0,0)*rx(i,1,0) + rx(i,0,1)*rx(i,1,1));  // c12 = 2.0*(rx*sx + ry*sy) NGH: added factor of 2
                        coefArray[ind2(face,5,faceNum)](i1,i2,i3) = 0.0;                                        // c0
                        i++;
                    }
                }
                else
                {
                    // ---- three dimensions -----
                    RealArray rxx(3,3), rxy(3,3), rxz(3,3);
                    int i=0;
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        // -- compute rx.x, rx.y, etc. ---
                        // we actually compute more values than we use below
                        for( int m2=0; m2<numberOfDimensions; m2++ )
                        {
                            for( int m1=0; m1<numberOfDimensions; m1++ )
                            {
#if(COEF_COMPUTE_HIGHER_DIM==1)
                                Real rxr = (-rxa[1](i,m1,m2)+9.*rxa[2](i,m1,m2)-45.*rxa[3](i,m1,m2)+45.*rxa[4](i,m1,m2)-9.*rxa[5](i,m1,m2)+rxa[6](i,m1,m2))/(60.*dr[0]);
                                Real rxs = (-rxa[7](i,m1,m2)+9.*rxa[8](i,m1,m2)-45.*rxa[9](i,m1,m2)+45.*rxa[10](i,m1,m2)-9.*rxa[11](i,m1,m2)+rxa[12](i,m1,m2))/(60.*dr[1]);
                                Real rxt = (-rxa[13](i,m1,m2)+9.*rxa[14](i,m1,m2)-45.*rxa[15](i,m1,m2)+45.*rxa[16](i,m1,m2)-9.*rxa[17](i,m1,m2)+rxa[18](i,m1,m2))/(60.*dr[2]);
#endif
#if(COEF_COMPUTE_HIGHER_DIM==0)
                                Real rxr = (rxa[1](i,m1,m2)-rxa[2](i,m1,m2))/(2.*delta);  // finite difference approx.
                                Real rxs = (rxa[3](i,m1,m2)-rxa[4](i,m1,m2))/(2.*delta);
                                Real rxt = (rxa[5](i,m1,m2)-rxa[6](i,m1,m2))/(2.*delta);
#endif
                                
                                rxx(m1,m2) = rx(i,0,0)*rxr + rx(i,1,0)*rxs + rx(i,2,0)*rxt;
                                rxy(m1,m2) = rx(i,0,1)*rxr + rx(i,1,1)*rxs + rx(i,2,1)*rxt;
                                rxz(m1,m2) = rx(i,0,2)*rxr + rx(i,1,2)*rxs + rx(i,2,2)*rxt;
                            }
                        }
                        // **CHECK ME**
                        coefArray[ind2(face,0,faceNum)](i1,i2,i3) = SQR(rx(i,0,0)) + SQR(rx(i,0,1)) + SQR(rx(i,0,2)); // c11 = (r_x)^2 + (r_y)^2 _ (r_z)^2
                        coefArray[ind2(face,1,faceNum)](i1,i2,i3) = SQR(rx(i,1,0)) + SQR(rx(i,1,1)) + SQR(rx(i,1,2)); // c22 = (s_x)^2 + (s_y)^2 _ (s_z)^2
                        coefArray[ind2(face,2,faceNum)](i1,i2,i3) = SQR(rx(i,2,0)) + SQR(rx(i,2,1)) + SQR(rx(i,2,2)); // c33 = (t_x)^2 + (t_y)^2 _ (t_z)^2
                        
                        /* Notation: rxx = [rx]_x = [rxx rxy; sxx sxy] */
                        coefArray[ind2(face,3,faceNum)](i1,i2,i3) = rxx(0,0) + rxy(0,1) + rxz(0,2); // c1 = rxx + ryy + rzz
                        coefArray[ind2(face,4,faceNum)](i1,i2,i3) = rxx(1,0) + rxy(1,1) + rxz(1,2); // c2 = sxx + syy + szz
                        coefArray[ind2(face,5,faceNum)](i1,i2,i3) = rxx(2,0) + rxy(2,1) + rxz(2,2); // c3 = txx + tyy + tzz
                        
                        /* Nour: coefficients fixed */
                        coefArray[ind2(face,6,faceNum)](i1,i2,i3) = 2.0*(rx(i,0,0)*rx(i,1,0) + rx(i,0,1)*rx(i,1,1) + rx(i,0,2)*rx(i,1,2)); // c12 : 2*(rx*sx + ry*sy + rz*sz) NGH: added factor of 2
                        coefArray[ind2(face,7,faceNum)](i1,i2,i3) = 2.0*(rx(i,0,0)*rx(i,2,0) + rx(i,0,1)*rx(i,2,1) + rx(i,0,2)*rx(i,2,2)); // c13 : 2*(rx*tx + ry*ty + rz*tz) NGH: added factor of 2
                        coefArray[ind2(face,8,faceNum)](i1,i2,i3) = 2.0*(rx(i,1,0)*rx(i,2,0) + rx(i,1,1)*rx(i,2,1) + rx(i,1,2)*rx(i,2,2)); // c23 : 2*(sx*tx + sy*ty + sz*tz)  NGH: added factor of 2
                        
                        coefArray[ind2(face,9,faceNum)](i1,i2,i3) = 0.0;  // c0
                        
                        
//                        printf("Curvilinear case\n") ;
//                        getchar();
//                        for(int kap = 0; kap<9; kap++)
//                            printf("coef %d = %f\n", kap, coefArray[ind2(face,kap,faceNum)](i1,i2,i3));
                        
                        i++;
                    }
                }
                
                delete [] rxa;
            }// end else if not Rectangular
            
            for(int coefNum = 0; coefNum<maxCoefNum; coefNum++)
            {
                LcbcCoef[face].Fn[coefNum] = coefArray[ind2(face,coefNum,faceNum)].getDataPointer();
            }
            
        } // end if mg.boundaryCondition(side,axis)>0
    }// end of for Boundary
}


// ========================================================================================================
/// \brief Assign data for LCBC routines
// ========================================================================================================
void CgWave::getLcbcData( MappedGrid & mg, LcbcData *&fn, LcbcData *&gn, RealArray tmpGn[], RealArray tmpFn[], Real t )
{
    
    // printF(">>> getLcbcData: evaluate BC data at t=%9.3e\n",t);
    
    const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const bool twilightZone = forcingOption==twilightZoneForcing;
    
    const int numberOfDimensions       = mg.numberOfDimensions();
    const int & orderOfAccuracy        = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInSpace = orderOfAccuracy;
    const int p = orderOfAccuracyInSpace/2;
    const int degreeInTime = dbase.get<int>("degreeInTime");
    const TwilightZoneEnum & tzType = dbase.get<TwilightZoneEnum>("twilightZone");
    

    const bool & evaluateAllFaces = dbase.get<bool>("evaluateAllFaces");
    
    OGFunction & e  = *dbase.get<OGFunction*>("tz");
    
    OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);             // array of grid points xLocal(i1,i2,i3,0:nd-1)
    
    int faceNum = 2*numberOfDimensions;
    Index ib[3];
    Range ibf[3];
    Range ibf_ext[3];
    int includeGhost = p;
    int extraTangentialGhostPoints = (p-1);
    Range C = 1;
    int q = 2;
    int lengthAxis[3], ghostAxis[3] = {p,p,dimBasedValue(numberOfDimensions,0,p)};
    int gn_ghostAxis[3] = {p,p,dimBasedValue(numberOfDimensions,0,p)};
    int gn_lengthAxis[3];
    
    ForBoundary(side,axis)
    {
        int face = side + 2*axis;
        int bcType = mg.boundaryCondition(side,axis);
        if( mg.boundaryCondition(side,axis)>0 || evaluateAllFaces )
        {
            // faceEval[face] = true;
            //---- prepare the boundary grid functions for LCBC ----
            
            getBoundaryIndex(mg.gridIndexRange(),side,axis,ib[0],ib[1],ib[2],includeGhost);
            
            for(int ax = 0; ax<3; ax++)
            {
                if(ax<numberOfDimensions)
                {
                    if(ax == axis)
                    {
                        ibf[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                        ibf_ext[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                        ghostAxis[ax] = p;
                    }
                    else
                    {
                        ibf[ax] = Range(ib[ax].getBase(),ib[ax].getBound());
                        ibf_ext[ax] = Range((ib[ax].getBase() - extraTangentialGhostPoints),(ib[ax].getBound() + extraTangentialGhostPoints));
                        ghostAxis[ax] = p + extraTangentialGhostPoints;
                    }
                }
                else
                {
                    ibf_ext[ax] = Range(0,0);
                    ibf[ax] = Range(0,0);
                    ghostAxis[ax] = 0;
                }
                
                lengthAxis[ax] =  ibf_ext[ax].getBound() - ibf_ext[ax].getBase() + 1;
                gn_lengthAxis[ax] =  ib[ax].getBound() - ib[ax].getBase() + 1;
            }
            
            bool fixData = true;
            gn[face].initialize(gn_lengthAxis,gn_ghostAxis,(p+1));
            fn[face].initialize(lengthAxis,ghostAxis,(p+1),NULL,fixData);
            
            if( twilightZone )
            {
                // --- assign LCBS data for twilight-zone ----
                RealArray sol(ib[0],ib[1],ib[2]);    sol = 0.;
                RealArray utt(ibf[0],ibf[1],ibf[2]); utt = 0.;
                RealArray uxx(ibf[0],ibf[1],ibf[2]); uxx = 0.;
                RealArray uyy(ibf[0],ibf[1],ibf[2]); uyy = 0.;
                RealArray uzz(ibf[0],ibf[1],ibf[2]); uzz = 0.;
                
                // Evaluate time derivatives of the solution and forcing
                for(int nu = 0; nu<(p+1); nu++)
                {
                    if(!(tzType==polynomial && (q*(nu+1))>degreeInTime)){
                        e.gd( utt,xLocal,numberOfDimensions,0,(q*(nu+1)),0,0,0,ibf[0],ibf[1],ibf[2],C,t);
                    }else{
                        utt = 0.;
                    }
                    if(!(tzType==polynomial && (q*nu)>degreeInTime)){

                        int normalDeriv[3] = {0,0,0};
                        normalDeriv[axis] = ((bcType == 1)?(0):(1));
                        e.gd(sol,xLocal,numberOfDimensions,0,(q*nu),normalDeriv[0],normalDeriv[1],normalDeriv[2],ib[0],ib[1],ib[2],C,t);

                        e.gd( uxx,xLocal,numberOfDimensions,0,(q*nu),2,0,0,ibf[0],ibf[1],ibf[2],C,t);
                        e.gd( uyy,xLocal,numberOfDimensions,0,(q*nu),0,2,0,ibf[0],ibf[1],ibf[2],C,t);
                        if(numberOfDimensions == 3)
                        {
                            e.gd( uzz,xLocal,numberOfDimensions,0,(q*nu),0,0,2,ibf[0],ibf[1],ibf[2],C,t);
                        }
                    }else{
                        sol = 0.;
                        uxx = 0.;
                        uyy = 0.;
                        uzz = 0.;
                    }

                    tmpGn[ind2(face,nu,faceNum)].redim(ib[0],ib[1],ib[2]);
                    tmpGn[ind2(face,nu,faceNum)] = sol(ib[0],ib[1],ib[2]);
                    gn[face].Fn[nu] = tmpGn[ind2(face,nu,faceNum)].getDataPointer();

                    tmpFn[ind2(face,nu,faceNum)].redim(ibf_ext[0],ibf_ext[1],ibf_ext[2]); // created over extended range
                    tmpFn[ind2(face,nu,faceNum)](ibf[0],ibf[1],ibf[2]) = utt(ibf[0],ibf[1],ibf[2]) - uxx(ibf[0],ibf[1],ibf[2]) - uyy(ibf[0],ibf[1],ibf[2]) - uzz(ibf[0],ibf[1],ibf[2]);
                    fn[face].Fn[nu] = tmpFn[ind2(face,nu,faceNum)].getDataPointer();
                    
                } // end for nu
            }
            else
            {
                // --- fill in zero data ---
                
                // RealArray sol(ib[0],ib[1],ib[2]);    sol = 0.;
                // RealArray utt(ibf[0],ibf[1],ibf[2]);  utt  = 0.;
                for(int nu = 0; nu<(p+1); nu++)
                {
                    tmpGn[ind2(face,nu,faceNum)].redim(ib[0],ib[1],ib[2]);
                    tmpGn[ind2(face,nu,faceNum)] = 0.;
                    gn[face].Fn[nu] = tmpGn[ind2(face,nu,faceNum)].getDataPointer();
                    
                    
                    tmpFn[ind2(face,nu,faceNum)].redim(ibf_ext[0],ibf_ext[1],ibf_ext[2]); // created over extended range
                    tmpFn[ind2(face,nu,faceNum)]=0;
                    fn[face].Fn[nu] = tmpFn[ind2(face,nu,faceNum)].getDataPointer();
                }
                
            }
        }// end of if statement
    }// end of ForBoundary
}


// ===============================================================================
/// \brief Initialize the LCBC objects -- local compatibility boundary condtions
// ===============================================================================
int CgWave::initializeLCBC()
{
    Real cpu0 = getCPU();
    
    const int & debug                                 = dbase.get<int>("debug");
    const int & orderOfAccuracy                       = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInSpace                = orderOfAccuracy;
    const int & orderOfAccuracyInTime                 = dbase.get<int>("orderOfAccuracyInTime");
    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");
    const ForcingOptionEnum & forcingOption           = dbase.get<ForcingOptionEnum>("forcingOption");
    const bool twilightZone = forcingOption==twilightZoneForcing;
    
    const int numberOfDimensions = cg.numberOfDimensions();
    
    if( !dbase.has_key("evaluateAllFaces") )
        dbase.put<bool>("evaluateAllFaces");
    
    bool & evaluateAllFaces = dbase.get<bool>("evaluateAllFaces");
    evaluateAllFaces=false; // if true, evaluate all faces independent of the bc
    
    
    int p = orderOfAccuracyInSpace/2;
    
    
    if( !dbase.has_key("LCBC") )
    {
        Lcbc *& pLcbc = dbase.put<Lcbc*>("LCBC");
        pLcbc = new Lcbc [cg.numberOfComponentGrids()];
    }
    
    Lcbc *& lcbc = dbase.get<Lcbc*>("LCBC");
    
    // --- Create the twilightzone exact solution ----
    // OGFunction *tz=NULL;
    // createTwilightZoneSolution( twilightZone, tz, numberOfDimensions, degreeInSpace, degreeInTime, trigFreq );
    // assert( tz!=NULL );
    // OGFunction & e = *tz;
    
    OGFunction *tz  = dbase.get<OGFunction*>("tz");
    OGFunction & e = *tz;
    
    /* NEW: Adding this line to create function pointers for LCBC */
    ogFun = tz;
    LCBCFn coefFnPtr=NULL, forcingFnPtr=NULL, bdryFnPtr=NULL; // specify PDE coeff etc with functions
    
    const int faceNum = 2*numberOfDimensions;
    if( !dbase.has_key("faceEval") )
    {
        IntegerArray & faceEval = dbase.put<IntegerArray>("faceEval");  // Nour expects this array to remain around

            faceEval.redim(6,cg.numberOfComponentGrids());
        
        
        faceEval=0;
    }
    
    
    Index I1,I2,I3;    // index values for the interior
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        // build geometry arrays
        mg.update(MappedGrid::THEmask | MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEinverseVertexDerivative );
        
        OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);                   // array of grid points xLocal(i1,i2,i3,0:nd-1)
        // OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal); // metrics: rxLocal(i1,i2,i3,0:*) rx, ry, sx, sy
        
        getIndex(cg[grid].dimension(),I1,I2,I3);
        
        printF("initializeLCBC: grid=%d orderOfAccuracyInSpace=%d, orderOfAccuracyInTime=%d\n",
               grid,orderOfAccuracyInSpace,orderOfAccuracyInTime);
        
        // ::display(mg.dimension(),"dimension","%5i ");             // index bounds for all grid points including ghost
        // ::display(mg.gridIndexRange(),"gridIndexRange","%5i ");   // index bounds for grid points including boundaries
        
        
        /* ===== Grid information needed for LCBC ===== */
        
        printf("\n======= p = %d ======\n",p);
        
        int userNumGhost = -mg.dimension(0,0); // **FIX ME**
        int Nx[3]={mg.gridIndexRange(1,0),mg.gridIndexRange(1,1),mg.gridIndexRange(1,2)};
        if( debug & 4 )
            printf("initializeLCBC:userNumGhost = %d and Nx = %d,%d,%d\n", userNumGhost,  Nx[0], Nx[1], Nx[2]);
        
        IntegerArray & faceEval = dbase.get<IntegerArray>("faceEval");
        ForBoundary(side,axis)
        {
            int face = side + 2*axis;
            
//            faceEval(face,grid) = mg.boundaryCondition(side,axis);
            
            if( mg.boundaryCondition(side,axis)>0 || evaluateAllFaces )
            {
                faceEval(face,grid) = mg.boundaryCondition(side,axis);
            }
            else if(mg.boundaryCondition(side,axis)<0)
            {
                faceEval(face,grid) = -1;
            }
            else
            {
                faceEval(face,grid) = 0;
            }
        }// end of ForBoundary
        
        /* ===== BUILD the coefficient grid ===== */
        
        int maxCoefNum = ((numberOfDimensions+1)*(numberOfDimensions + 2))/2;
        
        // RealArray coefArray[(maxCoefNum*faceNum)];
        RealArray *coefArray = new RealArray [(maxCoefNum*faceNum)];                   // *** WDH delete me
        LcbcData *LcbcCoef = new LcbcData[faceNum];
        //************************************************************
        
        getLcbcCoef( mg, LcbcCoef, coefArray, p, evaluateAllFaces );
        
        /* ===================================== */
        
        
        /* Call the LCBC initialize */
        printF("initializeLCBC: Call the LCBC constructor...\n");
        
        
        bool LCBC_cstCoef, LCBC_noForcing, *LCBC_zeroBC;
//        if(twilightZone){
           LCBC_zeroBC = NULL;
           LCBC_noForcing = false;
           LCBC_cstCoef = false;
//        }
//        else
//        {
//            LCBC_cstCoef = false;
//            LCBC_noForcing = true;
//
//            LCBC_zeroBC = new bool[(2*numberOfDimensions)];
//            for(int face = 0; face<(2*numberOfDimensions); face++)
//                LCBC_zeroBC[face] = true;
//        }
        
        lcbc[grid].initialize(numberOfDimensions,orderOfAccuracyInSpace,orderOfAccuracyInTime,Nx,userNumGhost,
                              &faceEval(0,grid), LcbcCoef, coefFnPtr, bdryFnPtr, forcingFnPtr, NULL, LCBC_cstCoef, LCBC_zeroBC, LCBC_noForcing);
        
        delete [] LCBC_zeroBC;
        
        // OLD:
        // lcbc[grid].initialize(numberOfDimensions,orderOfAccuracyInSpace,orderOfAccuracyInTime,Nx,userNumGhost, coef, coefFnPtr, bdryFnPtr, forcingFnPtr);
        printF("... done LCBC constructor\n");
        
        //   delete [] coef;
        //   coef = NULL;
        
        printf("initializeLCBC: Finished grid = %d\n", grid);
        
    } // end for grid
    
    Real cpu = getCPU() - cpu0;
    printF("initializeLCBC: time to initialize LCBC = %9.3e (s)\n",cpu);
    
    
    return 0;
}

// ===============================================================================
/// \brief Assign the LCBC boundary conditions
// ===============================================================================
int CgWave::assignLCBC( realMappedGridFunction & u, Real t, Real dt, int grid )
{
    Real cpu0 = getCPU();
    
    const int & debug                                 = dbase.get<int>("debug");
    const int & orderOfAccuracy                       = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInSpace                = orderOfAccuracy;
    const int & orderOfAccuracyInTime                 = dbase.get<int>("orderOfAccuracyInTime");
    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");
    const ForcingOptionEnum & forcingOption           = dbase.get<ForcingOptionEnum>("forcingOption");
    const bool twilightZone = forcingOption==twilightZoneForcing;
    
    const int numberOfDimensions = cg.numberOfDimensions();
    
    int p = orderOfAccuracyInSpace/2;
    
    
    // if( !dbase.has_key("LCBC") )
    // {
    //   Lcbc *& pLcbc = dbase.put<Lcbc*>("LCBC");
    //   pLcbc = new Lcbc [cg.numberOfComponentGrids()];
    // }
    
    Lcbc *& lcbc = dbase.get<Lcbc*>("LCBC");
    
    // --- Create the twilightzone exact solution ----
    // OGFunction *tz=NULL;
    // createTwilightZoneSolution( twilightZone, tz, numberOfDimensions, degreeInSpace, degreeInTime, trigFreq );
    // assert( tz!=NULL );
    // OGFunction & e = *tz;
    
    OGFunction *tz  = dbase.get<OGFunction*>("tz");
    OGFunction & e = *tz;
    
    /* NEW: Adding this line to create function pointers for LCBC */
    ogFun = tz;
    
    
    Index I1,I2,I3;    // index values for the interior
    
    MappedGrid & mg = *u.getMappedGrid();
    
    // build geometry arrays
    // mg.update(MappedGrid::THEmask | MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEinverseVertexDerivative );
    
    OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);                   // array of grid points xLocal(i1,i2,i3,0:nd-1)
    // OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal); // metrics: rxLocal(i1,i2,i3,0:*) rx, ry, sx, sy
    
    OV_GET_SERIAL_ARRAY(real,u,uLocal);
    
    // ***** Nour added mask here *****
    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
    int *pmask = maskLocal.getDataPointer();
    // ***** end of mask addition *****
    
    getIndex( mg.dimension(),I1,I2,I3);
    
//     printF("initializeLCBC: grid=%d orderOfAccuracyInSpace=%d, orderOfAccuracyInTime=%d\n",
//           grid,orderOfAccuracyInSpace,orderOfAccuracyInTime);
//
//     ::display(mg.dimension(),"dimension","%5i ");             // index bounds for all grid points including ghost
//     ::display(mg.gridIndexRange(),"gridIndexRange","%5i ");   // index bounds for grid points including boundaries
    
    /* ===== Grid information needed for LCBC ===== */
    
    // printf("\n======= p = %d ======\n",p);
    
    int userNumGhost = -mg.dimension(0,0);
    int Nx[3]={mg.gridIndexRange(1,0),mg.gridIndexRange(1,1),mg.gridIndexRange(1,2)};
    if( debug & 4 )
        printf("userNumGhost = %d and Nx = %d,%d,%d\n", userNumGhost,  Nx[0], Nx[1], Nx[2]);
    
    
    const int faceNum = 2*numberOfDimensions;
    
    //  =========== EVALUATE BOUNDARY DATA ================
    
    
    int NU = (p+1);
    LcbcData *gn = new LcbcData[faceNum];
    LcbcData *fn = new LcbcData[faceNum];
    
    RealArray tmpFn[(faceNum*NU)];
    RealArray tmpGn[(faceNum*NU)];
    
    // bool & evaluateAllFaces = dbase.get<bool>("evaluateAllFaces");
    getLcbcData( mg, fn, gn, tmpGn, tmpFn, t );
    
    /*========== END OF EVALUATE BOUNDARYDATA  ================ */
    
    
    /* Update the solution at the ghost points */
    double *un = uLocal.getDataPointer();
    // printf("assignLCBC: updating %d ghost points normal to each boundary point\n",(orderOfAccuracyInSpace/2));
    
#if DEBUG_SET_EXACT == 1
//    if(grid == 0){
        Range all;
        realCompositeGridFunction ue(cg,all,all,all);   // holds true solution
        ue.setName("ue",0);
        ue=0.;
        
        OV_GET_SERIAL_ARRAY(real,ue[grid],ueLocal);
        
        Index Ie1,Ie2,Ie3;
        getIndex(cg[grid].dimension(),Ie1,Ie2,Ie3);// assign all points including ghost points
        
        int numberOfComponents=1;
        Range C=numberOfComponents;
        int isRectangular=0;
        int ntd=0, nxd=0, nyd=0, nzd=0;
        e.gd( ueLocal ,xLocal,numberOfDimensions,isRectangular,ntd,nxd,nyd,nzd,Ie1,Ie2,Ie3,C,t); // assign exact solution at time t
        double *ue_LCBC = ueLocal.getDataPointer();
        lcbc[grid].setExactGhost(un, ue_LCBC, pmask);
//    }
#endif

    // ============== CALL LCBC =========================
    lcbc[grid].updateGhost( un, pmask, t, dt, gn, fn );
    
    /* NOUR: for debugging */
//    getchar();
//    printf("grid = %d\n", grid);
//    printIndexRange(lcbc[grid].maskRange);
////    lcbc[grid].printBdryRange();
////    lcbc[grid].printMaskRange(pmask);
//    getchar();
    /* ------------------- */
    
    if( debug & 2 )
    {
        Real cpu = getCPU() - cpu0;
        printF("assignLCBC: t=%9.3e, cpu = %9.3e (s)\n",t,cpu);
    }
    
    return 0;
}


void printArray(double *arr, int L1, int L2){
    printf("Printing Array of size (%d,%d)\n", L1, L2);
    for(int i = 0; i<L1; i++){
        for(int j = 0; j<L2; j++){
            printf("%f|", arr[(i + L1*j)]);
        }// end of j loop
        printf("\n");
    }// end of i loop
    printf("\n");
}// end of print2Darray
    

#undef DEBUG_SET_EXACT
