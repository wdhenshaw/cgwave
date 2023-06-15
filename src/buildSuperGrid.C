#include "CgWave.h"
#include "ParallelUtility.h"
#include "PlotStuff.h"
#include "display.h"

// ***********************************************************************************************
//    Build SuperGrid Layer functions
// 
// Note: This file started from cg/mx/src/buildSuperGrid.C
// ***********************************************************************************************

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  
#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )


// ================================================================================================
///  \brief  Adjust index bounds to account for the an absorbing BC (.e.g. supergrid)
//   when computing errors or plotting 
//
/// /param extra (input) : an additional offset 
/// /Iv (input/output) :
/// /Return value: true if the Index Iv was changed.
// ================================================================================================
bool CgWave::
adjustBoundsForAbsorbingLayer( MappedGrid & mg, Index Iv[3], int extra /* =0 */ )
{
  bool useAbsorbing = (mg.boundaryCondition(0,0)==absorbing || mg.boundaryCondition(1,0)==absorbing ||
                       mg.boundaryCondition(0,1)==absorbing || mg.boundaryCondition(1,1)==absorbing ||
                       mg.boundaryCondition(0,2)==absorbing || mg.boundaryCondition(1,2)==absorbing);
  
  if( !useAbsorbing ) return false;

  const bool isRectangular=mg.isRectangular();
  assert( isRectangular );  // 
  real dx[3]={1.,1.,1.};
  mg.getDeltaX( dx );
  const IntegerArray & gid = mg.gridIndexRange();
      
  const real & superGridWidth = dbase.get<real>("superGridWidth");

  //printF(" +++ adjustBoundsForAbsorbingLayer: superGridWidth=%g, extra=%d\n",superGridWidth,extra);
  //printF(" Input: Iv=[%d,%d][%d,%d]\n",Iv[0].getBase(),Iv[0].getBound(),Iv[1].getBase(),Iv[1].getBound());
  for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
  {
    int offset =  int( superGridWidth/dx[axis] + .5 ) + extra;  // what should this be ? 
    //printF("[axis=%d dx=%e offset= %d]",axis,dx[axis],offset);

    int na=Iv[axis].getBase();
    if( mg.boundaryCondition(0,axis)==absorbing )
    {
      na = max(na,gid(0,axis) + offset);
      na = min(na, int( gid(0,axis)+(gid(0,axis)+gid(1,axis)-1)/2  ));
    }
    
    int nb=Iv[axis].getBound();
    if( mg.boundaryCondition(1,axis)==absorbing )
    {
      nb = min(nb,gid(1,axis) - offset);
      nb = max(nb, int( gid(0,axis)+ (gid(0,axis)+gid(1,axis)+1)/2  ));
    }
    

    Iv[axis]=Range(na,nb);
  }
  //printF("\n");
  //printF(" Output: Iv=[%d,%d][%d,%d]\n",Iv[0].getBase(),Iv[0].getBound(),Iv[1].getBase(),Iv[1].getBound());

  
  return useAbsorbing;
}




// =================================================================================================
/// \brief Build the super-grid absorbing layer functions 
// =================================================================================================
int CgWave::
buildSuperGrid( )
{
  int & useSuperGrid = dbase.get<int>("useSuperGrid");

  // useSuperGrid=true;  // *** TEST ****
  
  if( useSuperGrid )
  {
    printF("\n ++++++Build Supergrid absorbing layer functions +++++++\n");

    // Here is the CompositeGrid: 
    // assert( cgp!=NULL );
    // CompositeGrid & cg = *cgp;
    const int numberOfComponentGrids = cg.numberOfComponentGrids();
    const int numberOfDimensions = cg.numberOfDimensions();

    // --- Build the Super-grid layer functions ----

    if( !dbase.has_key("etaxSuperGrid") )
    {
      dbase.put<RealArray*>("etaxSuperGrid" )=NULL;  // ** delete me ** 
      dbase.put<RealArray*>("etaySuperGrid" )=NULL;
      dbase.put<RealArray*>("etazSuperGrid" )=NULL;
      
    }
    RealArray *& etaxSuperGrid = dbase.get<RealArray*>("etaxSuperGrid" );
    RealArray *& etaySuperGrid = dbase.get<RealArray*>("etaySuperGrid" );
    RealArray *& etazSuperGrid = dbase.get<RealArray*>("etazSuperGrid" );
    
    etaxSuperGrid = new RealArray[numberOfComponentGrids];
    etaySuperGrid = new RealArray[numberOfComponentGrids];
    etazSuperGrid = new RealArray[numberOfComponentGrids];

    real epsL,pSG,qSG;
    real ax,bx,ay,by;
    
#define sigmaSG(z) ( (1.-epsL)*pow( 1 - pow( (1 -(z)/superGridWidth),pSG),qSG) )
#define etaSG(z)  (1. -sigmaSG(z)) 

#define sigmaz(z) ( (1.-epsL)* qSG*pow( 1 - pow( (1 -(z)/superGridWidth),pSG),qSG-1) * (pSG/superGridWidth) * pow( (1 -(z)/superGridWidth),pSG-1) )
#define etaSGz(z) ( -(1.-epsL)* qSG*pow( 1 - pow( (1 -(z)/superGridWidth),pSG),qSG-1) * (pSG/superGridWidth) * pow( (1 -(z)/superGridWidth),pSG-1) )

//  sigmaz = @(z)  ( (1-epsL)* q*( 1 - (1 -z/l).^p ).^(q-1) ).* ( (p/l)*(1 -z/l).^(p-1) );  % d(sigma)/dz    

    epsL=1e-4;
    pSG=4;
    qSG=4;
    // superGridWidth=.2;
    const real & superGridWidth = dbase.get<Real>("superGridWidth");

    if( 1==1 )
    {
      // -- check the formual for the derivative of eta

      Real x= superGridWidth/2;
      Real delta = pow( REAL_EPSILON,1./3.);
      Real etap = etaSG(x+delta);
      Real etam = etaSG(x-delta);
      Real etaPrime = (etap-etam)/(2.*delta);

      Real etaDeriv, err;
      etaDeriv = etaSGz(x); 

      err = ( etaDeriv - etaPrime)/etaPrime;
      printF(">> SUPER GRID: x=%12.4e, etaz = %12.6e (analytic), etaPrime = %12.6e (FD), rel-err=%9.2e\n",x,etaDeriv,etaPrime,err );

    }
   
    // useAbsorbingLayer(axis,grid) = true or false if we use an absorbing layer for this axis and grid 
    IntegerArray & useAbsorbingLayer = dbase.get<IntegerArray>("useAbsorbingLayer");
    useAbsorbingLayer.redim(3,numberOfComponentGrids);
    useAbsorbingLayer=false;
    
    for( int grid=0; grid<numberOfComponentGrids; grid++ )
    {
      MappedGrid & mg = cg[grid];
      const bool isRectangular=mg.isRectangular();
      const IntegerArray & dimension = mg.dimension();
      const IntegerArray & gid = mg.gridIndexRange();
      const IntegerArray & bc = mg.boundaryCondition();

      bool buildLayersThisGrid=false;
      // bool buildLayersThisAxis[3] ={false,false,false }; // 
      for( int axis=0; axis<numberOfDimensions; axis++)
      {
        for( int side=0; side<=1; side++ )
        {
          if( bc(side,axis)==absorbing )   
          {
            buildLayersThisGrid=true;
            // buildLayersThisAxis[axis]=true;
            useAbsorbingLayer(axis,grid)=true;
          }
        }
      }
      if( !buildLayersThisGrid ) continue;

      // // *** TEMP: For now build layers in all directions if any direction is needed -- need to optimize implementation
      // //  to allow only some directions 
      //        for( int axis=0; axis<numberOfDimensions; axis++)
      //        {
      //          buildLayersThisAxis[axis]=buildLayersThisGrid;
      //        }
      

      // assert( isRectangular );  // assume this for now 

      RealArray & etax = etaxSuperGrid[grid];
      RealArray & etay = etaySuperGrid[grid];
      RealArray & etaz = etazSuperGrid[grid];
      
      Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
      getIndex( mg.dimension(),I1,I2,I3 );          // all points including ghost points.
      OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
      bool ok = ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,I1,I2,I3,1);   
      if( !ok ) return 0;  // no points on this processor (NOTE: no communication should be done after this point)

      // For Cartesian grids we save the coefficients we need: 
      //   eta(:,0) = etax^2
      //   eta(:,1) = etaxx 
      //   u.xx = (r.x)^2 D+xD-x + r.xx Dzx 
      const int numEtaComponents = isRectangular ? 2 : 1;
      
      if( useAbsorbingLayer(0,grid) )
      {
        Range D1=maskLocal.dimension(0);
        etax.redim(D1,numEtaComponents);         // local dimensions of etax should match dimensions of uLocal (for call to advBA)
        etax(D1,0)=1.; 
        if( numEtaComponents==2 )
          etax(D1,1)=0.; 
      }
      if( useAbsorbingLayer(1,grid) )
      {
        Range D2=maskLocal.dimension(1);
        etay.redim(D2,numEtaComponents);
        etay(D2,0)=1.;
        if( numEtaComponents==2 )
          etay(D2,1)=0.;
      }
      if( useAbsorbingLayer(2,grid) )
      {
        Range D3=maskLocal.dimension(2);
        etaz.redim(D3,numEtaComponents);
        etaz(D3,0)=1.;
        if( numEtaComponents==2 )
          etaz(D3,1)=0.;
      }
      

      // // -- we optimize for Cartesian grids (we can avoid creating the vertex array)
      // if( !isRectangular )
      // {
      //   mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
      //   OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);

      //   OV_ABORT("buildSuperGrid: FINISH ME FOR CURVILINEAR GRIDS");
      // }
      
      // real dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
      // int iv0[3]={0,0,0}; //
      // int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];  // NOTE: iv[0]==i1, iv[1]==i2, iv[2]==i3
      // real xv[3]={0.,0.,0.};
      // if( isRectangular )
      // {
      //   mg.getRectangularGridParameters( dvx, xab );
      //   for( int dir=0; dir<mg.numberOfDimensions(); dir++ )
      //   {
      //     iv0[dir]=mg.gridIndexRange(0,dir);
      //     if( mg.isAllCellCentered() )
      //       xab[0][dir]+=.5*dvx[dir];  // offset for cell centered
      //   }
      // }
      // This macro defines the grid points for rectangular grids:
      // #undef XC
      // #define XC(iv,axis) (xab[0][axis]+dvx[axis]*(iv[axis]-iv0[axis]))
   

      // -- Parameterize by the unit square coordinates ----

      int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];  // NOTE: iv[0]==i1, iv[1]==i2, iv[2]==i3
      Real dr[3]={1.,1.,1.};
      for( int axis=0; axis<numberOfDimensions; axis++)
        dr[axis] = mg.gridSpacing(axis);

      #define RC(iv,axis) (iv[axis]*dr[axis])

      for( int axis=0; axis<numberOfDimensions; axis++)
      {
        if( useAbsorbingLayer(axis,grid) )
        {
          // ---- Build layer function for this axis ----
          RealArray & eta = axis==0 ? etax : axis==1 ? etay : etaz;
          for( int dir=0; dir<3; dir++ ){ iv[dir]=gid(0,dir); } // maybe not needed 

          for( int i=Iv[axis].getBase(); i<=Iv[axis].getBound(); i++ )
          {
            assert( i>=eta.getBase(0) && i<=eta.getBound(0) );

            iv[axis]=i;
            real x = RC(iv,axis);
            if( bc(0,axis) == absorbing ) 
            {
              // real z = xab[0][axis] +superGridWidth - x;
              real z = 0.0 +superGridWidth - x;
              if( z>superGridWidth )
              {
                eta(i,0)=epsL;
                if( isRectangular )
                {
                  eta(i,1) = 0.;          // etaz
                  eta(i,0) = eta(i,0)*eta(i,0); // eta^2 
                }
              }
              else if( z>0 )
              {
                eta(i,0)=etaSG(z);
                if( isRectangular )
                {
                  eta(i,1)  = - eta(i,0)*etaSGz(z);   // rxx = rx* (rx).r  : etaz note minus since defn of z 
                  eta(i,0)  =   eta(i,0)*eta(i,0);  // eta^2 
                }

              }

              // printF(" i=%d, x=%e, z=%e, eta=%e\n",i,x,z,eta(i));
            }
            if( bc(1,axis) ==absorbing )
            {

              // real z = x - ( xab[1][axis] -superGridWidth); 
              real z = x - ( 1.0 -superGridWidth); 
              if( z>superGridWidth )
              {
                eta(i,0)=epsL;
                if( isRectangular )
                {
                  eta(i,1) = 0.;                // etaz        
                  eta(i,0) = eta(i,0)*eta(i,0); // eta^2
                } 
              }
              else if( z>0 )
              {
                eta(i,0)=etaSG(z);
                if( isRectangular )
                {
                  eta(i,1)  = eta(i,0)*etaSGz(z);    // rxx = rx* (rx).r  : etaz
                  eta(i,0)  = eta(i,0)*eta(i,0);   // eta^2 = rx^2                    
                }

              }
            }
          }
        } // end for buildLayers
      } // end for axis
    
      if( !isRectangular )
      {
        // === CURVILINEAR GRID : SCALE METRICS =====
        printF("++++++ SUPERGRID : SCALE METRICS FOR GRID=%d ++++++\n",grid);

        mg.update(MappedGrid::THEinverseVertexDerivative );
        OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rxLocal);  
        #define RXLOCAL(i1,i2,i3,dir,axis) rxLocal(i1,i2,i3,(dir)+numberOfDimensions*(axis))      

        getIndex( mg.dimension(),I1,I2,I3); 
        for( int axis=0; axis<numberOfDimensions; axis++)
        {
          FOR_3D(i1,i2,i3,I1,I2,I3)
          {
            // D_x =   r_x * eta_r * D_eta
            RXLOCAL(i1,i2,i3,0,axis) *= etax(i1);
            RXLOCAL(i1,i2,i3,1,axis) *= etay(i2);
            if( numberOfDimensions==3 )
              RXLOCAL(i1,i2,i3,2,axis) *= etaz(i3);
          }
        }
 
        // OV_ABORT("buildSuperGrid: FINISH ME FOR CURVILINEAR GRIDS");
      }    

      bool plotLayers=false;

      if( plotLayers )
      {
        assert( useAbsorbingLayer(0,grid) );
        
        // assert( gip !=NULL );
        // GenericGraphicsInterface & gi = *gip;
        // bool plotOption=true;  // by default we plot interactively
        GL_GraphicsInterface & gi = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface());

        PlotStuffParameters psp;
        aString title=sPrintF("Super-grid functions: width=%g p=%g q=%g",superGridWidth,pSG,qSG);
        // int numFields = 1;
        // int numPoints = dimension(1,0)-dimension(0,0)+1;
        // RealArray fields(numPoints,numFields);

        aString names[]={"etax"};

        RealArray xv(I1);
        int axis=0;
        for( int i=dimension(0,axis); i<=dimension(1,axis); i++ )
        {
          iv[axis]=i;
          xv(i)=RC(iv,axis);
        }
         
        if( false )
        {
          ::display(xv,"xv","%8.2e ");
          ::display(etax,"etax","%8.2e ");
        }

        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);

        #ifndef USE_PPP
           PlotIt::plot(gi, xv, etax, title, "super-grid", names,psp );
           names[0]="etay";
           PlotIt::plot(gi, xv, etay, title, "super-grid", names,psp );
        #else
         printF("SuperGrid: FINISH PLOTTING IN PARALLEL\n");
        #endif

      }
      


    }  // end for grid

  }
  
  return 0;
}
