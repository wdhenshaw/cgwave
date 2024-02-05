
#include "CgWave.h"
#include "OGFunction.h"
#include "ParallelUtility.h"
#include "gridFunctionNorms.h"
#include "Integrate.h"
#include "CompositeGridOperators.h"   

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  
#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

// ======================================================================================================
/// \brief Compute an approximation to the energy
// ======================================================================================================
Real CgWave::
getEnergyNorm( int cur, Real t )
{


  Real energyNorm=0.;

  const Real & dt = dbase.get<real>("dt");
  const Real & c  = dbase.get<real>("c");

  const int numberOfDimensions = cg.numberOfDimensions();


  const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
  const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
  const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

  realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
  CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

  realCompositeGridFunction & uc = u[cur];     // current time 
  realCompositeGridFunction & up = u[prev];    // previous time

  if( !dbase.has_key("integrate") )
  {
    Integrate & integrate = dbase.put<Integrate>("integrate");
    integrate.updateToMatchGrid(cg);
  }
  Integrate & integrate = dbase.get<Integrate>("integrate");

  RealCompositeGridFunction w(cg);
  Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];

  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    const bool isRectangular = mg.isRectangular();

    Real dx[3]={1.,1.,1.};
    Real dr[3]={1.,1.,1.};
    if( isRectangular )
    {
      mg.getDeltaX(dx);    
    }
    else
    {
      mg.update(MappedGrid::THEinverseVertexDerivative);
      for( int axis=0; axis<numberOfDimensions; axis++ )
        dr[axis]=mg.gridSpacing(axis);
    }
    OV_GET_SERIAL_ARRAY(Real,mg.inverseVertexDerivative(),rxLocal); 
    #define RX(i1,i2,i3,m1,m2) rxLocal(i1,i2,i3,(m1)+numberOfDimensions*(m2))     

    int extra=0; // do this for now
    getIndex(mg.gridIndexRange(),I1,I2,I3,extra); 
    for( int axis=0; axis<numberOfDimensions; axis++ )
    {
      Iv[axis] = Range( Iv[axis].getBase()+1,Iv[axis].getBound() ); // we use D- below so shift range
    }

    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
    OV_GET_SERIAL_ARRAY(Real,w[grid],wLocal);      
    OV_GET_SERIAL_ARRAY(Real,uc[grid],ucLocal);      
    OV_GET_SERIAL_ARRAY(Real,up[grid],upLocal);      

    // operators[grid].derivative(MappedGridOperators::laplacianOperator,qLocal,lapLocal,I1,I2,I3,ia);
    Real cSq = c*c;
    wLocal=0.;


    int i1,i2,i3;
    Real ut=0.,ux,uy,uz=0., ur1,ur2,ur3; 
    FOR_3D(i1,i2,i3,I1,I2,I3)
    {
      if( maskLocal(i1,i2,i3) > 0 )
      {   

        if( t>0. ) // fix me for t==0 
          ut = (ucLocal(i1,i2,i3)-upLocal(i1,i2,i3))/dt;   // D-t 

        if( isRectangular )
        {
          ux = (ucLocal(i1,i2,i3)-ucLocal(i1-1,i2,i3))/dx[0];    // D-x 
          uy = (ucLocal(i1,i2,i3)-ucLocal(i1,i2-1,i3))/dx[1];    // D-y 
          if( numberOfDimensions==3 )
            uz = (ucLocal(i1,i2,i3)-ucLocal(i1,i2,i3-1))/dx[2];  // D-z        
        }
        else
        {
          ur1 = (ucLocal(i1,i2,i3)-ucLocal(i1-1,i2,i3))/dr[0];
          ur2 = (ucLocal(i1,i2,i3)-ucLocal(i1,i2-1,i3))/dr[1];
          if( numberOfDimensions==2 )
          {
            ux =  RX(i1,i2,i3,0,0)*ur1 + RX(i1,i2,i3,1,0)*ur2;
            uy =  RX(i1,i2,i3,0,1)*ur1 + RX(i1,i2,i3,1,1)*ur2;                 
          }
          else
          {
            ur3 = (ucLocal(i1,i2,i3)-ucLocal(i1,i2,i3-1))/dr[2];
            ux = RX(i1,i2,i3,0,0)*ur1 + RX(i1,i2,i3,1,0)*ur2 + RX(i1,i2,i3,2,0)*ur3;   
            uy = RX(i1,i2,i3,0,1)*ur1 + RX(i1,i2,i3,1,1)*ur2 + RX(i1,i2,i3,2,1)*ur3;  
            uz = RX(i1,i2,i3,0,2)*ur1 + RX(i1,i2,i3,1,2)*ur2 + RX(i1,i2,i3,2,2)*ur3;                                           
          }
        }

        wLocal(i1,i2,i3) = .5*( ut*ut + cSq*( ux*ux + uy*uy + uz*uz ) );

        // Real ucx = (ucLocal(i1,i2,i3)-ucLocal(i1-1,i2,i3))/dx[0];
        // Real upx = (upLocal(i1,i2,i3)-ucLocal(i1-1,i2,i3))/dx[0];

        // Real ucy = (ucLocal(i1,i2,i3)-ucLocal(i1,i2-1,i3))/dx[1];
        // Real upy = (upLocal(i1,i2,i3)-upLocal(i1,i2-1,i3))/dx[1];
        // if( numberOfDimensions==3 )
        // {
        //   ucz = (ucLocal(i1,i2,i3)-ucLocal(i1,i2,i3-1))/dx[2];
        //   upz = (upLocal(i1,i2,i3)-upLocal(i1,i2,i3-1))/dx[2];
        // }
        // // discrete energy: see max2006b
        // wLocal(i1,i2,i3) = .5*( ut*ut + cSq*( ucx*upx + ucy*upy + ucz*upz ) );

      }
    }
  } // end for grid


  energyNorm = integrate.volumeIntegral(w); 
  // energyNorm = sqrt(energyNorm);

  //   printF("$$$$$$$$$$ getEnergyNorm: t=%9.2e: energyNorm=%9.2e\n",t,energyNorm);
  
  return energyNorm;

}