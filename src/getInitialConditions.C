#include "CgWave.h"
#include "display.h"
#include "ParallelUtility.h"
#include "OGPolyFunction.h"
#include "OGTrigFunction.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  

#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )


// ======================================================================================================
/// \brief Assign the initial conditions
/// \param current (input)  : assign the grid function u[current]                                 
/// \param getTimeDerivative (input) : if true, return the time derivative of the initial conditions                                  
// ======================================================================================================
void CgWave::
getInitialConditions( int current, real t, bool getTimeDerivative /* = false */ )
{
  real cpu0 =getCPU();
  

  // const int myid=Communication_Manager::My_Process_Number;

  real & dt = dbase.get<real>("dt");
  printF("++++++++++++ getInitialConditions current=%d, t=%9.3e, dt=%9.3e ++++++++++++++ \n",current,t,dt);
  
  // enum InitialConditionOptionEnum
  // {
  //   smoothPulse,
  //   pulse
  // };

  const int numberOfDimensions = cg.numberOfDimensions();

  const InitialConditionOptionEnum & initialConditionOption = dbase.get<InitialConditionOptionEnum>("initialConditionOption");

  // InitialConditionOptionEnum option=smoothPulse;

  const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");

  const int & addForcing = dbase.get<int>("addForcing");
  const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");


  realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
  real & c         = dbase.get<real>("c");

  // Remove these: *wdh* Oct 31, 2021
  // real & beta = dbase.get<real>("beta");
  // real & x0   = dbase.get<real>("x0");
  // real & y0   = dbase.get<real>("y0");
  // real & z0   = dbase.get<real>("z0");

  Index I1,I2,I3;

  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    
    int cur = current;

    OV_GET_SERIAL_ARRAY(real,u[cur ][grid],ucLocal);
    // OV_GET_SERIAL_ARRAY(real,u[prev][grid],upLocal);

    getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points.
    const bool isRectangular=cg[grid].isRectangular();

    if( initialConditionOption == zeroInitialCondition )
    {
      ucLocal = 0.;
      // upLocal = 0.;
    }
    else if( addForcing && forcingOption==twilightZoneForcing )
    {
      // ---- Twilight Zone Initial conditions ----

      printF("$$$$ Get TZ initial conditions: t=%9.3e, getTimeDerivative=%d $$$\n",t,getTimeDerivative);

      assert( dbase.get<OGFunction*>("tz")!=NULL );
      OGFunction & e = *dbase.get<OGFunction*>("tz");

      mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
      OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);

      const int includeGhost=1;
      bool ok=ParallelUtility::getLocalArrayBounds(u[0][grid],ucLocal,I1,I2,I3,includeGhost);
      if( ok )
      {
        int numberOfComponents=1;
        Range C=numberOfComponents;
        int isRectangular=0;
        int ntd = getTimeDerivative ? 1 : 0; // number of time derivatives
        printF("$$$$ Get TZ initial conditions: ntd=%d $$$\n",ntd);
        e.gd( ucLocal ,xLocal,numberOfDimensions,isRectangular,ntd,0,0,0,I1,I2,I3,C,t);

      }
      
    }
    else if( knownSolutionOption == "userDefinedKnownSolution" )
    {
      // -- User defined known solution ---
      int numberOfTimeDerivatives = getTimeDerivative ? 1 : 0;
      getUserDefinedKnownSolution( t, grid, u[cur ][grid], I1,I2,I3, numberOfTimeDerivatives );
        
    }
    else if( initialConditionOption == pulseInitialCondition  )
    {

      real beta =20, x0=0., y0=0.; // do this for now : *wdh* Oct 31, 2021

      // Pulse parameters:
      real alpha=beta; // 200.;
      real pulsePow=2.; // 20
      real a0=3.;
      real xPulse=x0;
      real yPulse=y0;



      #define U0(x,y,t) exp( - alpha*( SQR((x)-(xPulse-c*t)) + SQR((y)-yPulse) ) )
      // time derivative: (check me)
      #define U1(x,y,t) exp( - alpha*( SQR((x)-(xPulse-c*t)) + SQR((y)-yPulse) ) )*( (-2.*c*alpha)*( ((x)-(xPulse-c*t)) ) )

      // restrict the bounds (I1,I2,i3) to the local array bounds (including parallel ghost pts):
      const int includeGhost=1;
      bool ok=ParallelUtility::getLocalArrayBounds(u[0][grid],ucLocal,I1,I2,I3,includeGhost);

      if( isRectangular )
      {
        // for a rectangular grid we avoid building the array of verticies.
        // we assign the initial conditions with C-style loops

        if( !ok ) continue;  // nothing to do on this processor

        real dx[3]={0.,0.,0.}, xab[2][3]={0.,0.,0.,0.,0.,0.};
        if( cg[grid].isRectangular() )
          cg[grid].getRectangularGridParameters( dx, xab );

        const real xa=xab[0][0], dx0=dx[0];
        const real ya=xab[0][1], dy0=dx[1];
        const real za=xab[0][2], dz0=dx[2];

        const int i0a=cg[grid].gridIndexRange(0,0);
        const int i1a=cg[grid].gridIndexRange(0,1);
        const int i2a=cg[grid].gridIndexRange(0,2);

#define VERTEX0(i0,i1,i2) xa+dx0*(i0-i0a)
#define VERTEX1(i0,i1,i2) ya+dy0*(i1-i1a)
#define VERTEX2(i0,i1,i2) za+dz0*(i2-i2a)

        // Here we grab a pointer to the data of the array so we can index it as a C-array
        // real *upm= upLocal.Array_Descriptor.Array_View_Pointer3;
        real *up = ucLocal.Array_Descriptor.Array_View_Pointer3;
        const int uDim0=ucLocal.getRawDataSize(0);
        const int uDim1=ucLocal.getRawDataSize(1);
        const int d1=uDim0, d2=d1*uDim1; 
#define U(i0,i1,i2) up[(i0)+(i1)*d1+(i2)*d2]
// #define UM(i0,i1,i2) upm[(i0)+(i1)*d1+(i2)*d2]

        int i1,i2,i3;
        FOR_3(i1,i2,i3,I1,I2,I3) // loop over all points
        {
          // UM(i1,i2,i3)=U0(VERTEX0(i1,i2,i3),VERTEX1(i1,i2,i3),-dt);
          if( getTimeDerivative )
            U(i1,i2,i3) =U1(VERTEX0(i1,i2,i3),VERTEX1(i1,i2,i3),t);
          else
            U(i1,i2,i3) =U0(VERTEX0(i1,i2,i3),VERTEX1(i1,i2,i3),t);

        }
        
#undef VERTEX0
#undef VERTEX1
#undef VERTEX2
#undef U
#undef UM
      }
      else
      {
        cg[grid].update(MappedGrid::THEvertex | MappedGrid::THEcenter );  // build the array of vertices
        realArray & vertex = cg[grid].vertex();

        if( !ok ) continue;  // nothing to do on this processor

        // real *upm= upLocal.Array_Descriptor.Array_View_Pointer3;
        real *up = ucLocal.Array_Descriptor.Array_View_Pointer3;
        const int uDim0=ucLocal.getRawDataSize(0);
        const int uDim1=ucLocal.getRawDataSize(1);
        const int d1=uDim0, d2=d1*uDim1; 
#define U(i0,i1,i2) up[(i0)+(i1)*d1+(i2)*d2]
// #define UM(i0,i1,i2) upm[(i0)+(i1)*d1+(i2)*d2]

        OV_GET_SERIAL_ARRAY(real,vertex,vertexLocal);
        const real *vertexp = vertexLocal.Array_Descriptor.Array_View_Pointer3;
        const int vertexDim0=vertexLocal.getRawDataSize(0);
        const int vertexDim1=vertexLocal.getRawDataSize(1);
        const int vertexDim2=vertexLocal.getRawDataSize(2);
#define VERTEX(i0,i1,i2,i3) vertexp[i0+vertexDim0*(i1+vertexDim1*(i2+vertexDim2*(i3)))]

        int i1,i2,i3;
        FOR_3(i1,i2,i3,I1,I2,I3) // loop over all points
        {
          // UM(i1,i2,i3)=U0(VERTEX(i1,i2,i3,0),VERTEX(i1,i2,i3,1),-dt);
          if( getTimeDerivative )
            U(i1,i2,i3) =U1(VERTEX(i1,i2,i3,0),VERTEX(i1,i2,i3,1),t);
          else
            U(i1,i2,i3) =U0(VERTEX(i1,i2,i3,0),VERTEX(i1,i2,i3,1),t);

        }

      }
#undef U
#undef UM
#undef VERTEX

    }
    else
    {
      printF("CgWave::getInitialConditions: ERROR: unknown initialConditionOption=%d\n",(int)initialConditionOption);
      OV_ABORT("error");

    }
    // if( !plotOption ) 
    //   cg[grid].destroy(MappedGrid::THEvertex);  // vertices are no nolonger needed.

  }

  printF("CgWave: done initial conditions\n");
  timing(timeForInitialConditions) += getCPU()- cpu0;
}

