// This file automatically generated from userDefinedKnownSolution.bC with bpp.
#include "CgWave.h"
#include "GenericGraphicsInterface.h"
#include "ParallelUtility.h"


#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++)                                       for(i2=I2Base; i2<=I2Bound; i2++)                                     for(i1=I1Base; i1<=I1Bound; i1++)

typedef ::real LocalReal;

// ==========================================================================================
/// \brief  Evaluate a user defined known solution.
///
/// \param numberOfTimeDerivatives (input) : evaluate this many time-derivatives of the solution.
///     Normally  numberOfTimeDerivatives=0, but it can be 1 when the known solution is used
//      to define boundary conditions 
// ==========================================================================================
int CgWave::
getUserDefinedKnownSolution(real t,  int grid, 
                                                        realArray & ua, const Index & I1a, const Index &I2a, const Index &I3a, 
                                                        int numberOfTimeDerivatives /* = 0 */ )
{
    const real & c= dbase.get<real>("c");
    const real & dt= dbase.get<real>("dt");

    if( false && t<= 2.*dt )
        printF("--CgWave-- getUserDefinedKnownSolution at t=%9.3e \n",t);

    MappedGrid & mg = cg[grid];
    const int numberOfDimensions = cg.numberOfDimensions();

    
    if( ! dbase.has_key("userDefinedKnownSolutionData") )
    {
        printF("--MX-- getUserDefinedKnownSolution:ERROR: sub-directory `userDefinedKnownSolutionData' not found!\n");
        OV_ABORT("error");
    }
    DataBase & db =  dbase.get<DataBase>("userDefinedKnownSolutionData");

    const aString & userKnownSolution = db.get<aString>("userKnownSolution");

    real *rpar = db.get<real[20]>("rpar");
    int *ipar = db.get<int[20]>("ipar");
    
    OV_GET_SERIAL_ARRAY(real,ua,uLocal);

    Index I1=I1a, I2=I2a, I3=I3a;
    bool ok = ParallelUtility::getLocalArrayBounds(ua,uLocal,I1,I2,I3,1);   
    if( !ok ) return 0;  // no points on this processor (NOTE: no communication should be done after this point)

  // -- we optimize for Cartesian grids (we can avoid creating the vertex array)
    const bool isRectangular=mg.isRectangular();
    if( !isRectangular )
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
    OV_GET_SERIAL_ARRAY_CONDITIONAL(real,mg.center(),xLocal,!isRectangular );

    real dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
    int iv0[3]={0,0,0}; //
    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];  // NOTE: iv[0]==i1, iv[1]==i2, iv[2]==i3
    real xv[3]={0.,0.,0.};
    if( isRectangular )
    {
        mg.getRectangularGridParameters( dvx, xab );
        for( int dir=0; dir<mg.numberOfDimensions(); dir++ )
        {
            iv0[dir]=mg.gridIndexRange(0,dir);
            if( mg.isAllCellCentered() )
                xab[0][dir]+=.5*dvx[dir];  // offset for cell centered
        }
    }
  // This macro defines the grid points for rectangular grids:
#undef XC
#define XC(iv,axis) (xab[0][axis]+dvx[axis]*(iv[axis]-iv0[axis]))

    
    assert( numberOfTimeDerivatives==0 );  // this option currently not used 

    if( userKnownSolution=="planeWave" )
    {
        const real amp = rpar[0];
        real kx  = rpar[1]*twoPi;
        real ky  = rpar[2]*twoPi;
        real kz  = rpar[3]*twoPi;

        real k;
        if( numberOfDimensions==2 )
            k = sqrt( kx*kx + ky*ky );
        else      
            k = sqrt( kx*kx + ky*ky +kz*kz );
            
        const real omega = c*k;
        
        real x,y,z;
        if( numberOfDimensions==2 )
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( !isRectangular )
                {
                    x= xLocal(i1,i2,i3,0);
                    y= xLocal(i1,i2,i3,1);
                }
                else
                {
                    x=XC(iv,0);
                    y=XC(iv,1);
                }
                    
                uLocal(i1,i2,i3,0) = amp*sin( kx*x+ky*y - omega*t );
            }
        }
        else
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( !isRectangular )
                {
                    x= xLocal(i1,i2,i3,0);
                    y= xLocal(i1,i2,i3,1);
                    z= xLocal(i1,i2,i3,2);
                }
                else
                {
                    x=XC(iv,0);
                    y=XC(iv,1);
                    z=XC(iv,2);
                }

                uLocal(i1,i2,i3,0) = amp*sin( kx*x+ky*y+kz*z - omega*t );

            }
        }
        
    }
    else if(userKnownSolution=="boxHelmholtz" )
    {

        const real omega = rpar[0];
        const real kx  = rpar[1]*twoPi;
        const real ky  = rpar[2]*twoPi;
        const real kz  = rpar[3]*twoPi;

        if( false && t<= 2.*dt )
            printF("userDefinedKnownSolution: eval boxHelmholtz: omega=%g, kx=%g, ky=%g, kz=%g at t=%9.3e\n",omega,kx,ky,kz,t);

    // printF(" I1=[%i,%i] I2=[%i,%i] I3=[%i,%i]\n",I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound());
    // printF(" uLocal=[%i,%i][%i,%i][%i,%i]\n",
    //     uLocal.getBase(0),uLocal.getBound(0),
    //     uLocal.getBase(1),uLocal.getBound(1),
    //     uLocal.getBase(2),uLocal.getBound(2));

        const real amp=1.;
        const real coswt = cos(omega*t);
            
        real x,y,z;
        if( numberOfDimensions==2 )
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( !isRectangular )
                {
                    x= xLocal(i1,i2,i3,0);
                    y= xLocal(i1,i2,i3,1);
                }
                else
                {
                    x=XC(iv,0);
                    y=XC(iv,1);
                }
                uLocal(i1,i2,i3,0) = amp*sin( kx*x )*sin( ky*y )*coswt;
            }
        }
        else
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( !isRectangular )
                {
                    x= xLocal(i1,i2,i3,0);
                    y= xLocal(i1,i2,i3,1);
                    z= xLocal(i1,i2,i3,2);
                }
                else
                {
                    x=XC(iv,0);
                    y=XC(iv,1);
                    z=XC(iv,2);
                }
                uLocal(i1,i2,i3,0) = amp*sin( kx*x )*sin( ky*y )*sin( kz*z )*coswt;

            }
        }

    }
    else
    {
        printF("getUserDefinedKnownSolution:ERROR: unknown value for userDefinedKnownSolution=%s\n",
                      (const char*)userKnownSolution);
        OV_ABORT("ERROR");
    }
    
    return 0;
}


int CgWave::
updateUserDefinedKnownSolution()
// ==========================================================================================
/// \brief This function is called to set the user defined known solution.
/// 
// ==========================================================================================
{
  // Make  dbase.get<real >("a") sub-directory in the data-base to store variables used here
    if( ! dbase.has_key("userDefinedKnownSolutionData") )
          dbase.put<DataBase>("userDefinedKnownSolutionData");
    DataBase & db =  dbase.get<DataBase>("userDefinedKnownSolutionData");

    if( !db.has_key("userKnownSolution") )
    {
        db.put<aString>("userKnownSolution");
        db.get<aString>("userKnownSolution")="unknownSolution";
        
        db.put<real[20]>("rpar");
        db.put<int[20]>("ipar");
    }
    aString & userKnownSolution = db.get<aString>("userKnownSolution");
    real *rpar = db.get<real[20]>("rpar");
    int *ipar = db.get<int[20]>("ipar");


    const aString menu[]=
        {
            "no known solution",
            "plane wave",
            "box helmholtz",
            "done",
            ""
        }; 

    gi.appendToTheDefaultPrompt("userDefinedKnownSolution>");
    aString answer;
    for( ;; ) 
    {

        int response=gi.getMenuItem(menu,answer,"Choose a known solution");
        
        if( answer=="done" || answer=="exit" )
        {
            break;
        }
        else if( answer=="no known solution" )
        {
            userKnownSolution="unknownSolution";
        }
        else if( answer=="plane wave" ) 
        {
            userKnownSolution="planeWave";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution depends on time
            
            gi.inputString(answer,"Enter amp,kx,ky,kz");
            sScanF(answer,"%e %e %e %e",&rpar[0],&rpar[1],&rpar[2],&rpar[3]);
            printF(" Setting amp=%g, [kx,ky,kz]=[%g,%g,%g]\n",rpar[0],rpar[1],rpar[2],rpar[3]);



      // Save parameters in dbase so we can look them up in bcOptWave
            const Real amp=rpar[0], kx=rpar[1]*twoPi, ky=rpar[2]*twoPi, kz=rpar[3]*twoPi; 
            const Real & c= dbase.get<real>("c");
            real k;
            if( cg.numberOfDimensions()==2 )
                k = sqrt( kx*kx + ky*ky );
            else      
                k = sqrt( kx*kx + ky*ky +kz*kz );
            const real omega = c*k;      

            dbase.put<Real>("ampPlaneWave")   = amp;
            dbase.put<Real>("kxPlaneWave")    = kx;
            dbase.put<Real>("kyPlaneWave")    = ky;
            dbase.put<Real>("kzPlaneWave")    = kz;
            dbase.put<Real>("omegaPlaneWave") = omega;

            
        }
        else if( answer=="box helmholtz" ) 
        {
            printF("----------------- box helmholtz -----------------\n"
                          " Define a time-periodic solution for a square or box that depends on a forcing term.\n");
            printF("   We solve :   utt = c^2 * Delta(u) - f(x)*cos(omega*t) \n");
            printF("                u=0 on the boundary                    ) \n");
            printF("The solution is of the form: cos(omega*t)*sin(kx*2*pi*x)*sin(ky*2*pi*y)\n");
            
            userKnownSolution="boxHelmholtz";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution depends on time
            
            gi.inputString(answer,"Enter omega,kx,ky,kz");
            sScanF(answer,"%e %e %e %e",&rpar[0],&rpar[1],&rpar[2],&rpar[3]);
            printF(" Setting omega=%g, [kx,ky,kz]=[%g,%g,%g]\n",rpar[0],rpar[1],rpar[2],rpar[3]);

            dbase.get<real>("omega") = rpar[0]; // define the Helmholtz omega for advance 

      // Save parameters in dbase so we can look them up in bcOptWave
            const Real omega=rpar[0], kx=rpar[1]*twoPi, ky=rpar[2]*twoPi, kz=rpar[3]*twoPi; 
            dbase.put<Real>("omegaBoxHelmholtz") = omega;
            dbase.put<Real>("kxBoxHelmholtz")    = kx;
            dbase.put<Real>("kyBoxHelmholtz")    = ky;
            dbase.put<Real>("kzBoxHelmholtz")    = kz;
            
        }
        else
        {
            printF("unknown response=[%s]\n",(const char*)answer);
            gi.stopReadingCommandFile();
        }
        
    }

    gi.unAppendTheDefaultPrompt();
    bool knownSolutionChosen = userKnownSolution!="unknownSolution";
    return knownSolutionChosen;
}


