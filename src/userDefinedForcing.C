// This file automatically generated from userDefinedForcing.bC with bpp.
#include "CgWave.h"
#include "GenericGraphicsInterface.h"
#include "ParallelUtility.h"



#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)


//==============================================================================================
/// \brief Evaluate the user defined forcing.
/// \details This function is called to actually evaluate the user defined forcing
///   The function setupUserDefinedForcing is first 
///   called to assign the option and parameters. Rewrite or add new options to 
///   this function and to setupUserDefinedForcing to supply your own forcing option.
///
/// \param f (input/output) : add to this forcing function
/// \param iparam[] (input) : holds some integer parameters
/// \param rparam[] (input) : holds some real parameters
/// 
// NOTES:
//   (1) The forcing function is added to the right-hand side of Maxwell's equations
//        in second order form:
//               E_tt = c^2 ( E_xx + E_yy + E_zz ) + F(x,y,z,t) 
//==============================================================================================

int CgWave:: 
userDefinedForcing( realArray & f, int iparam[], real rparam[] )
{
  // Look for the userDefinedForcing sub-directory in the data-base
    if( !dbase.has_key("userDefinedForcingData") )
    {
    // if the directory is not there then assume that there is no user defined forcing
        return 0;
    }
    DataBase & db = dbase.get<DataBase>("userDefinedForcingData");

    aString & option= db.get<aString>("option");

    if( option=="none" ) // No user defined forcing has been specified
      return 0;

    const real t =rparam[0];       // Add the forcing at this time.
    const real dt =rparam[1];      // Current time step.
    const int & grid = iparam[0];  // Here is the grid we are on 
    const int & current = iparam[1];

    const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");    
    const RealArray & frequencyArray  = dbase.get<RealArray>("frequencyArray");


  //  ---- NOTE: change sign of forcing for Helmholtz since we want to solve ----
  //       ( omega^2 I + c^2 Delta) w = f      
  // const int & solveHelmholtz = dbase.get<int>("solveHelmholtz");
    
    ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const Real fSign = forcingOption==helmholtzForcing ? -1.0 : 1.0;  
  // const Real fSign = 1.;
    
    if( true )
    {
        printF("userDefinedForcing: t=%9.3e option=[%s] fSign=%g\n",t,(const char*)option,fSign);
    }

    const int numberOfComponentGrids = cg.numberOfComponentGrids();
    const int numberOfDimensions = cg.numberOfDimensions();

    assert( grid>=0 && grid<numberOfComponentGrids );
    MappedGrid & mg = cg[grid];

  // Here is the current solution: 
  // realCompositeGridFunction & ucg = cgfields[current];
  // realMappedGridFunction & u = ucg[grid];
    
  // Access the local arrays on this processor:
  // OV_GET_SERIAL_ARRAY(real,u,uLocal);
    OV_GET_SERIAL_ARRAY(real,f,fLocal);

  // -- we optimize for Cartesian grids (we can avoid creating the vertex array)
    const bool isRectangular=mg.isRectangular();
    if( true || !isRectangular ) // *** FIX ME *****************************
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
    OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);

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

    
    Index I1,I2,I3;
    getIndex( mg.dimension(),I1,I2,I3 );          // all points including ghost points.
  // getIndex( mg.gridIndexRange(),I1,I2,I3 );  // boundary plus interior points.
  // restrict bounds to local processor, include parallel ghost pts:
    bool ok = ParallelUtility::getLocalArrayBounds(f,fLocal,I1,I2,I3,1);   
    if( !ok ) return 0;  // no points on this processor (NOTE: no communication should be done after this point)


    if( option=="gaussianSources" )
    {
        fLocal=0.;
    
        const IntegerArray & numberOfGaussianSources = db.get<IntegerArray>("numberOfGaussianSources");
        const RealArray & gaussianParameters = db.get<RealArray>("gaussianParameters");

    // Add the Gaussian source terms to fLocal
        for( int freq=0; freq<numberOfFrequencies; freq++ )
        {
            for( int m=0; m<numberOfGaussianSources(freq); m++ )
            {
                const real a    = gaussianParameters(0,m,freq);
                const real beta = gaussianParameters(1,m,freq); 
        // const real omega= gaussianParameters(2,m,freq); 
                const real p    = gaussianParameters(3,m,freq);
                const real x0   = gaussianParameters(4,m,freq); 
                const real y0   = gaussianParameters(5,m,freq); 
                const real z0   = gaussianParameters(6,m,freq);
                const real t0   = gaussianParameters(7,m,freq);

        // Use frequencyArray so we get any adjusted frequencies: --> doesn't matter if only evaluated at t=0
                Real omega = frequencyArray(freq);

                if( true || t0 < 2.*dt )
                    printF(">>> Gaussian source %i, freq=%d: setting a=%8.2e, beta=%8.2e, omega=%8.2e, p=%8.2e, x0=%8.2e, y0=%8.2e, "
                                  "z0=%8.2e, t0=%8.2e, t=%8.2e\n", m,freq,a,beta,omega,p,x0,y0,z0,t0,t);

        // const real cost=cos(2.*Pi*omega*(t-t0));
        // const real sint=sin(2.*Pi*omega*(t-t0));
                const real cost=cos(omega*(t-t0));
        // const real sint=sin(omega*(t-t0));

                if( mg.numberOfDimensions()==2 )
                {
          // --- 2D ---
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        real x= xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1);
                        real rSq = SQR(x-x0)+SQR(y-y0);
            // real g = a*cost*exp( -beta*pow( rSq, p ) );
                        real aExp = a*exp( -beta*pow( rSq, p ) );
            // real g = aExp*cost;
                        real g = aExp*cost;
                        real rPow = p==1. ? 1 :  pow(rSq,p-1.);

                        fLocal(i1,i2,i3,freq)+= -fSign*rPow*g;
                        
                    }
                }
                else
                {
          // -- 3D ---
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        real x= xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1), z=xLocal(i1,i2,i3,2);
                        real rSq =  SQR(x-x0)+SQR(y-y0)+SQR(z-z0);
            // real g = a*cost*exp( -beta*pow( rSq, p ) );
                        real aExp = a*exp( -beta*pow( rSq, p ) );
                        real g = aExp*cost;
                        real rPow = pow(rSq,p-1.);
                        
                        fLocal(i1,i2,i3,freq)+= -fSign*rPow*g;

                    }
                }
            }
        }
    }
    else if( option=="boxHelmholtz" )
    {

        const real & c = dbase.get<real>("c");

    // --- Get data from the userDefinedKnownSolution ---
        
        DataBase & db =  dbase.get<DataBase>("userDefinedKnownSolutionData");
        const aString & userKnownSolution = db.get<aString>("userKnownSolution");
        assert( userKnownSolution=="boxHelmholtz" || userKnownSolution=="computedHelmholtz" );
        real *rpar = db.get<real[20]>("rpar");
        int *ipar = db.get<int[20]>("ipar");
    
    // if( omega != frequencyArray(0) )
    // {
    //   printF("userDefinedForcing: eval boxHelmholtz FORCING: Error - omega=%g is not equal to frequencyArray(0)=%g\n",omega,frequencyArray);
    //   OV_ABORT("error");
    // }
        for( int freq=0; freq<numberOfFrequencies; freq++ )
        {
      // const Real omega = frequencyArray(freq); 
      // // These next lines must match between userDefinedKnownSolution, userDefinedForcing and bcOptWave.bf90: 
      // const Real kx  = (rpar[1]+freq)*twoPi;
      // const Real ky  = (rpar[2]+freq)*twoPi;
      // const Real kz  = (rpar[3]+freq)*twoPi;

            Real omega,kx,ky,kz; 
          // This macro is used in: 
          //    userDefinedKnownSolution.bC
          //    userDefinedForcing.bC
          //    solveHelmholtz.bC 
                    omega = frequencyArray(freq);
                    kx  = rpar[1]*twoPi*(freq*.5+1.);
                    ky  = rpar[2]*twoPi*(freq*.5+1.);
                    kz  = rpar[3]*twoPi*(freq*.5+1.);  
          // kx  = (rpar[1]+freq)*twoPi;
          // ky  = (rpar[2]+freq)*twoPi;
          // kz  = (rpar[3]+freq)*twoPi; 

            printF("userDefinedForcing: eval boxHelmholtz FORCING: numberOfFrequencies=%d, freq=%d, omega=%g, kx=%g, ky=%g, kz=%g fSign=%g at t=%9.3e\n",
                        numberOfFrequencies,freq,omega,kx,ky,kz,fSign,t);
            

            if( mg.numberOfDimensions()==2 )
            {
                const real amp = fSign*( -SQR(omega) + c*c*( kx*kx + ky*ky ) );
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    real x= xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1);
                      
                    fLocal(i1,i2,i3,freq) = amp*sin(kx*x)*sin(ky*y);
                }
            }
            else
            {
        // --- 3D ---
                const real amp = fSign*( -SQR(omega) + c*c*( kx*kx + ky*ky + kz*kz ) );
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    real x= xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1), z=xLocal(i1,i2,i3,2);
                    fLocal(i1,i2,i3,freq) = amp*sin(kx*x)*sin(ky*y)*sin(kz*z);
                }
            }
        }// end for freq 
    }

    else if( option=="polyPeriodic" )
    {

        const real & c = dbase.get<real>("c");

    // --- Get data from the userDefinedKnownSolution ---
        const real & omega        = dbase.get<Real>("omegaPolyPeriodic");
        const int & degreeInSpace = dbase.get<int>("degreeInSpacePolyPeriodic");
        const real & a0           = dbase.get<Real>("a0PolyPeriodic");
        const real & a1           = dbase.get<Real>("a1PolyPeriodic");
        const real & b1           = dbase.get<Real>("b1PolyPeriodic");
        const real & c1           = dbase.get<Real>("c1PolyPeriodic"); 

        if( true || t<= 2.*dt )
            printF("userDefinedForcing: eval polyPeriodic omega=%g, a0=%g, a1=%g, b1=%g, c1=%g at t=%9.3e\n",omega,a0,a1,b1,c1,t);       
        
        if( mg.numberOfDimensions()==2 )
        {
            const real amp = fSign*( -SQR(omega) ); 
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                real x= xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1);
                  
                fLocal(i1,i2,i3) = ( a0 + a1*x + b1*y )*amp;
            }
        }
        else
        {
      // --- 3D ---
            const real amp = fSign*( -SQR(omega) );
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                real x= xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1), z=xLocal(i1,i2,i3,2);
                fLocal(i1,i2,i3) = ( a0 + a1*x + b1*y + c1*z )*amp;
            }
        }
    }  

    else
    {
        printF("CgWave::userDefinedForcing:ERROR: unknown option =[%s]\n",(const char*)option);
        OV_ABORT("error");
    }

    return 0;
}




//==============================================================================================
/// \brief This function is used to choose a user defined forcing and input parameters etc.
/// \details This function is used to setup and define the forcing to use.
/// The function userDefinedForcing (above) is called to actually assign the forcing.
///  Rewrite or add new options to  this routine to supply your own forcing.
/// Choose the "user defined forcing" option to have this routine called.
///
//==============================================================================================
int CgWave::
setupUserDefinedForcing()
{
  // here is a menu of possible forcing options
    aString menu[]=  
    {
        "no forcing",
        "gaussian sources",
        "box Helmholtz",
        "poly periodic",
        "exit",
        ""
    };
    aString answer,answer2;
    char buff[100];
    gi.appendToTheDefaultPrompt(">user defined forcing");

  // Make a sub-directory in the data-base to store variables used in userDefinedForcing
    if( !dbase.has_key("userDefinedForcingData") )
        dbase.put<DataBase>("userDefinedForcingData");

    DataBase & db = dbase.get<DataBase>("userDefinedForcingData");

  // option = the name of the user defined forcing.
    if( !db.has_key("option") )
    { // create option variable in the data base.
        db.put<aString>("option");
        db.get<aString>("option")="none"; // default option
    }
    aString & option = db.get<aString>("option");

    const int numberOfComponentGrids = cg.numberOfComponentGrids();
    const int numberOfDimensions = cg.numberOfDimensions();
    const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");    

    for( ;; )
    {
        gi.getMenuItem(menu,answer,"enter an option");
        
        if( answer=="exit" || answer=="done" )
        {
            break;
        }
        else if( answer=="no forcing" )
        {
            option="none";
        }
        else if( answer=="gaussian sources" )

        {
      // define a Gaussian forcing
            option="gaussianSources";

      // We save parameters in the data base 

            if( !db.has_key("numberOfGaussianSources") ) db.put<IntegerArray>("numberOfGaussianSources");
            IntegerArray & numberOfGaussianSources = db.get<IntegerArray>("numberOfGaussianSources");

            numberOfGaussianSources.redim(numberOfFrequencies);
            numberOfGaussianSources=0; 
  

            if( numberOfDimensions==2 )
            {
                printF("The Gaussian source in 2D is of the form:\n"
                              " u(x,y,t) = a*cos(omega*(t-t0) )*exp( -beta*[ (x-x0)^2 + (y-y0)^2 ]^p )\n"
                    );
            }
            else
            {
                printF("The Gaussian source in 3D is of the form:\n"
                              " g(x,y,z,t) = a*cos( omega*(t-t0) )*exp( -beta*[ (x-x0)^2 + (y-y0)^2 + (z-z0)^2 ]^p )\n"
                    );
            }
            
            if( !db.has_key("gaussianParameters") )
                db.put<RealArray>("gaussianParameters");
            RealArray & gaussianParameters = db.get<RealArray>("gaussianParameters");

            const int maxGaussianSources=10;
            gaussianParameters.redim(8,maxGaussianSources,numberOfFrequencies);
            gaussianParameters=0.;

      // --------- Enter Gaussian sources for each frequency ---------
            for( int freq=0; freq<numberOfFrequencies; freq++ )
            {

                gi.inputString(answer2,sPrintF("Enter the number of Gaussian sources for frequency %d (default = 1)",freq));
                sScanF(answer2,"%i",&numberOfGaussianSources(freq));

                if( numberOfGaussianSources(freq) > maxGaussianSources)
                {
                    printF("Error: requested more than %d Gaussian sources. Fix me Bill!\n");
                    OV_ABORT("error");
                }

        // if( freq==0 )
        // {
        //   gaussianParameters.redim(8,numberOfGaussianSources(freq),numberOfFrequencies);
        //   gaussianParameters=0.;
        //   maxGaussianSources = numberOfGaussianSources(freq);
        // }
        // else if( numberOfGaussianSources(freq)>maxGaussianSources )
        // { // **** THIS IS BROKEN ****
        //   RealArray gp;
        //   gp = gaussianParameters; // save existing
        //   gaussianParameters.redim(8,numberOfGaussianSources(freq),numberOfFrequencies);
        //   gaussianParameters=0.; 
        //   gaussianParameters(????, Range(maxGaussianSources),Range(numberOfFrequencies))=gp; // copy previous
        //   maxGaussianSources = numberOfGaussianSources(freq);
        // }

                for( int m=0; m<numberOfGaussianSources(freq); m++ )
                {
                    real a=1., beta=10., omega=1., p=1., x0=0., y0=0., z0=0., t0=0.;
                    gi.inputString(answer2,sPrintF("Source %i for freq %d: Enter a,beta,omega,p,x0,y0,z0,t0",m,freq));
                    sScanF(answer2,"%e %e %e %e %e %e %e %e",&a,&beta,&omega,&p,&x0,&y0,&z0,&t0);

                    printF("Gaussian source %i, freq=%d: setting a=%8.2e, beta=%8.2e, omega=%8.2e, p=%8.2e, x0=%8.2e, y0=%8.2e, "
                                  "z0=%8.2e, t0=%8.2e\n", m,freq,a,beta,omega,p,x0,y0,z0,t0);

                    gaussianParameters(0,m,freq)=a; 
                    gaussianParameters(1,m,freq)=beta; 
                    gaussianParameters(2,m,freq)=omega; 
                    gaussianParameters(3,m,freq)=p;
                    gaussianParameters(4,m,freq)=x0; 
                    gaussianParameters(5,m,freq)=y0; 
                    gaussianParameters(6,m,freq)=z0;
                    gaussianParameters(7,m,freq)=t0;
                }
            }

        }
        else if( answer=="box Helmholtz" )
        {
            printF("----------------- FORCING FOR the box helmholtz Solution-----------------\n"
                          " Define a time-periodic solution for a square or box that depends on a forcing term.\n");
            printF("   We solve :   utt = c^2 * Delta(u) - f(x)*cos(omega*t) \n");
            printF("                u=0 on the boundary                    ) \n");
            printF("The solution is of the form: cos(omega*t)*sin(kx*2*pi*x)*sin(ky*2*pi*y)\n");

            option="boxHelmholtz";

        }
        else if( answer=="poly periodic" )
        {
            printF("----------------- polynomial in space and periodic in time -----------------\n");
            printF("   We solve :   utt = c^2 * Delta(u) - f(x)*cos(omega*t) \n");
            printF("The solution is of the form: (a0 + a1*x + a2*x^2 + ... b1*y + b2*y^2 + ... + c1*z + c2*z^2 + ...)*cos(omega*t)\n");

            option="polyPeriodic";

        }

        else 
        {
            printF("Maxwell::setupUserDefinedForcing:ERROR: unknown option =[%s]\n",(const char*)answer);
            gi.stopReadingCommandFile();
        }
        
    }
    
    gi.unAppendTheDefaultPrompt();
    return 0;
}

