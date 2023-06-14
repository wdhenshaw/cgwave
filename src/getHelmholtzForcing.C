// This file automatically generated from getHelmholtzForcing.bC with bpp.
#include "CgWave.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "OGFunction.h"

// -------- function prototypes for Fortran routines --------
#define bcOptWave EXTERN_C_NAME(bcoptwave)
extern "C"
{

void bcOptWave( const int&nd, 
                                const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                                const int&gridIndexRange, const int& dimRange, const int &isPeriodic, real&u, const int&mask,
                                const real&rsxy, const real&xy, real &uTemp1, real & uTemp2, 
                                const int&boundaryCondition, const real & frequencyArray, 
                                const DataBase *pdb, const int&ipar, const real&rpar, int&ierr );

}


// The getBcOptParameters macro is defined here:
// =======================================================================================
// Macro: get the parameters for calling the optimized fortran BC routine
// =======================================================================================

// =======================================================================================================
/// \brief Fill in the forcing (right-hand-side) for the direct Helmholtz solver (solveHelmHoltz)
/// \details Fill in the interior forcing and boundary conditions 
/// \param f (output) 
// =======================================================================================================
int CgWave::getHelmholtzForcing( realCompositeGridFunction & f  )
{
    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::numberOfProcessors());

    const int & debug = dbase.get<int>("debug");

    const int & orderOfAccuracy = dbase.get<int>("orderOfAccuracy");
    const Real & c              = dbase.get<real>("c");
    const real & dt             = dbase.get<real>("dt");
    const int & upwind                = dbase.get<int>("upwind");
  // const real & ad4            = dbase.get<real>("ad4"); // coeff of the artificial dissipation.
  // bool useUpwindDissipation   = ad4  > 0.;
    bool useUpwindDissipation   = upwind!=0;
    const int & solveHelmholtz  = dbase.get<int>("solveHelmholtz");

    const int & numberOfFrequencies     = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray    = dbase.get<RealArray>("frequencyArray");

    const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption"); 

    const int & applyKnownSolutionAtBoundaries = dbase.get<int>("applyKnownSolutionAtBoundaries"); // by default, do NOT apply known solution at boundaries

    const int & addForcing                  = dbase.get<int>("addForcing");
    const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const bool twilightZone = forcingOption==twilightZoneForcing; 

    const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");

    Range all; 
    f.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
    
    Index I1,I2,I3;

    if( forcingOption==helmholtzForcing )
    {
        printF("CgWave::getHelmholtzForcing: ***FILL IN RHS FOR DIRECT HELMHOLTZ SOLVER***\n");
        
        int current=0; // not used 
        real t=0.;

        int ipar[10];
        real rpar[10];
        ipar[1]=current;
        rpar[0]=t;
        rpar[1]=dt;
    // -- evaluate the forcing for a Helmholtz solve ---
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            ipar[0]=grid;
            userDefinedForcing( f[grid], ipar,rpar );
        }
        
    }
    else
    {
        f=0.;
    }

  // -------- Fill in boundary conditions for implicit and direct Helmholtz solver ---------

    Real t=0.; // this should not matter
    int numGhost = orderOfAccuracy/2;
    if( useUpwindDissipation ) numGhost++;

    const int assignBCForImplicit = 1; 

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
    // const IntegerArray & gid = mg.gridIndexRange();
        
        OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

    // get parameters for calling fortran
            IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
            ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( f[grid],indexRangeLocal,dimLocal,bcLocal );
            const bool isRectangular=mg.isRectangular();
            real dx[3]={1.,1.,1.};
            if( isRectangular )
                mg.getDeltaX(dx);
      // int assignKnownSolutionAtBoundaries = 0;  // changed below 
            DataBase *pdb = &dbase;
      // Real cfl1 = pdb->get<real>("cfl");
      // printF(" CFL from pdb: cfl1=%g\n",cfl1);
            int knownSolutionOption=0; // no known solution
            if( dbase.has_key("userDefinedKnownSolutionData") )
            {
                DataBase & db =  dbase.get<DataBase>("userDefinedKnownSolutionData");
                const aString & userKnownSolution = db.get<aString>("userKnownSolution");
                if( userKnownSolution=="planeWave"  )
                {
                    knownSolutionOption=1;                   // this number must match in bcOptWave.bf90
          // assignKnownSolutionAtBoundaries=1;
                }
                else if( userKnownSolution=="gaussianPlaneWave"  ) 
                {
                    knownSolutionOption=2;                   // this number must match in bcOptWave.bf90
          // assignKnownSolutionAtBoundaries=1;  
                }    
                else if( userKnownSolution=="boxHelmholtz"  ) 
                {
                    knownSolutionOption=3;                   // this number must match in bcOptWave.bf90
          // assignKnownSolutionAtBoundaries=1;  // not needed for square or box but is needed for cic **fix me**
                }
                else if( userKnownSolution=="polyPeriodic"  ) 
                {
                    knownSolutionOption=4;                   // this number must match in bcOptWave.bf90
          // assignKnownSolutionAtBoundaries=1;  
                } 
                else
                {
                    knownSolutionOption=1000;  // all other user defined known solutions
                }   
            }
            int addForcingBC= forcingOption==noForcing ? 0 : 1;
            if( applyKnownSolutionAtBoundaries )
                addForcingBC=1; 
            int gridType = isRectangular ? 0 : 1;
            int gridIsImplicit = 0; 
            int numberOfComponents = 1;          // for now CgWave only has a single component 
            int uc = 0;                          // first component
            int ipar[] = {
                uc                  ,            // ipar( 0)
                numberOfComponents  ,            // ipar( 1)
                grid                ,            // ipar( 2)
                gridType            ,            // ipar( 3)
                orderOfAccuracy     ,            // ipar( 4)
                gridIsImplicit      ,            // ipar( 5)
                twilightZone        ,            // ipar( 6)
                np                  ,            // ipar( 7)
                debug               ,            // ipar( 8)
                myid                ,            // ipar( 9)
                applyKnownSolutionAtBoundaries,  // ipar(10)
                knownSolutionOption,             // ipar(11)
                addForcingBC,                    // ipar(12)
                forcingOption,                   // ipar(13)
                useUpwindDissipation,            // ipar(14)
                numGhost,                        // ipar(15)
                assignBCForImplicit,             // ipar(16)
                bcApproach,                      // ipar(17)
                numberOfFrequencies              // ipar(18)
                                      };
            real rpar[] = {
                t                , //  rpar( 0)
                dt               , //  rpar( 1)
                dx[0]            , //  rpar( 2)
                dx[1]            , //  rpar( 3)
                dx[2]            , //  rpar( 4)
                mg.gridSpacing(0), //  rpar( 5)
                mg.gridSpacing(1), //  rpar( 6)
                mg.gridSpacing(2), //  rpar( 7)
                (real &)(dbase.get<OGFunction* >("tz")) ,        //  rpar( 8) ! pointer for exact solution -- new : 110311 
                REAL_MIN,         //  rpar( 9)
                c                 //  rpar(10)
                                        };
            real *pu = fLocal.getDataPointer();
            int *pmask = maskLocal.getDataPointer();
            real temp, *pxy=&temp, *prsxy=&temp;
            if( !isRectangular )
            {
                mg.update(MappedGrid::THEinverseVertexDerivative); 
                #ifdef USE_PPP
                  prsxy=mg.inverseVertexDerivative().getLocalArray().getDataPointer();
                #else
                  prsxy=mg.inverseVertexDerivative().getDataPointer();
                #endif  
            }
            bool vertexNeeded = twilightZone || knownSolutionOption!=0;
            if( vertexNeeded )
            {
                mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter ); 
                #ifdef USE_PPP
                  pxy=mg.vertex().getLocalArray().getDataPointer();
                #else
                  pxy=mg.vertex().getDataPointer();
                #endif    
            }
            real *puTemp1=&temp, *puTemp2=&temp;
            if( orderOfAccuracy==4 && bcApproach==useCompatibilityBoundaryConditions )
            {
        // WE currently need work space for bcOptWave and standard CBC at order four ******************FIX ME: just make a stencil ****
                if( !dbase.has_key("uTemp1") )
                {
                    RealArray & uTemp1 = dbase.put<RealArray>("uTemp1");
                    RealArray & uTemp2 = dbase.put<RealArray>("uTemp2");
          // -- find the grid with physical boundaries with the most points ---
                    int numElements=0;
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cg[grid];
                        const IntegerArray & boundaryCondition = mg.boundaryCondition();
                        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                        if( max(boundaryCondition)>0 )
                            numElements = max( numElements, maskLocal.elementCount() );
                    }
                    printF(">>> INFO CgWave::applyBC: allocate uTemp1 and utemp2 for order 4 CBC: numElements=%d\n",numElements);
                    if( numElements>0 )
                    {
                        uTemp1.redim(numElements);
                        uTemp2.redim(numElements);
                    }
                }
                RealArray & uTemp1 = dbase.get<RealArray>("uTemp1");
                RealArray & uTemp2 = dbase.get<RealArray>("uTemp2");
                puTemp1 = uTemp1.getDataPointer();
                puTemp2 = uTemp2.getDataPointer();
            }

        int ierr=0;
        bcOptWave(mg.numberOfDimensions(),
                            fLocal.getBase(0),fLocal.getBound(0),fLocal.getBase(1),fLocal.getBound(1),
                            fLocal.getBase(2),fLocal.getBound(2),
                            indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),
                            *pu, *pmask, *prsxy, *pxy, *puTemp1, *puTemp2,
                            bcLocal(0,0), frequencyArray(0), 
                            pdb, ipar[0],rpar[0], ierr );

    // // ...swap periodic edges 
    // u[grid].periodicUpdate();
    // u[grid].updateGhostBoundaries();
        

    } // end for grid 


    return 0;
}