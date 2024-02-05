// =======================================================================================
// Macro: get the parameters for calling the optimized fortran BC routine
// =======================================================================================
#beginMacro getBcOptParameters(u,uLocal,unLocal)   

  IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
  ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( u[grid],indexRangeLocal,dimLocal,bcLocal );

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
  // int gridIsImplicit = 0; 
  int numberOfComponents = 1;          // for now CgWave only has a single component 
  int uc = 0;                          // first component
  int ipar[] = {
    uc                  ,            // ipar( 0)
    numberOfComponents  ,            // ipar( 1)
    grid                ,            // ipar( 2)
    gridType            ,            // ipar( 3)
    orderOfAccuracy     ,            // ipar( 4)
    gridIsImplicit(grid),            // ipar( 5)  // added Nov 22, 2023
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
    numberOfFrequencies,             // ipar(18)
    assignCornerGhostPoints          // ipar(19)
               };

  Real cEM2 = c;
  if( solveHelmholtz ) 
  {
    // Adjust c for the EM2 absorbing BC to account for time-discretization errors
    //   D+t (Dx ) w + A+( ... )
    if( frequencyArraySave(0)*dt>0 )
    {
      cEM2 = c*tan(frequencyArray(0)*dt/2.)/(frequencyArraySave(0)*dt/2.);
      // printF("\n XXXXXXX cEM2=%e XXXXXXX\n\n",cEM2);
    }
    else
    {
      if( frequencyArraySave(0)==0 )
        printF("WARNING: bcOpt: frequencyArraySave(0)=%12.4e. NOT ADJUSTING c for EM2 absorbing BC\n",frequencyArraySave(0));
      if( dt<=0 )
        printF("WARNING: bcOpt: dt<= 0 ! dt=%12.4e. NOT ADJUSTING c for EM2 absorbing BC\n",dt);
    }
  }           
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
    c,                //  rpar(10)
    cEM2              //  rpar(11)
                };

  real *pu = uLocal.getDataPointer();
  real *pun = unLocal.getDataPointer();
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

#endMacro                                  
