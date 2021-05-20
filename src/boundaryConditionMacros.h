// =======================================================================================
// Macro: get the parameters for calling the optimized fortran BC routine
// =======================================================================================
#beginMacro getBcOptParameters(u,uLocal)   

  IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
  ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( u[grid],indexRangeLocal,dimLocal,bcLocal );

  const bool isRectangular=mg.isRectangular();
  real dx[3]={1.,1.,1.};
  if( isRectangular )
    mg.getDeltaX(dx);

  int assignKnownSolutionAtBoundaries = 0;  // changed below 

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
      knownSolutionOption=1;
      assignKnownSolutionAtBoundaries=1;
    }
    else if( userKnownSolution=="boxHelmholtz"  ) 
    {
      knownSolutionOption=2;
      assignKnownSolutionAtBoundaries=1;  // not needed for square or box but is needed for cic **fix me**
    }

  }

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
    assignKnownSolutionAtBoundaries, // ipar(10)
    knownSolutionOption,             // ipar(11)
    addForcing,                      // ipar(12)
    forcingOption,                   // ipar(13)
    useUpwindDissipation,            // ipar(14)
    numGhost,                        // ipar(15)
    assignBCForImplicit              // ipar(16)
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
    REAL_MIN           //  rpar( 9)
                };

  real *pu = uLocal.getDataPointer();
  int *pmask = maskLocal.getDataPointer();
  real temp, *pxy=&temp, *prsxy=&temp;
  if( !isRectangular )
  {
    #ifdef USE_PPP
     prsxy=mg.inverseVertexDerivative().getLocalArray().getDataPointer();
    #else
     prsxy=mg.inverseVertexDerivative().getDataPointer();
    #endif    
  }
  bool vertexNeeded = twilightZone || knownSolutionOption!=0;
  if( vertexNeeded )
  {
    mg.update(MappedGrid::THEvertex);
    #ifdef USE_PPP
     pxy=mg.vertex().getLocalArray().getDataPointer();
    #else
     pxy=mg.vertex().getDataPointer();
    #endif    
  }
#endMacro                                  
