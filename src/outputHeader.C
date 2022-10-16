#include "CgWave.h"

// =======================================================================================
/// \brief Output the header banner with parameters and grid info.
// =======================================================================================
void CgWave::
outputHeader()
{
  int numberOfComponents=0;

  const int & np = dbase.get<int>("np");

  const real & c                            = dbase.get<real>("c");
  const real & cfl                          = dbase.get<real>("cfl");
  const real & tFinal                       = dbase.get<real>("tFinal");
  const real & tPlot                        = dbase.get<real>("tPlot");
      
  const int & upwind                        = dbase.get<int>("upwind");
  // const real & ad4                          = dbase.get<real>("ad4"); // coeff of the artificial dissipation.
  const int & dissipationFrequency          = dbase.get<int>("dissipationFrequency");
  const int & preComputeUpwindUt             = dbase.get<int>("preComputeUpwindUt");
  const AssignInterpolationNeighboursEnum & assignInterpNeighbours = 
                             dbase.get<AssignInterpolationNeighboursEnum>("assignInterpNeighbours");
  const int & implicitUpwind                = dbase.get<int>("implicitUpwind");
      
  const int & orderOfAccuracy               = dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime         = dbase.get<int>("orderOfAccuracyInTime");
  const real & dt                           = dbase.get<real>("dt");
  const int & modifiedEquationApproach      = dbase.get<int>("modifiedEquationApproach");
      
  int & addForcing                          = dbase.get<int>("addForcing");
  ForcingOptionEnum & forcingOption         = dbase.get<ForcingOptionEnum>("forcingOption");
  const IntegerArray & gridIsImplicit       = dbase.get<IntegerArray>("gridIsImplicit");
  const RealArray & bImp                    = dbase.get<RealArray>("bImp");
  const RealArray & cImp                    = dbase.get<RealArray>("cImp");
  const int & chooseImplicitTimeStepFromCFL = dbase.get<int>("chooseImplicitTimeStepFromCFL");

  const int & computeErrors                 = dbase.get<int>("computeErrors");
  const int & useKnownSolutionForFirstStep  = dbase.get<int>("useKnownSolutionForFirstStep"); 
	      
  const int & solveHelmholtz                = dbase.get<int>("solveHelmholtz");
  const int & computeTimeIntegral           = dbase.get<int>("computeTimeIntegral");
  const int & adjustOmega                   = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  const int & adjustHelmholtzForUpwinding   = dbase.get<int>("adjustHelmholtzForUpwinding");
  RealArray & dxMinMax                      = dbase.get<RealArray>("dxMinMax");

  TwilightZoneEnum & twilightZone = dbase.get<TwilightZoneEnum>("twilightZone");
  const int & degreeInSpace = dbase.get<int>("degreeInSpace");
  const int & degreeInTime =  dbase.get<int>("degreeInTime");

  const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");
  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

  const aString & nameOfGridFile= dbase.get<aString>("nameOfGridFile");
  real & numberOfGridPoints = dbase.get<real>("numberOfGridPoints");
  numberOfGridPoints=0.;

  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    // numberOfGridPoints+=cg[grid].mask().elementCount();  // old way

    // Oct 2, 2022 : just count interior and boundary points
    MappedGrid & mg = cg[grid];
    const IntegerArray & gid = mg.gridIndexRange();
    numberOfGridPoints += (gid(1,0)-gid(0,0)+1)*(gid(1,1)-gid(0,1)+1)*(gid(1,2)-gid(0,2)+1);

  }



  FILE *& logFile = dbase.get<FILE*>("logFile");

  for( int fileio=0; fileio<2; fileio++ )
  {
    FILE *file = fileio==0 ? logFile : stdout; 
    fPrintF(file,"\n"
	    "*********************************************************************************\n"
	    "           CgWave : Wave Equation Solver                    \n"
	    "           -----------------------------                  \n");

    fPrintF(file," tFinal=%f, dt=%9.3e, tPlot=%9.3e cfl=%3.2f\n",tFinal,dt,tPlot,cfl );
    fPrintF(file," timeSteppingMethod = %s\n",(timeSteppingMethod==explicitTimeStepping ? "explicit (modified equation)" :
                                               timeSteppingMethod==implicitTimeStepping ? "implicit" : "unknown") );
    fPrintF(file," modifiedEquationApproach = %s\n",
          (modifiedEquationApproach==0 ? "standard modified equation" : 
                                          "hierarchical modified equation" ));
    fPrintF(file," implicit time-stepping weights: cImp(-1,0)=%g, cImp(0,0)=%g, cImp(1,0)=%g (2nd-order term)\n",cImp(-1,0),cImp(0,0),cImp(1,0));
    fPrintF(file,"                               : cImp(-1,1)=%g, cImp(0,1)=%g, cImp(1,1)=%g (4th-order term)\n",cImp(-1,1),cImp(0,1),cImp(1,1));
    fPrintF(file," chooseImplicitTimeStepFromCFL=%d (1=choose implicit dt from cfl, 0=choose dt from dtMax)\n",chooseImplicitTimeStepFromCFL);
    fPrintF(file," orderOfAccuracy=%i, orderOfAccuracyInTime=%d \n",orderOfAccuracy, (orderOfAccuracyInTime==-1 ? orderOfAccuracy : orderOfAccuracyInTime) );

    fPrintF(file," upwind dissipation is %s, dissipationFrequency=%i\n",(upwind ? "on" : "off"),dissipationFrequency);
    fPrintF(file," upwind dissipation: preComputeUpwindUt=%i \n"
                 "                     true=precompute Ut in upwind dissipation,\n"
                 "                     false=compute Ut inline in Gauss-Seidel fashion)\n",preComputeUpwindUt);
    fPrintF(file," implicitUpwind=%d : if true, include upwinding in implicit matrix when implicit time-stepping.\n",implicitUpwind);
    fPrintF(file," assignInterpNeighbours = %s (for wider upwind stencil)\n",
        (assignInterpNeighbours==interpolateInterpNeighbours ? "interpolate" : 
         assignInterpNeighbours==extrapolateInterpNeighbours ? "extrapolate" :
                                                               "interpolate" )); // default 
    
    fPrintF(file," adjust Helmholtz for upwinding=%d\n",adjustHelmholtzForUpwinding);
    fPrintF(file," forcingOption=%s.\n",(forcingOption==noForcing           ? "noForcing"           :
                                         forcingOption==twilightZoneForcing ? "twilightZoneForcing" :
                                         forcingOption==userForcing         ? "userForcing"         :
                                         forcingOption==helmholtzForcing    ? "helmholtzForcing"    : 
                                                                              "unknown"));

    fPrintF(file," twilightZone = %s, degreeInSpace=%d, degreeInTime=%d\n",
            (twilightZone==polynomial ? "polynomial" : "trigonometric"), degreeInSpace,degreeInTime);
    fPrintF(file," useKnownSolutionForFirstStep = %d (use known solution for first step step, if possible)\n",
          useKnownSolutionForFirstStep);
    
    fPrintF(file," BC approach = %s. [useDefault|useOneSided|useCompatibility|useLocalCompatibility]\n",
                  ( bcApproach==defaultBoundaryConditionApproach        ? "useDefaultApproachForBCs"                :
                    bcApproach==useOneSidedBoundaryConditions           ? "useOneSidedBCs"                          :
                    bcApproach==useCompatibilityBoundaryConditions      ? "useCompatibilityBCs"                     :
                    bcApproach==useLocalCompatibilityBoundaryConditions ? "useLocalCompatibilityBoundaryConditions" : 
                                                                          "unknown" ));

    fPrintF(file," solveHelmholtz=%d (1= we are solving a Helmholtz problem).\n",solveHelmholtz);
    if( solveHelmholtz )
    {
      fPrintF(file," **** solveHelmholtz=true : Solving the Helmholtz problem, adjustOmega=%d (for discrete symbol of D+tD-t) ****\n",adjustOmega);
    }
      
    fPrintF(file," computeTimeIntegral=%d\n",computeTimeIntegral);
    fPrintF(file," computeErrors=%d\n",computeErrors);
    

    if( forcingOption==twilightZoneForcing )
      fPrintF(file," Twilightzone flow is on.");


    fPrintF(file," number of processors=%d\n",np);

    int maxNameLength=3;
    int grid;
    for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      maxNameLength=max( maxNameLength,cg[grid].getName().length());

    fPrintF(file,"\n");
    fPrintF(file," Grid: %s \n",(const char*)nameOfGridFile);
    aString blanks="                                                                           ";
    fPrintF(file,"               Grid Data\n"
	    "               ---------\n"
	    "grid     name%s  gridIndexRange(0:1,0:2)           gridPoints        hmx      hmn   time-stepping\n",
	    (const char *)blanks(0,min(maxNameLength-3,blanks.length()-1)));
    char buff[180];
    sPrintF(buff,"%%4i: %%%is   ([%%2i:%%5i],[%%2i:%%5i],[%%2i:%%5i])  %%12g   %%8.2e %%8.2e    %%s\n",maxNameLength);
    real maxMax=0.,maxMin=0.,minMin=REAL_MAX;
    // numberOfGridPoints=0.; 
    for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];
      real & hMin = dxMinMax(grid,0);
      real & hMax = dxMinMax(grid,1);
      
      maxMax=max(maxMax,hMax);
      maxMin=max(maxMin,hMin);
      minMin=min(minMin,hMin);
    
      // real numGridPoints = c.mask().elementCount();
      const IntegerArray & gid = mg.gridIndexRange();
      real numGridPoints = (gid(1,0)-gid(0,0)+1)*(gid(1,1)-gid(0,1)+1)*(gid(1,2)-gid(0,2)+1);
    
      // fPrintF(file,"%4i: %20s ([%2i:%5i],[%2i:%5i],[%2i:%5i])  %8i   %8.2e %8.2e \n",
      fPrintF(file,buff,grid, (const char *)cg[grid].getName(),
	      mg.gridIndexRange(Start,axis1),mg.gridIndexRange(End,axis1),
	      mg.gridIndexRange(Start,axis2),mg.gridIndexRange(End,axis2),
	      mg.gridIndexRange(Start,axis3),mg.gridIndexRange(End,axis3),
	      numGridPoints,hMax,hMin, (gridIsImplicit(grid) ? "implicit" : "explicit") );
      // numberOfGridPoints+=numGridPoints;
    }
    fPrintF(file," total number of grid points =%g, min(hmn)=%6.2e, max(hmn)=%6.2e, max(hmx)=%6.2e,  \n\n",
	    numberOfGridPoints,minMin,maxMin,maxMax);

    displayBoundaryConditions(file);

    fPrintF(file,"*********************************************************************************\n\n");
    
  }

  // -- Title line for check file ----
  FILE *& checkFile = dbase.get<FILE*>("checkFile");
  assert( checkFile != NULL );

  // Get the current date
  time_t *tp= new time_t;
  time(tp);
  // tm *ptm=localtime(tp);
  const char *dateString = ctime(tp);

//  fPrintF(checkFile,"# Check file for CgWave. date=%s",dateString);

  delete tp;

}


void CgWave::
displayBoundaryConditions(FILE *file /* = stdout */)
{
//===================================================================================
// /Description:
//   Print names for boundary conditions
//
// /file (input) : write to this file.
//===================================================================================
  assert( file!=NULL );

  const int & numberOfBCNames = dbase.get<int>("numberOfBCNames");
  aString *& bcNames = dbase.get<aString*>("bcNames");

  // const int & useSuperGrid = parameters.dbase.get<int>("useSuperGrid");

  int maxNameLength=3;
  int grid;
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    maxNameLength=max( maxNameLength,cg[grid].getName().length());

  char buff[80];
  sPrintF(buff," %%4i: %%%is     %%i    %%i    %%3i :",maxNameLength); // build a format string


  aString blanks="                                                                           ";
  fPrintF(file," grid   name%s side axis    boundary condition and name\n",
           (const char *)blanks(0,min(maxNameLength-3,blanks.length()-1)));
  fPrintF(file," ----   ----%s ---- ----    ---------------------------\n",
           (const char *)blanks(0,min(maxNameLength-3,blanks.length()-1)));
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    for( int axis=axis1; axis<cg.numberOfDimensions(); axis++ )
    for( int side=Start; side<=End; side++ )
    {
      int bc=cg[grid].boundaryCondition()(side,axis);
      fPrintF(file,buff,grid,(const char *)cg[grid].getName(),side,axis,bc);

      // if( bc==absorbing && useSuperGrid )
      // {
      //  fPrintF(file," super-grid absorbing layer.\n");
      // }
      if( bc > 0 && bc<numberOfBCNames)
        fPrintF(file," %s \n",(const char*)bcNames[bc]);
      else if( bc==0 )
        fPrintF(file," %s \n","none");
      else if( bc<0 )
        fPrintF(file," %s \n","periodic");
      else
        fPrintF(file," %s \n","unknown");
    }
  }
}
