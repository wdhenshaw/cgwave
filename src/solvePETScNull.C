// -- NULL VERSION OF CGWAVEHOLTZ PETSC SOLVER 
#include "CgWaveHoltz.h" 
#include "CgWave.h" 
#include "Overture.h"
#include "gridFunctionNorms.h"
#include "display.h"



// ============================================================================
/// \brief Solve for the Helholtz solution using PETSc
// ============================================================================
int CgWaveHoltz::
solvePETSc(int argc,char **args)
{
  printF("ERROR: CgWaveHoltz has not been compiled with PETSC.\n");
  printF("Use 'setenv usePETSc on' and rebuild.\n");


  return 0;
}
