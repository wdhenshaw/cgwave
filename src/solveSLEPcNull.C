// -- NULL VERSION OF solveSLEPc 
#include "CgWaveHoltz.h" 
#include "CgWave.h" 
#include "Overture.h"
// #include "gridFunctionNorms.h"
// #include "display.h"


// ============================================================================
/// \brief Solve for eigenmodes using SLEPc
// ============================================================================
int CgWaveHoltz::
solveSLEPc(int argc,char **args)
{

  printF("ERROR: CgWaveHoltz has not been compiled with SLEPc.\n");
  printF("Use 'setenv SLEPC_DIR <location of SLEPc folder>' and rebuild.\n");
  return 0;
}