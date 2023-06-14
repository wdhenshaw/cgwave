#
# CgWave: Composite Grid Wave Equation Solver 
# CgWaveHoltz: Composite Grid Helmholtz Solver using the WaveHoltz algorithm of Daniel Appelo 
#
# NOTE: To compile optimized:
#   setenv COMPILE [opt|dbg]
#


include ${Overture}/make.options

# here is Arpack : 
LIB_ARPACK = -Wl,-rpath,/home/henshw/software/arpack-ng/lib -L/home/henshw/software/arpack-ng/lib -larpack 
LIB_ARPACK =


usePETSc := on
# usePETSc := off
ifeq ($(usePETSc),on)
  usePETSc   = on
  petscSolver = obj/solvePETSc.o obj/solveSLEPc.o

  ifeq ($(OV_PARALLEL),parallel) 
    OGES_PETSC = buildEquationSolvers.o PETScSolver.o
  else
    OGES_PETSC = buildEquationSolvers.o PETScEquationSolver.o
  endif

  # PETSC_INCLUDE = -I$(PETSC_DIR)/include  -I$(PETSC_DIR)/$(PETSC_ARCH)/include -DOVERTURE_USE_PETSC -I$(PETSC_LIB)/include -I$(PETSC_DIR)/include/mpiuni
  PETSC_INCLUDE = -I$(PETSC_DIR)/include  -I$(PETSC_DIR)/$(PETSC_ARCH)/include -DOVERTURE_USE_PETSC -I$(PETSC_LIB)/include -I$(PETSC_DIR)/include/petsc/mpiuni
  PETSC_LIBS = -Wl,-rpath,$(PETSC_LIB) -L$(PETSC_LIB) -lpetsc

  SLEPC_INCLUDE = -I$(SLEPC_DIR) -I$(SLEPC_DIR)/linux-gnu-opt/include -I$(SLEPC_DIR)/include 

  SLEPC_LIBS = $(LIB_ARPACK) -Wl,-rpath,$(SLEPC_DIR)/linux-gnu-opt/lib -L$(SLEPC_DIR)/linux-gnu-opt/lib -lslepc

else
  usePETSc   = off
  petscSolver = obj/solvePETScNull.o

  OGES_PETSC = 
  PETSC_INCLUDE = 
  PETSC_LIBS = 

  SLEPC_INCLUDE =
  SLEPC_LIBS =
endif

# optimization flag:
#  -Ofast = -O3 + disregards strict standard compliance 
# OPTFLAG = -O
OPTFLAG = -O3
# OPTFLAG = -Ofast


CCFLAGS = $(OV_CXX_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include -I$(OpenGL)/include $(USE_PPP_FLAG)

CCFLAGS += $(SLEPC_INCLUDE)
CCFLAGS += $(PETSC_INCLUDE)

# for rtg1 which has an older compiler (Processor.h)
CCFLAGS += -std=c++11

CFLAGS = $(OV_CC_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include

# Some C++ files we compile optimized by default
ifeq ($(COMPILE),dbg)
  CCFLAGSO = $(OV_CXX_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include -I$(OpenGL)/include $(USE_PPP_FLAG) -g -w
else
  CCFLAGSO = $(OV_CXX_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include -I$(OpenGL)/include $(USE_PPP_FLAG) $(OPTFLAG)
endif

# Fortran flags, not optimzed unless COMPILE set to opt
FFLAGS  = $(OV_FORTRAN_FLAGS) $(OV_AUTO_DOUBLE_FLAGS) 

# Fortran flags, optimized by default 
FFLAGSO = $(OV_FORTRAN_FLAGS) $(OV_AUTO_DOUBLE_FLAGS)
#
# Fortran flags for debug always
FFLAGSG = $(OV_FORTRAN_FLAGS) $(OV_AUTO_DOUBLE_FLAGS) -g 

# no auto-double:
FFLAGSSO = $(OV_FORTRAN_FLAGS)

ifeq ($(COMPILE),opt)
  CCFLAGS += $(OPTFLAG)
  FFLAGS  += $(OPTFLAG)
  FFLAGSSO += $(OPTFLAG)
else
	ifeq ($(COMPILE),dbg)
    # debug:
    CCFLAGS += -w -g -finit-real=snan
    FFLAGS  += -g
    FFLAGSO += -g
    FFLAGSSO += g   
  else
    # default case:
    CCFLAGS += -w -g 
    FFLAGS  += -g 
    FFLAGSO += $(OPTFLAG) 
    FFLAGSG += -g
    FFLAGSSO += $(OPTFLAG)   
  endif
endif


# LAPACK_LIBS =  -Wl,-rpath,$(PGI)/linux86-64/7.0-4/lib -L$(PGI)/linux86-64/7.0-4/lib -llapack -lblas
LAPACK_LIBS =  -L$(LAPACK) -llapack -lblas
# LIBS = $(OVERTURE_LIBRARIES) $(APP_LIBRARIES) $(OV_HDF_LIBRARIES) $(OV_COMPILER_LIBS) 
# LIBS +=  $(OV_FORTRAN_LIBRARIES) $(OV_PERL_LIBRARIES) $(OV_OPENGL_LIBRARIES) $(OV_MOTIF_LIBRARIES) $(OV_X_LIBRARIES) $(LAPACK_LIBS) 

# Ubuntu -- move mesa libs forward
LIBS = $(OVERTURE_LIBRARIES) $(OV_OPENGL_LIBRARIES) $(APP_LIBRARIES) $(OV_HDF_LIBRARIES) $(OV_COMPILER_LIBS) 
# LIBS +=  $(PETSC_LIBS) $(OV_FORTRAN_LIBRARIES) $(OV_PERL_LIBRARIES) $(OV_MOTIF_LIBRARIES) $(OV_X_LIBRARIES) $(LAPACK_LIBS) 
LIBS += $(OV_FORTRAN_LIBRARIES) $(OV_PERL_LIBRARIES) $(OV_MOTIF_LIBRARIES) $(OV_X_LIBRARIES) $(LAPACK_LIBS) 


.SUFFIXES:
.SUFFIXES:.C .o .f .o .F .o .bf .f .c .o .cc .o
.C.o:; $(CXX) $(CCFLAGS) -c $<
.cc.o:; $(CXX) $(CCFLAGS) -c $<
.c.o:; $(CC) $(CFLAGS) -c $<
.f.o:; $(FC) $(FFLAGSO) -c $<
.F.o:; $(FC) $(FFLAGSO) -c $<
.bf.f:; bpp $*.bf
.bC.C:; bpp $*.bC
# .bf.o: $*.f ; 
.bC.o: $*.C ; 

# compile some fortran files optimized by default:
$(OBJO) : obj/%.o : %.f
	$(FC) $(FFLAGSO) -o $@ -c $<


# %.o : %.C ; $(CXX) $(CCFLAGS) -c $*.C

# .C: $(LIB_DEPENDENCIES)
#	 $(CXX) $(CCFLAGS) -o $@ $< $(CLIBS) $(FLIBS)  $(GLIBS)

BPP = bpp
%.C : %.bC
	$(BPP) -quiet -clean  $<
%.f : %.bf
	$(BPP) -quiet -clean $<


current = .
mapping = $(current)/../mapping
ogshow = $(current)/../ogshow

VPATH = src:obj

all = cgWave cgwh
# all = cgWave
all: $(all);


Oges = $(Overture)/Oges
linkFiles:
	ln -sf $(Oges)/PETScEquationSolver.C src/
	ln -sf $(Oges)/PETScSolver.C src/
	ln -sf $(Oges)/buildEquationSolvers.C src/

# This is needed with PETSc
buildEquationSolvers.o : $(Oges)/buildEquationSolvers.C; $(CXX) $(CCFLAGS) -DOVERTURE_USE_PETSC -c $(Oges)/buildEquationSolvers.C
PETScEquationSolver.o : $(Oges)/PETScEquationSolver.C; $(CXX) $(CCFLAGS) -DOVERTURE_USE_PETSC -c $(Oges)/PETScEquationSolver.C

OBJC = obj/CgWave.o obj/advance.o obj/plot.o obj/applyBoundaryConditions.o obj/userDefinedKnownSolution.o \
       obj/outputHeader.o obj/printStatistics.o obj/userDefinedForcing.o  obj/updateTimeIntegral.o \
       obj/getTimeStep.o obj/getHelmholtzForcing.o obj/implicit.o obj/getInitialConditions.o obj/saveShow.o obj/getErrors.o \
       obj/takeFirstStep.o obj/deflation.o obj/eigenModes.o \
       obj/rjbesl.o obj/rybesl.o
       
# LCBC files: 
OBJC += obj/initializeLCBC.o obj/LCBC.o obj/LCBC1.o obj/LCBC2.o obj/LCBC_data.o \
        obj/LagrangeDerivFunctions.o obj/LCBC_updateGhost.o obj/LCBC_deriv.o \
        obj/utility.o obj/LCBC_TzFnPointers.o obj/LCBC_annulusMap.o 

# OBJC += obj/initializeLCBC.o obj/LCBC.o obj/LCBC1.o obj/LCBC2.o obj/LCBC_corner.o obj/LCBC_vertex.o obj/LCBC_data.o \
#         obj/LCBC_newFun.o obj/LagrangeDerivFunctions.o \
#         obj/numericalDeriv.o obj/utility.o obj/LCBC_TzFnPointers.o obj/LCBC_annulusMap.o

# Fortran 90 (FN) object files: 
FNOBJO = obj/advWave.o\
         obj/advWave2dOrder2r.o obj/advWave3dOrder2r.o obj/advWave2dOrder2c.o obj/advWave3dOrder2c.o \
         obj/advWave2dOrder4r.o obj/advWave3dOrder4r.o obj/advWave2dOrder4c.o obj/advWave3dOrder4c.o \
         obj/advWave2dOrder6r.o obj/advWave3dOrder6r.o obj/advWave2dOrder6c.o obj/advWave3dOrder6c.o \
         obj/advWave2dOrder8r.o obj/advWave3dOrder8r.o obj/advWave2dOrder8c.o obj/advWave3dOrder8c.o \
         obj/bcOptWave.o \
         obj/bcOptWave2dOrder2.o obj/bcOptWave2dOrder4.o \
         obj/bcOptWave3dOrder2.o \
         obj/getWaveDerivatives.o obj/hierDeriv.o

# New high-order accurate modified equation versions
FNOBJO += obj/advWaveME.o \
          obj/advWaveME2dOrder2r.o obj/advWaveME2dOrder2c.o  obj/advWaveME3dOrder2r.o obj/advWaveME3dOrder2c.o   \
          obj/advWaveME2dOrder4r.o obj/advWaveME2dOrder4c.o  obj/advWaveME3dOrder4r.o obj/advWaveME3dOrder4c.o   \
          obj/advWaveME2dOrder6r.o obj/advWaveME2dOrder6c.o  obj/advWaveME3dOrder6r.o obj/advWaveME3dOrder6c.o   \
          obj/advWaveME2dOrder8r.o obj/advWaveME2dOrder8c.o  obj/advWaveME3dOrder8r.o obj/advWaveME3dOrder8c.o     

# Stencil versions
FNOBJO += obj/advWaveStencil.o \
          obj/advWaveStencil2dOrder2r.o  obj/advWaveStencil2dOrder2c.o \
          obj/advWaveStencil2dOrder4r.o  obj/advWaveStencil2dOrder4c.o \
          obj/advWaveStencil2dOrder6r.o  obj/advWaveStencil2dOrder6c.o \
          obj/advWaveStencil3dOrder2r.o  obj/advWaveStencil3dOrder2c.o \
          obj/advWaveStencil3dOrder4r.o  obj/advWaveStencil3dOrder4c.o \
          obj/advWaveStencil3dOrder6r.o  obj/advWaveStencil3dOrder6c.o \
          obj/advWaveStencil2dOrder8r.o  obj/advWaveStencil2dOrder8c.o \
          obj/advWaveStencil3dOrder8r.o 
          
#          obj/advWaveStencil2dOrder8r.o  obj/advWaveStencil2dOrder8c.o 

# FNOBJO += obj/advWaveStencil.o \
#           obj/advWaveStencil2dOrder2r.o obj/advWaveStencil2dOrder2c.o  obj/advWaveStencil3dOrder2r.o obj/advWaveStencil3dOrder2c.o   \
#           obj/advWaveStencil2dOrder4r.o obj/advWaveStencil2dOrder4c.o  obj/advWaveStencil3dOrder4r.o obj/advWaveStencil3dOrder4c.o   \
#           obj/advWaveStencil2dOrder6r.o obj/advWaveStencil2dOrder6c.o  obj/advWaveStencil3dOrder6r.o obj/advWaveStencil3dOrder6c.o   \
#           obj/advWaveStencil2dOrder8r.o obj/advWaveStencil2dOrder8c.o  obj/advWaveStencil3dOrder8r.o obj/advWaveStencil3dOrder8c.o     


# OBJO = obj/advWave.o\
#         obj/advWave2dOrder2r.o obj/advWave2dOrder2c.o obj/advWave2dOrder4r.o obj/advWave2dOrder4c.o 


# For regression tests: (Note: the master version of checkCheckFiles is in cg/common/src --> could move to Overture)
checkCheckFiles = obj/checkCheckFiles.o 
checkCheckFiles: $(checkCheckFiles); $(CXX) $(CCFLAGS) -o check/checkCheckFiles $(checkCheckFiles) $(LIBS)

# --- CgWave solver ---
# ogmgFiles = ${OvertureCheckout}/ogmg/smooth.o
cgWaveFiles = obj/cgWaveMain.o $(OBJC) $(OBJO) $(FNOBJO) $(petscSolver) $(OGES_PETSC)
cgWave: $(cgWaveFiles) ; $(CXX) $(CCFLAGS) -o bin/cgWave $(cgWaveFiles) $(SLEPC_LIBS) $(PETSC_LIBS) $(LIBS)


# ----- test of coeff matricies with wide stencils -----
# tcmWideStencilFiles = obj/tcmWideStencil.o $(OBJC) $(OBJO) $(FNOBJO)
tcmWideStencilFiles = obj/tcmWideStencil.o obj/cgesl1234.o $(OBJC) $(OBJO) $(petscSolver) $(FNOBJO) 
tcmWideStencil: $(tcmWideStencilFiles); $(CXX) $(CCFLAGS) -o bin/tcmWideStencil $(tcmWideStencilFiles) $(LIBS)

obj/cgesl1234.o : src/cgesl1234.F; $(FC) $(FFLAGSO) -I.  -DOV_USE_DOUBLE  -o $*.o -c $<
# cgesl1234.o:
#   gfortran -O  -fPIC  -fdefault-real-8 -fdefault-double-8   -I/home/henshw/Overture.g/include -I.   -DOV_USE_DOUBLE -c cgesl1234.F  

obj/tcmWideStencil.o : src/tcmWideStencil.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
# obj/solvePETScNull.o : src/solvePETScNull.C; $(CXX) $(CCFLAGS) -o $*.o -c $<


# ----- test quadrature formulae  -----
testQuadFiles = obj/testQuad.o 
# testQuadFiles = obj/testQuad.o obj/getQuadrartureWeights.o 
testQuad: $(testQuadFiles); $(CXX) $(CCFLAGS) -o bin/testQuad $(testQuadFiles) $(LIBS)

obj/testQuad.o : src/testQuad.C; $(CXX) $(CCFLAGS) -o $*.o -c $<  
obj/getQuadratureWeights.o : src/getQuadratureWeights.C; $(CXX) $(CCFLAGS) -o $*.o -c $<  

# ----- test high derivatives  -----
testHighDerivativesFiles = obj/testHighDerivatives.o obj/getWaveDerivatives.o obj/hierDeriv.o
testHighDerivatives: $(testHighDerivativesFiles); $(CXX) $(CCFLAGS) -o bin/testHighDerivatives $(testHighDerivativesFiles) $(LIBS)

obj/testHighDerivatives.o : src/testHighDerivatives.C; $(CXX) $(CCFLAGS) -o $*.o -c $<  

# ----- simple test before adding LCBC  -----
simpleTestFiles = obj/simpleTest.o
simpleTest: $(simpleTestFiles); $(CXX) $(CCFLAGS) -o bin/simpleTest $(simpleTestFiles) $(LIBS)

obj/simpleTest.o : src/simpleTest4.C; $(CXX) $(CCFLAGS) -o $*.o -c $<

# ----- test reordering routine
testReorderFiles = obj/testReorder.o 
testReorder: $(testReorderFiles); $(CXX) $(CCFLAGS) -o bin/testReorder $(testReorderFiles) $(LIBS)
obj/testReorder.o : src/testReorder.C; $(CXX) $(CCFLAGS) -o $*.o -c $<  

# ----- test reordering routine
testPlotFiles = obj/testPlot.o  $(OBJC) $(OBJO) $(FNOBJO) $(petscSolver) $(OGES_PETSC)
# testPlot: $(testPlotFiles); $(CXX) $(CCFLAGS) -o bin/testPlot $(testPlotFiles) $(LIBS)
testPlot: $(testPlotFiles) ; $(CXX) $(CCFLAGS) -o bin/testPlot $(testPlotFiles) $(SLEPC_LIBS) $(PETSC_LIBS) $(LIBS)
obj/testPlot.o : src/testPlot.C; $(CXX) $(CCFLAGS) -o $*.o -c $<    

# ----- test LCBC class : local compatibility boundary conditions -----
testLCBCFiles = obj/testLCBC.o obj/LCBC.o obj/LCBC1.o obj/LCBC2.o obj/LCBC_corner.o obj/LCBC_vertex.o obj/numericalDeriv.o obj/utility.o obj/LCBC_TzFnPointers.o obj/LCBC_annulusMap.o

testLCBC: $(testLCBCFiles); $(CXX) $(CCFLAGS) -o bin/testLCBC $(testLCBCFiles) $(LIBS)

obj/testLCBC.o : src/testLCBC.C; $(CXX) $(CCFLAGS) -o $*.o -c $<

# --- LCBC files - compile opt by default  ---
obj/LagrangeDerivFunctions.o : src/LagrangeDerivFunctions.C; $(CXX) $(CCFLAGSO) -o $*.o -c $< 
obj/LCBC.o : src/LCBC.C; $(CXX) $(CCFLAGSO) -o $*.o -c $< 
obj/LCBC1.o : src/LCBC1.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<
obj/LCBC2.o : src/LCBC2.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<
obj/LCBC_updateGhost.o : src/LCBC_updateGhost.C; $(CXX) $(CCFLAGSO) -o $*.o -c $< 
obj/LCBC_deriv.o : src/LCBC_deriv.C; $(CXX) $(CCFLAGSO) -o $*.o -c $< 
obj/LCBC_data.o : src/LCBC_data.C; $(CXX) $(CCFLAGSO) -o $*.o -c $< 
# obj/LCBC_newFun.o : src/LCBC_newFun.C; $(CXX) $(CCFLAGSO) -o $*.o -c $< 
# obj/numericalDeriv.o : src/numericalDeriv.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<
obj/utility.o : src/utility.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<
obj/LCBC_TzFnPointers.o : src/LCBC_TzFnPointers.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<
obj/LCBC_annulusMap.o : src/LCBC_annulusMap.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<

obj/initializeLCBC.o : src/initializeLCBC.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<
obj/deflation.o : src/deflation.C; $(CXX) $(CCFLAGSO) -o $*.o -c $<
obj/eigenModes.o : src/eigenModes.C; $(CXX) $(CCFLAGS) -o $*.o -c $<

# # ----- test LCBC class : local compatibility boundary conditions -----
# testLCBCFiles = obj/testLCBC.o obj/LCBC.o obj/LCBC1.o obj/LCBC2.o obj/LCBC_corner.o obj/LCBC_vertex.o obj/numericalDeriv.o obj/utility.o
# testLCBC: $(testLCBCFiles); $(CXX) $(CCFLAGS) -o bin/testLCBC $(testLCBCFiles) $(LIBS)

# obj/testLCBC.o : src/testLCBC.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 

# # --- LCBC files ---
# obj/LCBC.o : src/LCBC.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 
# obj/LCBC1.o : src/LCBC1.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 
# obj/LCBC2.o : src/LCBC2.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 
# obj/LCBC_corner.o : src/LCBC_corner.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 
# obj/LCBC_vertex.o : src/LCBC_vertex.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 
# obj/numericalDeriv.o : src/numericalDeriv.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 
# obj/utility.o : src/utility.C; $(CXX) $(CCFLAGS) -o $*.o -c $< 

# --------- CgWaveHoltz ----------

# test:
#   -@echo "usePETSc=[$(usePETSc)]"

#   -@echo "usePETSc=[$(usePETSc)]"
#   -@echo "petscSolver=$(petscSolver)"

test2:;  -@echo "CCFLAGS=$(CCFLAGS)"

test3:; -@echo "SLEPC_INCLUDE = $(SLEPC_INCLUDE)"

# cgwh = obj/cgwh.o obj/CgWaveHoltz.o $(petscSolver) $(OBJC) $(OBJO)
cgwh = obj/cgwh.o obj/CgWaveHoltz.o obj/solveHelmholtz.o obj/solveHelmholtzDirect.o obj/solveEigen.o $(petscSolver) $(OGES_PETSC) $(OBJC) $(OBJO) $(FNOBJO)
cgwh: $(cgwh) 
	$(CXX) $(CCFLAGS) -o bin/cgwh $(cgwh) $(PETSC_LIBS) $(SLEPC_LIBS) $(LIBS)


# bpp files: 
src/cgwh.C: src/cgwh.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include cgwh.bC
src/advance.C: src/advance.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include advance.bC
src/implicit.C: src/implicit.bC boundaryConditionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include implicit.bC
src/applyBoundaryConditions.C: src/applyBoundaryConditions.bC boundaryConditionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include applyBoundaryConditions.bC
src/plot.C: src/plot.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include plot.bC
src/userDefinedKnown.C: src/userDefinedKnown.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include userDefinedKnown.bC
src/userDefinedForcing.C: src/userDefinedForcing.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include userDefinedForcing.bC
src/getHelmholtzForcing.C: src/getHelmholtzForcing.bC boundaryConditionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include getHelmholtzForcing.bC
src/solveHelmholtzDirect.C: src/solveHelmholtzDirect.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include solveHelmholtzDirect.bC
src/solveHelmholtz.C: src/solveHelmholtz.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include solveHelmholtz.bC
src/solveEigen.C: src/solveEigen.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include solveEigen.bC
src/solveSLEPc.C: src/solveSLEPc.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include solveSLEPc.bC

src/tcmWideStencil.C: src/tcmWideStencil.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include tcmWideStencil.bC

src/takeFirstStep.C: src/takeFirstStep.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include takeFirstStep.bC

src/deflation.C: src/deflation.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include deflation.bC
src/eigenModes.C: src/eigenModes.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include eigenModes.bC

# -- optimized advance routines
src/advWave.f90: src/advWave.bf90 
	      @cd src; $(BPP) -clean -quiet -I$(Overture)/include advWave.bf90

src/advWave2dOrder2r.f90 : src/advWave.f90
src/advWave3dOrder2r.f90 : src/advWave.f90
src/advWave2dOrder2c.f90 : src/advWave.f90
src/advWave3dOrder2c.f90 : src/advWave.f90

src/advWave2dOrder4r.f90 : src/advWave.f90
src/advWave3dOrder4r.f90 : src/advWave.f90
src/advWave2dOrder4c.f90 : src/advWave.f90
src/advWave3dOrder4c.f90 : src/advWave.f90

src/advWave2dOrder6r.f90 : src/advWave.f90
src/advWave3dOrder6r.f90 : src/advWave.f90
src/advWave2dOrder6c.f90 : src/advWave.f90
src/advWave3dOrder6c.f90 : src/advWave.f90

src/advWave2dOrder8r.f90 : src/advWave.f90
src/advWave3dOrder8r.f90 : src/advWave.f90
src/advWave2dOrder8c.f90 : src/advWave.f90
src/advWave3dOrder8c.f90 : src/advWave.f90

src/advWaveME.f90: src/advWaveME.bf90 include/update2dOrder6Curvilinear.h include/update3dOrder6Curvilinear.h include/update2dOrder8Curvilinear.h include/update3dOrder8Curvilinear.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include advWaveME.bf90

src/advWaveME2dOrder2r.f90 : src/advWaveME.f90
src/advWaveME3dOrder2r.f90 : src/advWaveME.f90
src/advWaveME2dOrder2c.f90 : src/advWaveME.f90
src/advWaveME3dOrder2c.f90 : src/advWaveME.f90

src/advWaveME2dOrder4r.f90 : src/advWaveME.f90
src/advWaveME3dOrder4r.f90 : src/advWaveME.f90
src/advWaveME2dOrder4c.f90 : src/advWaveME.f90
src/advWaveME3dOrder4c.f90 : src/advWaveME.f90

src/advWaveME2dOrder6r.f90 : src/advWaveME.f90
src/advWaveME3dOrder6r.f90 : src/advWaveME.f90
src/advWaveME2dOrder6c.f90 : src/advWaveME.f90
src/advWaveME3dOrder6c.f90 : src/advWaveME.f90

src/advWaveME2dOrder8r.f90 : src/advWaveME.f90
src/advWaveME3dOrder8r.f90 : src/advWaveME.f90
src/advWaveME2dOrder8c.f90 : src/advWaveME.f90
src/advWaveME3dOrder8c.f90 : src/advWaveME.f90


src/advWaveStencil.f90: src/advWaveStencil.bf90; @cd src; $(BPP) -clean -quiet -I$(Overture)/include advWaveStencil.bf90

src/advWaveStencil2dOrder2r.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder2r.f90 : src/advWaveStencil.f90
src/advWaveStencil2dOrder2c.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder2c.f90 : src/advWaveStencil.f90

src/advWaveStencil2dOrder4r.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder4r.f90 : src/advWaveStencil.f90
src/advWaveStencil2dOrder4c.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder4c.f90 : src/advWaveStencil.f90

src/advWaveStencil2dOrder6r.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder6r.f90 : src/advWaveStencil.f90
src/advWaveStencil2dOrder6c.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder6c.f90 : src/advWaveStencil.f90

src/advWaveStencil2dOrder8r.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder8r.f90 : src/advWaveStencil.f90
src/advWaveStencil2dOrder8c.f90 : src/advWaveStencil.f90
src/advWaveStencil3dOrder8c.f90 : src/advWaveStencil.f90


# -- optimized BC routine
src/bcOptWave.f90: src/bcOptWave.bf90 src/knownSolutionMacros.h maple/defineGetDerivativesMacros.h
	      @cd src; $(BPP) -clean -quiet -I$(Overture)/include bcOptWave.bf90	


src/bcOptWave2dOrder2.f90 : src/bcOptWave.f90
src/bcOptWave2dOrder4.f90 : src/bcOptWave.f90
src/bcOptWave2dOrder6.f90 : src/bcOptWave.f90
src/bcOptWave2dOrder8.f90 : src/bcOptWave.f90
src/bcOptWave3dOrder2.f90 : src/bcOptWave.f90
src/bcOptWave3dOrder4.f90 : src/bcOptWave.f90
src/bcOptWave3dOrder6.f90 : src/bcOptWave.f90
src/bcOptWave3dOrder8.f90 : src/bcOptWave.f90  

# -- optimized derivatives evaluation
src/getWaveDerivatives.f90: src/getWaveDerivatives.bf90; @cd src; $(BPP) -clean -quiet -I$(Overture)/include getWaveDerivatives.bf90

# -- hierachical derivative evalulation
src/hierDeriv.f90: src/hierDeriv.bf90; @cd src; $(BPP) -clean -quiet -I$(Overture)/include hierDeriv.bf90


# dependencies
obj/getDt.o : src/getDt.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/cgwh.o : src/cgwh.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/CgWaveHoltz.o : src/CgWaveHoltz.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/solvePETSc.o : src/solvePETSc.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/solveSLEPc.o : src/solveSLEPc.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/solvePETScNull.o : src/solvePETScNull.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<	
obj/solveHelmholtz.o : src/solveHelmholtz.C src/CgWaveHoltz.h src/knownSolutionMacros.h; $(CXX) $(CCFLAGS) -o $*.o -c $<  
obj/solveHelmholtzDirect.o : src/solveHelmholtzDirect.C src/CgWaveHoltz.h src/knownSolutionMacros.h; $(CXX) $(CCFLAGS) -o $*.o -c $<  
obj/solveEigen.o : src/solveEigen.C src/CgWaveHoltz.h src/knownSolutionMacros.h; $(CXX) $(CCFLAGS) -o $*.o -c $<  

# --- cgWave ---
obj/cgWaveMain.o : src/cgWaveMain.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/CgWave.o : src/CgWave.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/advance.o : src/advance.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/plot.o : src/plot.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/applyBoundaryConditions.o : src/applyBoundaryConditions.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/implicit.o : src/implicit.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/userDefinedKnownSolution.o : src/userDefinedKnownSolution.C src/CgWave.h src/knownSolutionMacros.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/outputHeader.o : src/outputHeader.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/printStatistics.o : src/printStatistics.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/userDefinedForcing.o : src/userDefinedForcing.C src/CgWave.h src/knownSolutionMacros.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/updateTimeIntegral.o : src/updateTimeIntegral.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/getTimeStep.o : src/getTimeStep.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/getHelmholtzForcing.o : src/getHelmholtzForcing.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/getInitialConditions.o : src/getInitialConditions.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/saveShow.o : src/saveShow.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/getErrors.o : src/getErrors.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/takeFirstStep.o : src/takeFirstStep.C src/CgWave.h; $(CXX) $(CCFLAGS) -o $*.o -c $<

obj/rjbesl.o : src/rjbesl.f; $(FC) $(FFLAGSO) -o $@ -c $<  
obj/rybesl.o : src/rybesl.f; $(FC) $(FFLAGSO) -o $@ -c $<  

# For regression testing:
obj/checkCheckFiles.o : src/checkCheckFiles.C; $(CXX) $(CCFLAGS) -o $*.o -c $<


# obj/bcOptWave.o : src/bcOptWave.f90; $(FC) $(FFLAGSO) -ffree-line-length-none -o $*.o -c $<	

# compile f90 files optimized by default:
$(FNOBJO) : obj/%.o : %.f90
	$(FC) $(FFLAGSO) -ffree-line-length-none -finit-real=snan -o $@ -c $<	



clean:  
	rm -f *.o obj/*.o bin/cgwh bin/cgWave


.PRECIOUS: 


