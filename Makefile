#
# CgWave: Composite Grid Wave Equation Solver 
# CgWaveHoltz: Composite Grid Helmholtz Solver using the WaveHoltz algorithm of Daniel Appelo 
#
# NOTE: To compile optimized:
#   setenv COMPILE [opt|dbg]
#


include ${Overture}/make.options

usePETSc := on
# usePETSc := off
ifeq ($(usePETSc),on)
  usePETSc   = on
  petscSolver = obj/solvePETSc.o

  OGES_PETSC = buildEquationSolvers.o PETScEquationSolver.o
  PETSC_INCLUDE = -I$(PETSC_DIR)/include  -I$(PETSC_DIR)/$(PETSC_ARCH)/include -DOVERTURE_USE_PETSC -I$(PETSC_LIB)/include -I$(PETSC_DIR)/include/mpiuni
  PETSC_LIBS = -Wl,-rpath,$(PETSC_LIB) -L$(PETSC_LIB) -lpetsc

else
  usePETSc   = off
  petscSolver = obj/solvePETScNull.o

  OGES_PETSC = 
  PETSC_INCLUDE = 
  PETSC_LIBS = 
endif




CCFLAGS = $(OV_CXX_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include -I$(OpenGL)/include $(USE_PPP_FLAG)
CCFLAGS += $(PETSC_INCLUDE)
CFLAGS = $(OV_CC_FLAGS) -I. -I$(Overture)/include -I$(APlusPlus)/include

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
  CCFLAGS += -w -O 
  FFLAGS  += -O
  FFLAGSSO += -O
else
	ifeq ($(COMPILE),dbg)
    # debug:
    CCFLAGS += -w -g
    FFLAGS  += -g
    FFLAGSO += -g
    FFLAGSSO += g   
  else
    # default case:
    CCFLAGS += -w -g 
    FFLAGS  += -g 
    FFLAGSO += -O 
    FFLAGSG += -g
    FFLAGSSO += O   
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
.bf.o: $*.f ; 
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
       obj/takeFirstStep.o

# Fortran 90 (FN) object files: 
FNOBJO = obj/advWave.o\
         obj/advWave2dOrder2r.o obj/advWave3dOrder2r.o obj/advWave2dOrder2c.o obj/advWave3dOrder2c.o \
         obj/advWave2dOrder4r.o obj/advWave3dOrder4r.o obj/advWave2dOrder4c.o obj/advWave3dOrder4c.o \
         obj/advWave2dOrder6r.o obj/advWave3dOrder6r.o obj/advWave2dOrder6c.o obj/advWave3dOrder6c.o \
         obj/advWave2dOrder8r.o obj/advWave3dOrder8r.o obj/advWave2dOrder8c.o obj/advWave3dOrder8c.o \
         obj/bcOptWave.o obj/getWaveDerivatives.o

# OBJO = obj/advWave.o\
#         obj/advWave2dOrder2r.o obj/advWave2dOrder2c.o obj/advWave2dOrder4r.o obj/advWave2dOrder4c.o 


# For regression tests: (Note: the master version of checkCheckFiles is in cg/common/src --> could move to Overture)
checkCheckFiles = obj/checkCheckFiles.o 
checkCheckFiles: $(checkCheckFiles); $(CXX) $(CCFLAGS) -o check/checkCheckFiles $(checkCheckFiles) $(LIBS)

# --- CgWave solver ---
cgWaveFiles = obj/cgWaveMain.o $(OBJC) $(OBJO) $(FNOBJO)
cgWave: $(cgWaveFiles) ; $(CXX) $(CCFLAGS) -o bin/cgWave $(cgWaveFiles) $(LIBS)


# ----- test of coeff matricies with wide stencils -----
# tcmWideStencilFiles = obj/tcmWideStencil.o $(OBJC) $(OBJO) $(FNOBJO)
tcmWideStencilFiles = obj/tcmWideStencil.o obj/cgesl1234.o $(OBJC) $(OBJO) $(petscSolver) $(FNOBJO) 
tcmWideStencil: $(tcmWideStencilFiles); $(CXX) $(CCFLAGS) -o bin/tcmWideStencil $(tcmWideStencilFiles) $(LIBS)

obj/cgesl1234.o : src/cgesl1234.F; $(FC) $(FFLAGSO) -I.  -DOV_USE_DOUBLE  -o $*.o -c $<
# cgesl1234.o:
#   gfortran -O  -fPIC  -fdefault-real-8 -fdefault-double-8   -I/home/henshw/Overture.g/include -I.   -DOV_USE_DOUBLE -c cgesl1234.F  

obj/tcmWideStencil.o : src/tcmWideStencil.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/solvePETScNull.o : src/solvePETScNull.C; $(CXX) $(CCFLAGS) -o $*.o -c $<


# ----- test quadrature formulae  -----
testQuadFiles = obj/testQuad.o 
# testQuadFiles = obj/testQuad.o obj/getQuadrartureWeights.o 
testQuad: $(testQuadFiles); $(CXX) $(CCFLAGS) -o bin/testQuad $(testQuadFiles) $(LIBS)

obj/testQuad.o : src/testQuad.C; $(CXX) $(CCFLAGS) -o $*.o -c $<  
obj/getQuadratureWeights.o : src/getQuadratureWeights.C; $(CXX) $(CCFLAGS) -o $*.o -c $<  

# --------- CgWaveHoltz ----------

test: 
	-@echo "usePETSc=[$(usePETSc)]"
	-@echo "petscSolver=$(petscSolver)"

# cgwh = obj/cgwh.o obj/CgWaveHoltz.o $(petscSolver) $(OBJC) $(OBJO)
cgwh = obj/cgwh.o obj/CgWaveHoltz.o obj/solveHelmholtz.o $(petscSolver) $(OGES_PETSC) $(OBJC) $(OBJO) $(FNOBJO)
cgwh: $(cgwh) 
	$(CXX) $(CCFLAGS) -o bin/cgwh $(cgwh) $(PETSC_LIBS) $(LIBS)


# bpp files: 
src/advance.C: src/advance.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include advance.bC
src/implicit.C: src/implicit.bC boundaryConditionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include implicit.bC
src/applyBoundaryConditions.C: src/applyBoundaryConditions.bC boundaryConditionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include applyBoundaryConditions.bC
src/plot.C: src/plot.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include plot.bC
src/userDefinedKnown.C: src/userDefinedKnown.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include userDefinedKnown.bC
src/userDefinedForcing.C: src/userDefinedForcing.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include userDefinedForcing.bC
src/getHelmholtzForcing.C: src/getHelmholtzForcing.bC boundaryConditionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include getHelmholtzForcing.bC
src/solveHelmholtz.C: src/solveHelmholtz.bC src/knownSolutionMacros.h; @cd src; $(BPP) -clean -quiet -I$(Overture)/include solveHelmholtz.bC

src/tcmWideStencil.C: src/tcmWideStencil.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include tcmWideStencil.bC

src/takeFirstStep.C: src/takeFirstStep.bC; @cd src; $(BPP) -clean -quiet -I$(Overture)/include takeFirstStep.bC

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

# -- optimized BC routine
src/bcOptWave.f90: src/bcOptWave.bf90 src/knownSolutionMacros.h
	      @cd src; $(BPP) -clean -quiet -I$(Overture)/include bcOptWave.bf90	

# -- optimized derivatives evaluation
src/getWaveDerivatives.f90: src/getWaveDerivatives.bf90; @cd src; $(BPP) -clean -quiet -I$(Overture)/include getWaveDerivatives.bf90

# dependencies
obj/getDt.o : src/getDt.C; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/cgwh.o : src/cgwh.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/CgWaveHoltz.o : src/CgWaveHoltz.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/solvePETSc.o : src/solvePETSc.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<
obj/solvePETScNull.o : src/solvePETScNull.C src/CgWaveHoltz.h; $(CXX) $(CCFLAGS) -o $*.o -c $<	
obj/solveHelmholtz.o : src/solveHelmholtz.C src/CgWaveHoltz.h src/knownSolutionMacros.h; $(CXX) $(CCFLAGS) -o $*.o -c $<	

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

# For regression testing:
obj/checkCheckFiles.o : src/checkCheckFiles.C; $(CXX) $(CCFLAGS) -o $*.o -c $<


# obj/bcOptWave.o : src/bcOptWave.f90; $(FC) $(FFLAGSO) -ffree-line-length-none -o $*.o -c $<	

# compile f90 files optimized by default:
$(FNOBJO) : obj/%.o : %.f90
	$(FC) $(FFLAGSO) -ffree-line-length-none -o $@ -c $<	



clean:  
	rm -f obj/*.o bin/cgwh bin/cgWave


.PRECIOUS: 


