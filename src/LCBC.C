#include "LCBC.h"
#include "utility.h"
#include "LCBCmacros.h"

void extrapolateAllGhost(double *&unp1, int *mask, int unp1Dim[3], int indexRange[3][2], int *faceEval, int numGhost, int userNumGhost, int p, int dim);

void Default_bn(double bn[3], double *arg);
double Default_Gn(double *arg);
double Default_Fn(double *arg);
double Default_Cn(double *arg);
int countKnownBdryFaces(int vec[3]);
            
/// \brief Initialize the LCBC object defined in LCBC.h
/// \param dimIn (input): the dimensions in space
/// \param orderInSpace (input): the order of accuracy in space
/// \param orderInTimeIn (input): the order of accuracy in time
/// \param numGridPoints (input): a vector holding the number of grid points in each direction
/// \param numGhostIn (input): the number of ghost points
/// \param faceEvalIn (input): a vector holding the type of BCs on each face
/// faceEvalIn[face] = -1: Periodic BC
/// faceEvalIn[face] =  0: No BC
/// faceEvalIn[face] =  1: Dirichlet BC
/// faceEvalIn[face] =  2: Neumann BC
/// face = side + 2*axis
/// \param coefIn (input): an LcbcData object holding PDE coefficients information (see LCBC_data.C/h for more information about the class)
/// \param CnIn (input): a function handle for the PDE coefficients
/// \param GnIn (input): a function handle for the BCs on each face
/// \param FnIn (input): a function handle for the forcing
/// Each of the function handles above is defined using F1 type (typedef double (*F1)(double *), defined in LCBC.h) where the double argument, call it 'arg', is such that arg[0]: face (boundary) or coefficient number (PDE coefficient) or empty (forcing), arg[1]: x, arg[2]: y, arg[3]: z or t (depending on the dimension), arg[4]: t (if dimension = 3)
/// \param bnIn (input): a function handle for the mapped Neumann boundary operator
/// The function bnIn is defined using F2 type (typedef void (*F2)(double[3], double *), defined in LCBC.h) where the first argument, call it bn[3], holds the coefficients of u_x, u_y, u_z respectively and the second argument, call it arg, is such that arg[0]: face, arg[1]: x, arg[2]: y, arg[3]: z or t (depending on the dimension) and arg[4]: t (if dimension = 3)
/// \param cstCoefIn (input): a boolean to tell if the problem has constant coefficients
/// \param zeroBCIn (input): a boolean to tell if the problem has zero BC data on each face
/// \param noForcingIn (input): a boolean to tell if there is no forcing
void Lcbc::initialize(int dimIn, int orderInSpace, int orderInTimeIn, int *numGridPoints, int numGhostIn, int *faceEvalIn, LcbcData *coefIn, F1 CnIn, F1 GnIn, F1 FnIn, F2 bnIn, bool cstCoefIn, bool *zeroBCIn, bool noForcingIn){
    
    isInitialized = true;
    
    dim = dimIn;
    p = orderInSpace/2;
    extraDataGhost = (p-1); // extra ghost points needed for data grid functions
    faceNum = 2*dim;
    int n = (2*p+1);
    totalVarNum = n*n*dimBasedValue(dim, 1, n);
    auxiliaryEqNum = 0;
    getDerivCoef(p);
    numGhost = p; // number of ghost points needed
    orderInTime = orderInTimeIn;
    getLagrangeData(); // get the data needed from the Lagrange polynomial bases
    G.prepareGrid(numGridPoints,numGhost,dim); // prepare the spatial grid
    updateNeededParameters(); // update bdryTypeParam objects defined in LCBC.h
    
    /* ensure the number of ghost points used by the user is valid */
    if(numGhostIn<numGhost){
        printf("ERROR: need at least %d ghost points to use LCBC\n",numGhost);
        exit(-1);
    }
    else{
        userNumGhost = numGhostIn;
    }
    
    /* Set faceEval */
    faceEval = faceEvalIn;
    if(faceEval == NULL){
        faceEvalNewed = true;
        faceEval = new int[(2*dim)];
        for(int face = 0; face<(2*dim); face++){
            faceEval[face] = 1; // Set the faces to Dirichlet BCs as default option
        }// end of face loop
    }// end of faceEval condition
    
    /* set the forcing function */
    if(FnIn == NULL){
        Fn = Default_Fn; // default forcing to zero if not given
        defaultedFn = true;
    }else{
        Fn = FnIn;
    }
    noForcing = noForcingIn;
    
    /* set the BC functions */
    if(GnIn == NULL){
        Gn = Default_Gn; // default the BCs to zero if not given
        defaultedGn = true;
    }else{
        Gn = GnIn;
    }
    if(bnIn == NULL){
        bn = Default_bn; // default for bn if not given by user (see Default_bn)
    }else{
        bn = bnIn;
    }
    if(zeroBCIn == NULL){
        zeroBC = new bool[(2*dim)];
        zeroBCnewed = true;
        for(int face = 0; face<(2*dim); face++){
            zeroBC[face] = false;
        }// end of face
    }else{
        zeroBC = zeroBCIn;
    }
    
    /* set the PDE coefficient functions */
    maxCoefNum = ((dim+1)*(dim+2)/2); // the maximum number of coefficients expected based on the dimension
    coef = new LcbcData[(2*dim)]; // objects to carry the PDE coefficient data
    if(coefIn == NULL){
        if(CnIn == NULL){
            Cn = Default_Cn; // default to the Laplacian operator if the coefficients are not given
            cstCoef = true;
        }else{
            Cn = CnIn;
            cstCoef = cstCoefIn;
        }
        getCoef(); // put the coefficients in the corresponding LcbcData objects
    }else{
        /* If user adds coefficient information, we check the input to make sure it is entered correctly.
           We extrapolate any extra values if needed. */
        cstCoef = cstCoefIn;
        bool fix = true; // ensure we extrapolate coefficient values on extra ghost points for data if none are given from the user
        bool overRideCheck = (cstCoef)?(true):(false);
        checkUserData(coef, coefIn, G.Nx, G.Ngx, faceEval, dim, p, extraDataGhost, coefficient, fix, overRideCheck);
    }// end of else
    
    /* prepare the LCBC data functions */
    forcingData = new LcbcData[(2*dim)];
    bdryData = new LcbcData[(2*dim)];
    
    /* prepare memory for the 3D 8th-order accurate scheme here */
    if((dim == 3) && (p==4)){
        memoryComponents = 33 + 25;
        int componentSize = (4*p + 1)*(4*p + 1)*(2*p + 1);
        preallocateMemory(memoryComponents, componentSize);
        preallocatedMemory = true;
    }
    
}// end of initialize

/// \brief Analyze the input data from the user and
/// \param gn (input): LcbcData object containing bdry information from the user
/// \param fn (input): LcbcData object containing forcing information from the user
void Lcbc::analyzeUserInput(LcbcData *gn, LcbcData *fn){
    if(!isInitialized){
        printf("ERROR: must initialize the Lcbc object first\n");
        exit(-1);
    }
    
    if(fn == NULL){
        if(defaultedFn){
            noForcing = true;
        }
    }else if(!noForcing){
        checkUserData(forcingData, fn, G.Nx, G.Ngx, faceEval, dim, p, extraDataGhost, forcing);
        initializedForcingData = true;
    }
    if(gn == NULL){
        if(defaultedGn){
            if(!zeroBCnewed){
                zeroBC = new bool[(2*dim)];
                zeroBCnewed = true;
            }
            for(int face = 0; face<(2*dim); face++)
                zeroBC[face] = true;
        }
    }else{
        checkUserData(bdryData, gn, G.Nx, G.Ngx, faceEval, dim, p, extraDataGhost, boundary);
        initializedBdryData = true;
    }
    analyzedUserData = true;
}// end of analyze user data

/// \brief update the values of the grid solution at the ghost points
/// \param unp1 (input/output): the grid solution at time t
/// \param mask (input): the mask given from cgWave
/// \param t (input): the time t
/// \param dt (input): the time-step
/// \param gn (input): an LcbcData object carrying boundary information from the user (see LCBC_data.h)
/// \param fn (input): an LcbcData object carrying forcing information from the user (see LCBC_data.h)
void Lcbc::updateGhost(double *&unp1, int *mask, double t, double dt, LcbcData *&gn, LcbcData *&fn){
    
#if ( EXTRAPOLATE_ALL_GHOST==1 )
    /* An option to extrapolate all ghost values if needed for debugging */
    extrapolateAllGhost(unp1, mask, G.Ngx, G.indexRange, faceEval, numGhost, userNumGhost, p, dim);
#endif
    
    /* analyze the data coming in from the user and fix it if needed */
    if(!analyzedUserData)
        analyzeUserInput(gn, fn);
    /* ================================================== */

    int NU = p+1; // the maximum number of primary CBCs
    double **R = new double*[(NU*faceNum)]; // a vector that will hold boundary and forcing values from the right-hand-side of each BC and primary CBC on each boundary point. This will be used to form the vector R(t) in the lcbc.pdf in the stencil approach (page 20).
    
    /* prepare the derivatives of data functions over the stencil and update face ghost */
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                
#if PRINT_PROCESS == 1
                    printf("updating face %d\n", face);
#endif
                /* Fill or fix forcing and boundary information from user in LcbcData objects */
                if(initializedForcingData){
                    if(forcingData[face].fill){
                        fillData(forcingData[face], fn[face]);
                    }else{
                        forcingData[face].Fn = fn[face].Fn;
                    }
                    
                    if(forcingData[face].fix){
                        fixData(forcingData[face], p, dim, axis);
                    }
                }
                if(initializedBdryData && (!zeroBC[face])){
                    bdryData[face].Fn = gn[face].Fn;
                }
                
                if(!(faceParam[face].exists)){
                    faceParam[face].initialize(axis, side, mask, numGhost, userNumGhost, G.indexRange, G.Ngx, p, dim, faceEval, addAuxEqns);
                }
                
                /* Get the vector R defined and described above */
                getDataDeriv(R, t, dt, axis, side, bdryData[face], forcingData[face]);
                
                /* Update the solution at faces (away from edges or corners) */
                updateBdryGhost(unp1, mask, R, oneAxis, FaceMat[face], t, dt, &axis, &side);
                
            }// end if faceEval
        }// end side loop
    }// end axis loop
    
#if SKIP_LCBC_EDGE==0
    //Update the 2D corners or 3D edges
    for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<3; varAxis++){
        
        int fixedAxis[2]; getOtherAxes(fixedAxis, varAxis);
        
        for(int side1 = 0; side1<2; side1++){
            for(int side2 = 0; side2<2; side2++){
                int edgeFace1 = side1 + 2*fixedAxis[0];
                int edgeFace2 = side2 + 2*fixedAxis[1];
                
                int fixedSide[] = {side1,side2};
                int corner = ind3(side1,side2,dimBasedValue(dim, 0, varAxis),2,2,3);
                
                /* Case of an edge along non-periodic faces */
                if((faceEval[edgeFace1]>0) && (faceEval[edgeFace2]>0)){
#if PRINT_PROCESS == 1
                    printf("updating edge (%d,%d)\n", edgeFace1, edgeFace2);
#endif
                    updateBdryGhost(unp1, mask, R, twoAxis, EdgeMat[corner], t, dt, fixedAxis, fixedSide);
                }
                /* Case of an edge with one interpolation boundary */
                else if(((faceEval[edgeFace1]==0) && (faceEval[edgeFace2]>0))||((faceEval[edgeFace1]>0) && (faceEval[edgeFace2]==0))){
#if SKIP_LCBC_NEAR_CORNER==0
#if PRINT_PROCESS == 1
                    printf("updating near edge (%d,%d)\n", edgeFace1, edgeFace2);
#endif
                    updateGhostNearZeroBdry(unp1, mask, R, nearEdge, EdgeMat[corner], t, dt, fixedAxis, fixedSide);
#endif

                }
            }// end of side2
        }// end of side1
    }// end of varAxis loop
#endif
    
#if SKIP_LCBC_EDGE==0
    //Update the 2D corners or 3D edges that are periodic
    for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<3; varAxis++){
        int fixedAxis[2]; getOtherAxes(fixedAxis, varAxis);
        for(int side1 = 0; side1<2; side1++){
            for(int side2 = 0; side2<2; side2++){
                int edgeFace1 = side1 + 2*fixedAxis[0];
                int edgeFace2 = side2 + 2*fixedAxis[1];
                if((faceEval[edgeFace1]<0 && faceEval[edgeFace2]>0)||(faceEval[edgeFace1]>0 && faceEval[edgeFace2]<0)){
#if PRINT_PROCESS == 1
                    printf("updating periodic edge (%d,%d)\n", edgeFace1, edgeFace2);
#endif
                    updateEdgeGhostPeriodic(unp1, mask, side1, side2, varAxis, fixedAxis);
                }
            }// end of side2
        }// end of side1
    }// end of varAxis loop
#endif
    
#if SKIP_LCBC_CORNER==0
    //     Update the 3D vertices
    if(dim == 3){
        int fixedAxis[] = {0,1,2};
        for(int side2 = 0; side2<2; side2++){
            for(int side1 = 0; side1<2; side1++){
                for(int side0 = 0; side0<2; side0++){
                    
                    int fixedSide[] = {side0,side1,side2};
                    
                    int vertex = ind3(side0,side1,side2,2,2,2);
                    
                    int cornerFace0 = side0;
                    int cornerFace1 = side1 + 2;
                    int cornerFace2 = side2 + 4;
                    /* Case of a vertex formed by three physical boundaries */
                    if((faceEval[cornerFace0]>0) && (faceEval[cornerFace1]>0) && (faceEval[cornerFace2]>0)){
#if PRINT_PROCESS == 1
                        printf("updating vertex (%d,%d,%d)\n", cornerFace0, cornerFace1,cornerFace2);
#endif
                        updateBdryGhost(unp1, mask, R, threeAxis, VertexMat[vertex], t, dt, fixedAxis, fixedSide);
                    }
                    else{
                        /* Case of a vertex where at least one boundary is an interpolation boundary */
                        int vec[] = {faceEval[cornerFace0],faceEval[cornerFace1],faceEval[cornerFace2]};
                        int knownFaces = countKnownBdryFaces(vec);
                        
                        if (knownFaces>0 && !((faceEval[cornerFace0]<0) || (faceEval[cornerFace1]<0) || (faceEval[cornerFace2]<0))&& ((faceEval[cornerFace0]==0) || (faceEval[cornerFace1]==0) || (faceEval[cornerFace2]==0))){
#if PRINT_PROCESS == 1
                            printf("updating near vertex (%d,%d,%d)\n", cornerFace0, cornerFace1,cornerFace2);
#endif
                            int vec[] = {faceEval[cornerFace0],faceEval[cornerFace1],faceEval[cornerFace2]};
                            int knownFaces = countKnownBdryFaces(vec);
                            
                            updateGhostNearZeroBdry(unp1, mask, R, nearVertex[(knownFaces-1)], VertexMat[vertex], t, dt, fixedAxis, fixedSide);
                        }
                    }
                }// end of side0 loop
            }// end of side1 loop
        }// end of side2 loop
    }// end of if statement
#endif
    
#if SKIP_LCBC_CORNER==0
    //     Update the 3D vertices that are periodic
    if(dim == 3){
        for(int side2 = 0; side2<2; side2++){
            for(int side1 = 0; side1<2; side1++){
                for(int side0 = 0; side0<2; side0++){
                    
                    int cornerFace0 = side0;
                    int cornerFace1 = side1 + 2;
                    int cornerFace2 = side2 + 4;
                    
                    int vec[] = {faceEval[cornerFace0],faceEval[cornerFace1],faceEval[cornerFace2]};
                    int knownFaces = countKnownBdryFaces(vec);
                    
                    if ((knownFaces>0)&&((faceEval[cornerFace0]<0) || (faceEval[cornerFace1]<0) || (faceEval[cornerFace2]<0)) && !((faceEval[cornerFace0]<0) && (faceEval[cornerFace1]<0) && (faceEval[cornerFace2]<0))){
                        
#if PRINT_PROCESS == 1
                            printf("updating periodic vertex (%d,%d,%d)\n", cornerFace0, cornerFace1,cornerFace2);
#endif
                        
                        updateVertexGhostPeriodic(unp1, mask, side0, side1, side2);
                    }
                }// end of side0 loop
            }// end of side1 loop
        }// end of side2 loop
    }// end of if statement
#endif
    
    /* free the R variable */
    for(int face = 0; face<faceNum; face++){
        for(int nu = 0; nu<NU; nu++){
            
            if(faceEval[face]>0 && R[(nu+face*NU)]!=NULL){
                delete [] R[(nu+face*NU)];
                R[(nu+face*NU)] = NULL;
            }
        }// end of nu loop
    }// end of face loop
    delete [] R; R = NULL;
    
}// end of updateGhost

/// \brief prepare data coming in from Lagrange polynomials evaluated on integer values needed in the LCBC procedure
void Lcbc::getLagrangeData(){
    /* This function evaluate L_{i}(z) where z in [-w,w] and i in [-p,p] */
    int n = 2*p+1;
    int w = 2*p+1;
    int m = 2*w+1;
    LagrangeData = new double[(n*m)];
    
    int ci = 0, cz = 0;
    for (int i=-p; i<=p; i++) {
        for (int z=-w; z<=w; z++) {
            LagrangeData[ind2(ci,cz,n,m)] = LagrangePoly(z,i);
            cz = cz + 1;
        }
        cz = 0;
        ci = ci + 1;
    }
    LagrangeData_center = w;
}

/// \brief evaluate the Lagrange polynomial with index 'index' at 'z'
/// \param z (input): argument of the Lagrange polynomial
/// \param index (input): index of the Lagrange polynomial
double Lcbc::LagrangePoly(int z,int index){
    /* This function evaluates L_{index}(z) where L is the Lagrange polynomial */
    
    double product = 1;
    for (int k = -p; k<=p; k++) {
        if (k!=index) {
            product = product*((double)(z-k))/((double)(index-k));
        }
    }
    return product;
}

/// \brief The default coefficient functions of the PDE set to the Laplacian operator
/// \param arg (input): arg[0]: coefficient number (see numbering within function), arg[1]: x, arg[2]: y, arg[3]: z or t (depending on dimension), arg[4]:t (if dimension = 3)
double Default_Cn(double *arg){
    /* These are the default coefficient functions */
    /* 2D arg = {Coef#,x,y,z}, Coef# = {0,1,2,3,4,5} for {c11,c22,c1,c2,c12,c0} */
    /* 3D arg = {Coef#,x,y,z}, Coef# = {0,1,2,3,4,5,...} for {c11,c22,c33,c1,c2,c3,c12,c12,c23,c0} */
    
    double C = 0;
    if(arg[0] == 0){
        C = 1;
    }else if(arg[0] == 1){
        C = 1;
    }else{
        C = 0;
    }
    return C;
}

/// \brief the default boundary data set to zero
/// \param arg (input): arg[0]: face number, arg[1]: x, arg[2]: y, arg[3]: z or t (depending on dimension), arg[4]:t (if dimension = 3)
double Default_Gn(double *arg){
    /* This is the default boundary conditions functions */
    /* arg = {side,x,y,t}, side = {0,1,2,3} for {left,bottom,right,top} */
    
    double G = 0;
    return G;
}

/// \brief the default mapped Neumann boundary operator
/// \param bn (input): bn[0], bn[1], bn[2]: coefficients of u_x, u_y, and u_z, respectively
/// \param arg (input): arg[0]: face number, arg[1]: x, arg[2]: y, arg[3]: z or t (depending on dimension), arg[4]:t (if dimension = 3)
void Default_bn(double bn[3], double *arg){
    int face = arg[0];
    if(face == 0 || face == 1){ // axis = 0
        bn[0] = 1; bn[1] = 0; bn[2] = 0;
    }
    else if(face == 2 || face == 3){ // axis = 1
        bn[0] = 0; bn[1] = 1; bn[2] = 0;
    }
    else if(face == 4 || face == 5){// axis = 2
        bn[0] = 0; bn[1] = 0; bn[2] = 1;
    }
}// end of default bn

/// \brief the default forcing data set to zero
/// \param arg (input): arg[0]: empty, arg[1]: x, arg[2]: y, arg[3]: z or t (depending on dimension), arg[4]:t (if dimension = 3)
double Default_Fn(double *arg){
    /* This is the default forcing function */
    /* arg = {0,x,y,z,t}                    */
    
    double F = 0;
    return F;
}

/// \brief A function to fill the LcbcData object carrying the PDE coefficients
void Lcbc::getCoef(){
    int maxCoefNum = (((dim+1)*(dim + 2))/2);
    
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                int bdryRange[3][2],lth[3], wth[3]; getCoefGridLth(G.indexRange, lth, wth, bdryRange, axis, side, dim, p);
                
                coef[face].initialize(lth, wth, maxCoefNum);
                int i[3];
                for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                    for(i[2] = 0; i[2]<lth[2]; i[2]++){
                        for(i[1] = 0; i[1]<lth[1]; i[1]++){
                            for(i[0] = 0; i[0]<lth[0]; i[0]++){
                                
                                double arg[] = {((double) coefNum),(G.x[0][bdryRange[0][0]] + ((- wth[0] + i[0])*G.dx[0])),
                                    (G.x[1][bdryRange[1][0]] + ((- wth[1] + i[1])*G.dx[1])),
                                    (G.x[2][bdryRange[2][0]] + ((- wth[2] + i[2])*G.dx[2]))};
                                
                                coef[face].Fn[coefNum][ind(i,lth)] = Cn(arg);
                                
                            }// end of i[0] loop s
                        }// end of i[1] loop
                    }// end of i[2] loop
                    
                }// end of coefNum
            }// end of if faceEval
        }// end of side loop
    }// end of axis loop
    
}// getBdryGridFunctions

/// \brief Get the dimensions of the grid functions carrying the PDE coefficient values on the grid
/// \param indexRange (input): range of indices on each face
/// \param lth (output): the dimensions of the grid function in each axis
/// \param wth (output): the center point of the coefficient grid functions in each axis
/// \param bdryRange (output): the range of indices on a given boundary
/// \param fixedAxis (input): the fixed axis at a given boundary
/// \param fixedSide (input): the fixed side at a given boundary
/// \param dim (input): the dimension
/// \param p (input): orderInSpace/2
void Lcbc::getCoefGridLth(int indexRange[3][2], int lth[3], int wth[3], int bdryRange[3][2], int fixedAxis, int fixedSide, int dim, int p){
    
    if(cstCoef){
        for(int axis = 0; axis<3; axis++){
            wth[axis] = 0;
            lth[axis] = 1;
            for(int side = 0; side<2; side ++){
                bdryRange[axis][side] = 0;
            }// end of side loop
        }// end of axis loop
    }// end of if cstCoef
    else{
        for(int axis = 0; axis<3; axis++){
            wth[axis] = ((axis==fixedAxis)?(p):(p + extraDataGhost));
            if(axis==dim){
                wth[axis] = 0;
            }
            for(int side = 0; side<2; side++){
                if(axis == fixedAxis){
                    bdryRange[axis][side] = indexRange[fixedAxis][fixedSide];
                }
                else{
                    bdryRange[axis][side] = indexRange[axis][side];
                }// end of if axis
            }// end of side
            int bdryRangeLth = bdryRange[axis][1] - bdryRange[axis][0] + 1;
            lth[axis] = (2*wth[axis] + bdryRangeLth);
        }// end of axis
    }// end of else (non cst coef)
}

/// \brief allocates space in 'memory' double pointer to be used in the interior schemes to prevent dynamic allocation at each time-step
/// \param components (input): number of memory components needed
/// \param componentSize (input): the size of each memory component to be preallocated
void Lcbc::preallocateMemory(int components, int componentSize){
    memory = new char*[components];
    for(int i = 0; i<components; i++){
        memory[i] = new char[(sizeof(double)*componentSize)];
    }
}

/// \brief deletes the preallocated memory when time-stepping is finished
/// \param components (input): number of memory components used
void Lcbc::deleteMemory(int components){
    for(int i = 0; i<components; i++){
        delete [] memory[i];
    }// end of i loop
    
    delete [] memory;
}

/// \brief Free general variables used in the LCBC procedure
void Lcbc::freeVariables(){
    
    delete [] LagrangeData; LagrangeData = NULL;
    
    delete [] coef;
    delete [] forcingData;
    delete [] bdryData;
    
    if(faceEvalNewed)
        delete [] faceEval;
    
    if(zeroBCnewed)
        delete [] zeroBC;
    
    if(preallocatedMemory)
        deleteMemory(memoryComponents);
    
}// end of freeVariables

/// \brief Free face variables used in the LCBC procedure
void Lcbc::freeFaceVariables(){
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = ind2(side,axis,2,3);
            int bdryNg = (cstCoef)?(1):(faceParam[face].bdryNg);
            if(FaceMat[face].flag){
                for(int bdryPoint = 0; bdryPoint<bdryNg; bdryPoint++){
                    delete [] FaceMat[face].CaVec[bdryPoint]; FaceMat[face].CaVec[bdryPoint] = NULL;
                    delete [] FaceMat[face].CbVec[bdryPoint]; FaceMat[face].CbVec[bdryPoint] = NULL;
                }
                
                delete [] FaceMat[face].eqNum; FaceMat[face].eqNum = NULL;
                delete [] FaceMat[face].CaVec; FaceMat[face].CaVec = NULL;
                delete [] FaceMat[face].CbVec; FaceMat[face].CbVec = NULL;
            }
        }// end of side
    }// end of axis
}// end of freeFaceVariables

/// \brief Free corner (2D) or edge (3D) variables used in the LCBC procedure
void Lcbc::freeEdgeVariables(){
    for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<=2; varAxis++){
        
        int bdryNg = (cstCoef)?(1):(dimBasedValue(dim, 1, G.Ngx[varAxis]));
        
        for(int side1 = 0; side1<2; side1++){
            for(int side2 = 0; side2<2; side2++){
                
                int corner = ind3(side1,side2,dimBasedValue(dim, 0, varAxis),2,2,3);
                
                if(EdgeMat[corner].flag){
                    
                    for(int bdryPt = 0; bdryPt<bdryNg; bdryPt++){
                        delete [] EdgeMat[corner].CaVec[bdryPt];
                        delete [] EdgeMat[corner].CbVec[bdryPt];
                    }// end of bdryPt loop
                    
                    delete [] EdgeMat[corner].CaVec;
                    delete [] EdgeMat[corner].CbVec;
                    
                    delete [] EdgeMat[corner].eqNum;
                }// end of if flag
            }// end of side2
        }// end of side1
    }// end of varAxis
}// end of freeEdgeVariables

/// \brief Free 3D vertex variables used in the LCBC procedure
void Lcbc::freeVertexVariables(){
    for(int side2 = 0; side2<2; side2++){
        for(int side1 = 0; side1<2; side1++){
            for(int side0 = 0; side0<2; side0++){
                int vertex = ind3(side0,side1,side2,2,2,2);
                
                if(VertexMat[vertex].flag){
                    delete [] VertexMat[vertex].eqNum;
                    
                    delete [] VertexMat[vertex].CaVec[0];
                    delete [] VertexMat[vertex].CbVec[0];
                    
                    delete [] VertexMat[vertex].CaVec;
                    delete [] VertexMat[vertex].CbVec;
                }// end of if flag
                
            }// end of side0
        }// end of side1
    }// end of side2
}// end of freeVertexVariables

/// \brief a function to get the starting and ending indices of the boundary range
/// \param bdryRange (input/output): carries the indices at each axis and side
/// \param indexRange (input): carries the range of indices of the spatial grid without the ghost points
/// \param fixedAxis (input): the axis that is fixed at the given boundary
/// \param fixedSide (input): the side that is fixed at the given boundary
/// \param addOnSide0 (input): an option to add a number of grid points at side 0 of the boundary
/// \param addOnSide1 (input): an option to add a number of grid points at side 1 of the boundary
/// \param dim (input): the dimension in space
void getBdryRange(int bdryRange[3][2], int indexRange[3][2], int fixedAxis, int fixedSide, int addOnSide0, int addOnSide1, int dim){
    for(int axis = 0; axis<3; axis++){
        for(int side = 0; side<2; side++){
            if(axis == fixedAxis){
                bdryRange[axis][side] = indexRange[fixedAxis][fixedSide];
            }
            else{
                if(axis<dim){
                    bdryRange[axis][side] = indexRange[axis][side] + (1-side)*addOnSide0 + addOnSide1*side;
                }else{
                    bdryRange[axis][side] = indexRange[axis][side];
                }
            }// end of if axis
        }// end of side
    }// end of axis
}// end of getDataBdryRange

/// \brief A function to extrapolate the values of the solution at the ghost points
/// \param unp1 (input/output): the grid solution at the new time
/// \param mask (input): the mask to determine type of grid points
/// \param unp1Dim (input): the number of grid points in each spatial dimension of the grid solution
/// \param indexRange (input): the range of indices of the spatial grid
/// \param faceEval (input): carries the type of boundary at each face (defined above)
/// \param numGhost (input): number of ghost points needed
/// \param userNumGhost (input): number of ghost points utilized by the user
/// \param p (input): orderInSpace/2
/// \param dim (input): dimension in space
void extrapolateAllGhost(double *&unp1, int *mask, int unp1Dim[3], int indexRange[3][2], int *faceEval, int numGhost, int userNumGhost, int p, int dim){
    int order = 2*p;
    int i[3];
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = side + 2*axis;
            
            if(faceEval[face]>0){
                int bdryRange[3][2]; getBdryRange(bdryRange, indexRange, axis, side, (-p), (p), dim);
                
                for(i[2] = bdryRange[2][0]; i[2]<=bdryRange[2][1]; i[2]++){
                    for(i[1] = bdryRange[1][0]; i[1]<=bdryRange[1][1]; i[1]++){
                        for(i[0] = bdryRange[0][0]; i[0]<=bdryRange[0][1]; i[0]++){
                            
                            int index[3] = {i[0],i[1],i[2]};
                            
                            if(side == 0){
                                for(int j = (-1); j>=(-p); j--){
                                    index[axis] = i[axis] + j;
                                    
                                    if(mask[solInd(index,unp1Dim)]>0)
                                        unp1[solInd(index,unp1Dim)] = extrapolate_data(unp1, index, unp1Dim, axis, (1), order);
                                }
                            }
                            else{
                                for(int j=1; j<=p; j++){
                                    index[axis] = i[axis] + j;
                                    
                                    if(mask[solInd(index,unp1Dim)]>0)
                                        unp1[solInd(index,unp1Dim)] = extrapolate_data(unp1, index, unp1Dim, axis, (-1), order);
                                }
                            }
                            
                        }// end of i0 loop
                    }// end of i1 loop
                }// end of i2 loop
            }
        }// end of side loop
    }// end of axis loop
}// end of setBCsNEW

/// \brief A function to set the exact values of the solution on interpolation boundaries used for testing and debugging
/// \param un (input/output): the solution grid function
/// \param mask (input): the mask to determine type of grid points
/// \param ue (input): grid function carrying the exact solution values on the grid
void Lcbc::setExactGhost(double *&un, double *ue, int *mask){
    int i[3];
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            
            int face = side + 2*axis;
            if(faceEval[face] == 0){
                int j0[2] = {sideBasedValue(side, (-p), 0),sideBasedValue(side, 0, p)};
                
                if(!(faceParam[face].exists)){
                    faceParam[face].initialize(axis, side, mask, numGhost, userNumGhost, G.indexRange, G.Ngx, p, dim, faceEval, addAuxEqns);
                }
                
                for(i[2] = faceParam[face].solnBdryRange[2][0]; i[2]<=faceParam[face].solnBdryRange[2][1]; i[2]++){
                    for(i[1] = faceParam[face].solnBdryRange[1][0]; i[1]<=faceParam[face].solnBdryRange[1][1]; i[1]++){
                        for(i[0] = faceParam[face].solnBdryRange[0][0]; i[0]<=faceParam[face].solnBdryRange[0][1]; i[0]++){
                            
                            for(int k0 = j0[0]; k0<=j0[1]; k0++){
                                int I[3] = {i[0],i[1],i[2]};
                                I[axis] = i[axis] + k0;
                                un[solInd(I,G.Ngx)] = ue[solInd(I,G.Ngx)];
                            }
                        }// end of i0 loop
                    }// end of i1 loop
                }// end of i2 loop
            }
        }// end of side loop
    }// end of axis loop
}// end of setBCsNEW

/// \brief Update the bdryTypeParams objects. These are objects that carry information which depends on the boundary type.
/// The boundary type is determined based on the number of fixed axes and the number of known axes.
/// Fixed axes: axes that are fixed on a given boundary
/// Known axes: fixed axes where the boundary data is known from physical boundary conditions
void Lcbc::updateNeededParameters(){
    
    int n = (2*p+1);
    
    /* one fixed axis */
    oneAxis.knownAxesNum = 1;
    oneAxis.ghostPointNum = p;
    oneAxis.interiorEqNum = (p+1)*n*dimBasedValue(dim, 1, n);
    oneAxis.unknownVarNum = (totalVarNum - oneAxis.interiorEqNum);
    oneAxis.fixedAxesNum = 1;

    /* two fixed axes */
    twoAxis.fixedAxesNum = 2;
    twoAxis.knownAxesNum = 2;
    twoAxis.interiorEqNum = (p+1)*(p+1)*dimBasedValue(dim, 1,n);
    twoAxis.unknownVarNum = totalVarNum - twoAxis.interiorEqNum;
    twoAxis.ghostPointNum = twoAxis.unknownVarNum - 2*p;


    /* three fixed axes */
    threeAxis.fixedAxesNum = 3;
    threeAxis.knownAxesNum = 3;
    threeAxis.interiorEqNum = (p+1)*(p+1)*(p+1);
    threeAxis.unknownVarNum = totalVarNum - threeAxis.interiorEqNum;
    threeAxis.ghostPointNum = threeAxis.unknownVarNum - 3*p;

    
    /* two fixed axes and one known axes */
    nearEdge.fixedAxesNum = 2;
    nearEdge.knownAxesNum = 1;
    nearEdge.ghostPointNum = p*p;
    nearEdge.interiorEqNum = (p+1)*n*dimBasedValue(dim, 1, n);
    nearEdge.unknownVarNum = (totalVarNum - nearEdge.interiorEqNum);

    nearVertex[0].fixedAxesNum = 3;
    nearVertex[0].knownAxesNum = 1;
    nearVertex[0].ghostPointNum = p*p*p;
    nearVertex[0].interiorEqNum = (p+1)*n*dimBasedValue(dim, 1, n);
    nearVertex[0].unknownVarNum = (totalVarNum - nearVertex[0].interiorEqNum);
    
    nearVertex[1].fixedAxesNum = 3;
    nearVertex[1].knownAxesNum = 2;
    nearVertex[1].interiorEqNum = (p+1)*(p+1)*dimBasedValue(dim, 1,n);
    nearVertex[1].unknownVarNum = (totalVarNum - nearVertex[1].interiorEqNum);
    nearVertex[1].ghostPointNum = twoAxis.ghostPointNum - p*(p+1);
}

/// \brief A function to count the number of faces where the boundary data is known at a 3D vertex
/// \param vec (input): a vector that contains the boundary types (physical, periodic or interpolation) at each face 
int countKnownBdryFaces(int vec[3]){
    int known = 0;
    
    for(int face = 0; face<3; face++){
        if(vec[face]>0){
            known++;
        }
    }
    return known;
}
