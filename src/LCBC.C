#include "LCBC.h"
#include "utility.h"
#include "LCBCmacros.h"

bool debugc = false;

void Default_bn(double bn[3], double *arg);

void Lcbc::initialize(int dimIn, int orderInSpace, int orderInTimeIn, int *numGridPoints, int numGhostIn, int *faceEvalIn, LcbcData *coefIn, F1 CnIn, F1 GnIn, F1 FnIn, F2 bnIn, bool cstCoefIn, bool *zeroBCIn, bool noForcingIn){
    
    isInitialized = true;
    
    p = orderInSpace/2;
    extraDataGhost = (p-1);
    
    numGhost = p;
    if(numGhostIn<numGhost){
        printf("ERROR: need at least %d ghost points to use LCBC\n",numGhost);
        exit(-1);
    }
    else{
        userNumGhost = numGhostIn;
    }
    
    orderInTime = orderInTimeIn;
    dim = dimIn;
    param.initialize(dim, p, addAuxEqns);
    getLagrangeData();
    G.prepareGrid(numGridPoints,numGhost,dim);
    
    /* Set faceEval */
    
    faceEval = faceEvalIn;
    if(faceEval == NULL){
        faceEvalNewed = true;
        faceEval = new int[(2*dim)];
        for(int face = 0; face<(2*dim); face++){
            faceEval[face] = 1;
        }// end of face loop
    }// end of faceEval condition
    
    /* set the forcing function */
    if(FnIn == NULL){
        Fn = Default_Fn;
        defaultedFn = true;
    }else{
        Fn = FnIn;
    }
    noForcing = noForcingIn;
    
    /* set the BC functions */
    if(GnIn == NULL){
        Gn = Default_Gn;
        defaultedGn = true;
    }else{
        Gn = GnIn;
    }
    if(bnIn == NULL){
        bn = Default_bn;
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
    maxCoefNum = ((dim+1)*(dim+2)/2);
    coef = new LcbcData[(2*dim)];
    if(coefIn == NULL){
        if(CnIn == NULL){
            Cn = Default_Cn;
            cstCoef = true;
        }else{
            Cn = CnIn;
            cstCoef = cstCoefIn;
        }
        getCoef();
    }else{
        cstCoef = cstCoefIn;
        bool fix = true;
        bool overRideCheck = (cstCoef)?(true):(false);
        checkUserData(coef, coefIn, G.Nx, G.Ngx, faceEval, dim, p, extraDataGhost, coefficient, fix, overRideCheck);
    }// end of else
    
    /* prepare the LCBC data functions */
    forcingData = new LcbcData[(2*dim)];
    bdryData = new LcbcData[(2*dim)];
    
    /* prepare memory for the 8th-order accurate scheme here */
    if(dim == 3){
        memoryComponents = 32 + 26;
        int componentSize = (4*p + 1)*(4*p + 1)*(2*p + 1);
        preallocateMemory(memoryComponents, componentSize);
        preallocatedMemory = true;
    }
    
}// end of initialize

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


void Lcbc::updateGhost(double *&unp1, double t, double dt, LcbcData *&gn, LcbcData *&fn){
    
    if(!analyzedUserData)
        analyzeUserInput(gn, fn);
    /* ================================================== */

    int NU = param.NU;
    int faceNum = param.faceNum;
    double **R = new double*[(NU*faceNum)];
    
    /* prepare the derivatives of data functions over the stencil and update face ghost */
    
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                
                if(initializedForcingData){
                    if(forcingData[face].fill){
                        if(debugc){printf("filling Fn face %d\n",face); }
                        fillData(forcingData[face], fn[face]);
                    }else{
                        forcingData[face].Fn = fn[face].Fn;
                    }
                    
                    if(forcingData[face].fix){
                        if(debugc){printf("fixing Fn face %d\n",face); }
                        fixData(forcingData[face], p, dim, axis);
                    }
                }
                if(initializedBdryData && (!zeroBC[face])){
                    bdryData[face].Fn = gn[face].Fn;
                }
                
                if(!(faceParam[face].exists))
                    faceParam[face].initialize(axis, side, G.indexRange, G.Ngx, p, dim, faceEval, addAuxEqns);
                
                getDataDeriv(R, t, dt, axis, side, bdryData[face], forcingData[face]);
                updateFaceGhost(unp1, R, t, dt, axis, side);
                
            }// end if faceEval
        }// end side loop
    }// end axis loop
    
//     Update the 2D corners or 3D edges
        for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<3; varAxis++){

            int fixedAxis[2]; getFixedAxis(fixedAxis, varAxis);

            for(int side1 = 0; side1<2; side1++){
                for(int side2 = 0; side2<2; side2++){
                    int edgeFace1 = side1 + 2*fixedAxis[0];
                    int edgeFace2 = side2 + 2*fixedAxis[1];

                    /* Case of an edge along non-periodic faces */
                    if((faceEval[edgeFace1]>0) && (faceEval[edgeFace2]>0)){
                        updateEdgeGhost(R, unp1, t, dt, side1, side2, varAxis);
                    }

                    /* Case of an edge along a periodic face */
                    else if((faceEval[edgeFace1]*faceEval[edgeFace2])==(-1)){
                        updateEdgeGhostPeriodic(unp1, side1, side2, varAxis, fixedAxis);
                    }
                }// end of side2
            }// end of side1
        }// end of varAxis loop

//     Update the 3D vertices
    if(dim == 3){
        for(int side2 = 0; side2<2; side2++){
            for(int side1 = 0; side1<2; side1++){
                for(int side0 = 0; side0<2; side0++){
                        int cornerFace0 = side0;
                        int cornerFace1 = side1 + 2;
                        int cornerFace2 = side2 + 4;
                        /* Case of DDD vertex */
                        if((faceEval[cornerFace0]>0) && (faceEval[cornerFace1]>0) && (faceEval[cornerFace2]>0)){
                            updateVertexGhost(R, unp1, t, dt, side0, side1, side2);
                        }
                        /* Case of a vertex with a periodic face */
                        else if((faceEval[cornerFace0]==-1)||(faceEval[cornerFace1]==-1)||(faceEval[cornerFace2]==-1)){
                            updateVertexGhostPeriodic(unp1, side0, side2, side2);
                        }
                }// end of side0 loop
            }// end of side1 loop
        }// end of side2 loop
    }// end of if statement
    
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

void Lcbc::getLagrangeData(){
    /* This function evaluate L_{i}(z) where z in [-w,w] and i in [-p,p] */
    int n = 2*p+1;
    int w = 2*p;
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

double Lcbc::Default_Cn(double *arg){
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

double Lcbc::Default_Gn(double *arg){
    /* This is the default boundary conditions functions */
    /* arg = {side,x,y,t}, side = {0,1,2,3} for {left,bottom,right,top} */
    
    double G = 0;
    return G;
}

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

double Lcbc::Default_Fn(double *arg){
    /* This is the default forcing function */
    /* arg = {0,x,y,z,t}                    */
    
    double F = 0;
    return F;
}

void Lcbc::getCoef(){
    int faceNum = param.faceNum;
    int maxCoefNum = ((dim+1)*(dim + 2))/2;

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

void Lcbc::getCoefGridLth(int indexRange[3][2], int lth[3], int wth[3], int fixedAxis, int fixedSide, int dim, int p){
    int bdryRange[3][2];
    getCoefGridLth(indexRange, lth, wth, bdryRange, fixedAxis, fixedSide, dim, p);
}

/// preallocateMemory: allocates space in 'memory' double pointer to be used in the interior schemes to prevent dynamic allocation at each time-step
void Lcbc::preallocateMemory(int components, int componentSize){
    memory = new char*[components];
    for(int i = 0; i<components; i++){
        memory[i] = new char[(sizeof(double)*componentSize)];
    }
}

/// deleteMemory: deletes the preallocated memory when time-stepping is finished
void Lcbc::deleteMemory(int components){
    for(int i = 0; i<components; i++){
       delete [] memory[i];
    }// end of i loop
    
    delete [] memory;
}



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

void Lcbc::freeCornerVariables(){
    
    int interiorEqNum = (p+1)*(p+1)*dimBasedValue(dim, 1, (2*p+1));
    int totalVarNum   = param.totalVarNum;
    int unknownVarNum = totalVarNum - interiorEqNum;
    
    for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<=2; varAxis++){
        
        int bdryNg = (cstCoef)?(1):(dimBasedValue(dim, 1, G.Ngx[varAxis]));
        
        for(int side1 = 0; side1<2; side1++){
            for(int side2 = 0; side2<2; side2++){
                
                int corner = ind3(side1,side2,dimBasedValue(dim, 0, varAxis),2,2,3);
                
                if(CornerMat[corner].flag){
                    
                    for(int bdryPt = 0; bdryPt<bdryNg; bdryPt++){
                        delete [] CornerMat[corner].CaVec[bdryPt];
                        delete [] CornerMat[corner].CbVec[bdryPt];
                    }// end of bdryPt loop
                    
                    delete [] CornerMat[corner].CaVec;
                    delete [] CornerMat[corner].CbVec;
                    
                    delete [] CornerMat[corner].eqNum;
                }// end of if flag
            }// end of side2
        }// end of side1
    }// end of varAxis
}// end of freeCornerVariables

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

