#include "LCBC.h"
#include "utility.h"
#include "LCBCmacros.h"

void Lcbc::initialize(int dimIn, int orderInSpace, int orderInTimeIn, int *numGridPoints, int numGhostIn, int *faceEvalIn, double **coefIn, F1 CnIn, F1 GnIn, F1 FnIn){
    
    isInitialized = true;
    
    p = orderInSpace/2;
    
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
        noForcing = true;
    }else{ Fn = FnIn; noForcing = false; }
    
    /* set the BC functions */
    if(GnIn == NULL){
        Gn = Default_Gn;
        zeroBC = true;
    }else{ Gn = GnIn; zeroBC = false; }
    
    /* set the PDE coefficient functions */
    maxCoefNum = ((dim+1)*(dim+2)/2);
    if(coefIn == NULL){
        if(CnIn == NULL){
            Cn = Default_Cn;
            cstCoef = true;
        }else{
            Cn = CnIn;
            cstCoef = false;
        }
        getCoef();
        allocatedCoef = true;
    }else{
        coef = coefIn;
        fixCoefGridFunctions(coef, G.indexRange, dim, p, faceEval);
        allocatedCoef = false;
    }
    
    der.expandSpace(4);
}// end of initialize

void Lcbc::updateGhost(double *&unp1, double t, double dt, double **gn, double **fn){
    
    if(!isInitialized){
        printf("ERROR: must initialize the Lcbc object first\n");
        exit(-1); 
    }
    
    if(fn!=NULL && numGhost<(2*p)){
        fixForcingGridFunctions(fn, G.indexRange, dim, p, faceEval);
    }
        
    int n = (2*p+1);
    int NU = (p+1);
    int auxiliaryEqNum  = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int compCondNum     = (p+1)*n*dimBasedValue(dim,1,n);
    int approxEqNum     = compCondNum + auxiliaryEqNum;
    int faceNum = 2*dim;

    double **R = new double*[(NU*faceNum)];
    
    /* prepare the derivatives of data functions over the stencil and update face ghost */
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                getDataDeriv(R, t, dt, axis, side, gn, fn);
                updateFaceGhost(unp1, R, t, dt, approxEqNum, axis, side);
            }// end if faceEval
        }// end side loop
    }// end axis loop
    
    /* Update the 2D corners or 3D edges */
    approxEqNum = 2*compCondNum + auxiliaryEqNum;

    for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<3; varAxis++){
        int fixedAxis[2]; getFixedAxis(fixedAxis, varAxis);
        for(int side1 = 0; side1<2; side1++){
            for(int side2 = 0; side2<2; side2++){
                int edgeFace1 = side1 + 2*fixedAxis[0];
                int edgeFace2 = side2 + 2*fixedAxis[1];
                
                /* Case of DD edge or corner */
                if((faceEval[edgeFace1]==1) && (faceEval[edgeFace2]==1)){
                    updateEdgeGhost(R, unp1, t, dt, approxEqNum, side1, side2, varAxis);
                }
                 /* Case of a DP or PD 3D edge or 2D corner */
                else if((faceEval[edgeFace1]*faceEval[edgeFace2])==(-1)){
                    updateEdgeGhostPeriodic(unp1, side1, side2, varAxis, fixedAxis);
                }

            }// end of side2
        }// end of side1
    }// end of varAxis loop

    /* Update the 3D vertices */
    if(dim == 3){
        approxEqNum = 3*compCondNum + auxiliaryEqNum;
        for(int side2 = 0; side2<2; side2++){
            for(int side1 = 0; side1<2; side1++){
                for(int side0 = 0; side0<2; side0++){
                    int cornerFace0 = side0;
                    int cornerFace1 = side1 + 2;
                    int cornerFace2 = side2 + 4;
                    /* Case of DDD vertex */
                    if((faceEval[cornerFace0]==1) && (faceEval[cornerFace1]==1) && (faceEval[cornerFace2]==1)){
                        updateVertexGhost(R, unp1, t, dt, approxEqNum, side0, side1, side2);
                    }
                    /* Case of a vertex with a periodic face */
                    else if((faceEval[cornerFace0]==-1)||(faceEval[cornerFace1]==-1)||(faceEval[cornerFace2]==-1)){
                        updateVertexGhostPeriodic(unp1, side0, side2, side2);
                    }
                    /* Case when */
                }// end of side0 loop
            }// end of side1 loop
        }// end of side2 loop
    }// end of if statement
    
    /* free the R variable */
    for(int face = 0; face<faceNum; face++){
        for(int nu = 0; nu<NU; nu++){
            
            if(faceEval[face]>0){
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

double Lcbc::Default_Fn(double *arg){
    /* This is the default forcing function */
    /* arg = {0,x,y,z,t}                    */
    
    double F = 0;
    return F;
}

void Lcbc::getCoef(){
    int faceNum = 2*dim;
    int maxCoefNum = ((dim+1)*(dim + 2))/2;

    coef = new double*[(faceNum*maxCoefNum)];

    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                int bdryRange[3][2],lth[3], wth[3]; getCoefGridLth(G.indexRange, lth, wth, bdryRange, axis, side, dim, p);
                int i[3];
                for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                    coef[ind2(face,coefNum,faceNum, maxCoefNum)] = new double[(lth[0]*lth[1]*lth[2])];

                    for(i[2] = 0; i[2]<lth[2]; i[2]++){
                        for(i[1] = 0; i[1]<lth[1]; i[1]++){
                            for(i[0] = 0; i[0]<lth[0]; i[0]++){
                                
                                double arg[] = {((double) coefNum),(G.x[0][bdryRange[0][0]] + ((- wth[0] + i[0])*G.dx[0])),
                                    (G.x[1][bdryRange[1][0]] + ((- wth[1] + i[1])*G.dx[1])),
                                    (G.x[2][bdryRange[2][0]] + ((- wth[2] + i[2])*G.dx[2]))};
                                
                                coef[ind2(face,coefNum,faceNum,maxCoefNum)][ind(i,lth)] = Cn(arg);

                            }// end of i[0] loop s
                        }// end of i[1] loop
                    }// end of i[2] loop

                }// end of coefNum
            }// end of if faceEval
        }// end of side loop
    }// end of axis loop

}// getBdryGridFunctions


void Lcbc::deleteCoef(){
    int faceNum = 2*dim;
    int maxCoefNum = ((dim+1)*(dim + 2))/2;

    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                    
                    delete [] coef[ind2(face,coefNum,faceNum, maxCoefNum)];
                    
                }// end of coefNum
            }// end of if faceEval
        }// end of side loop
    }// end of axis loop
    
    delete [] coef;

}// getBdryGridFunctions

void Lcbc::getCoefGridLth(int indexRange[3][2], int lth[3], int fixedAxis, int fixedSide, int dim, int p){
    int wth[3], bdryRange[3][2];
    for(int axis = 0; axis<3; axis++){
        wth[axis] = ((axis==fixedAxis)?(p):(2*p));
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
}

void Lcbc::getCoefGridLth(int indexRange[3][2], int lth[3], int wth[3], int bdryRange[3][2], int fixedAxis, int fixedSide, int dim, int p){
    for(int axis = 0; axis<3; axis++){
        wth[axis] = ((axis==fixedAxis)?(p):(2*p));
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
}


void Lcbc::freeVariables(){
    
    delete [] LagrangeData; LagrangeData = NULL;
    
    if(allocatedCoef == true)
        deleteCoef();

    if(faceEvalNewed)
        delete [] faceEval;
}// end of freeVariables

void Lcbc::freeFaceVariables(){
    for(int axis = 0; axis<dim; axis++){
        int bdryNg = getBdryPointsNum(axis);
        for(int side = 0; side<2; side++){
            int face = ind2(side,axis,2,3);
            if(FaceMat[face].flag){
                for(int bdryPoint = 0; bdryPoint<bdryNg; bdryPoint++){
                    for(int vecNum = 0; vecNum<p; vecNum++){
                        delete [] FaceMat[face].CaVec[bdryPoint][vecNum];
                        delete [] FaceMat[face].CbVec[bdryPoint][vecNum];
                        FaceMat[face].CaVec[bdryPoint][vecNum] = NULL;
                        FaceMat[face].CbVec[bdryPoint][vecNum] = NULL;
                    }
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
    int totalVarNum   = (2*p+1)*(2*p+1)*dimBasedValue(dim, 1, (2*p+1));
    int unknownVarNum = totalVarNum - interiorEqNum;
    int gpNum = unknownVarNum - 2*p;
    
    for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<=2; varAxis++){
        
        int bdryNg = dimBasedValue(dim, 1, G.Ngx[varAxis]);
        
        for(int side1 = 0; side1<2; side1++){
            for(int side2 = 0; side2<2; side2++){
                
                int corner = ind3(side1,side2,dimBasedValue(dim, 0, varAxis),2,2,3);
                
                if(CornerMat[corner].flag){
                    
                    for(int bdryPt = 0; bdryPt<bdryNg; bdryPt++){
                        
                        for(int gp = 0; gp<gpNum; gp++){
                            delete [] CornerMat[corner].CaVec[bdryPt][gp];
                            delete [] CornerMat[corner].CbVec[bdryPt][gp];
                        }// end of gp loop
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
    
    int totalVarNum   = (2*p+1)*(2*p+1)*(2*p+1);
    int knownVarNum   = (p+1)*(p+1)*(p+1);
    int unknownVarNum = (totalVarNum - knownVarNum);
    int ghostPointsNum = unknownVarNum - 3*p;
    
    for(int side2 = 0; side2<2; side2++){
        for(int side1 = 0; side1<2; side1++){
            for(int side0 = 0; side0<2; side0++){
                int vertex = ind3(side0,side1,side2,2,2,2);
                
                if(VertexMat[vertex].flag){
                    delete [] VertexMat[vertex].eqNum;
                    
                    for(int gp = 0; gp<ghostPointsNum; gp++){
                        delete [] VertexMat[vertex].CaVec[0][gp];
                        delete [] VertexMat[vertex].CbVec[0][gp];
                    }// end of gp loop
                    
                    delete [] VertexMat[vertex].CaVec[0];
                    delete [] VertexMat[vertex].CbVec[0];
                    
                    delete [] VertexMat[vertex].CaVec;
                    delete [] VertexMat[vertex].CbVec;
                }// end of if flag
                
            }// end of side0
        }// end of side1
    }// end of side2
}// end of freeVertexVariables
