#include "LCBC.h"
#include "LCBCmacros.h"

void Lcbc::updateGhost(double *unp1, double t, double dt){

    int n = (2*p+1);
    int NU = (p+1);
    int auxiliaryEqNum  = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int compCondNum     = (p+1)*n*dimBasedValue(dim,1,n);
    int approxEqNum     = compCondNum + auxiliaryEqNum;
    int faceNum = 2*dim;
    
    double *R[(NU*faceNum)];
    
    /* prepare the derivatives of data functions over the stencil and update face ghost */
    for(int axis = 0; axis<dim; axis++){
        for(int side = 0; side<2; side++){
            prepDataDeriv(R, t, dt, axis, side);
            updateFaceGhost(R, unp1, t, dt, approxEqNum, axis, side);
        }// end side loop
    }// end axis loop
    
    /* Update the 2D corners or 3D edges */
    approxEqNum = 2*compCondNum + auxiliaryEqNum;

    for(int varAxis = dimBasedValue(dim, 2, 0); varAxis<3; varAxis++){
        for(int side1 = 0; side1<2; side1++){
            for(int side2 = 0; side2<2; side2++){
                updateEdgeGhost(R, unp1, t, dt, approxEqNum, side1, side2, varAxis);
            }// end of side2
        }// end of side1
    }// end of varAxis loop
    
    /* Update the 3D vertices */
    if(dim == 3){
        approxEqNum = 3*compCondNum + auxiliaryEqNum;
        for(int side2 = 0; side2<2; side2++){
            for(int side1 = 0; side1<2; side1++){
                for(int side0 = 0; side0<2; side0++){
                    updateVertexGhost(R, unp1, t, dt, approxEqNum, side0, side1, side2);
                }// end of side0 loop
            }// end of side1 loop
        }// end of side2 loop
    }// end of if statement
    
    /* free the R variable */
    for(int nu = 0; nu<(NU*faceNum); nu++){
        delete [] R[nu];
        R[nu] = NULL;
    }
}

void Lcbc::getLagrangeData(){
    /* This function evaluate L_{i}(z) where z in [-w,w] and i in [-p,p] */
    int n = 2*p+1;
    int w = 2*p;
    int m = 2*w+1;
    LagrangeData = (double*) malloc((n*m)*sizeof(double));

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
    /* arg = {Coef#,x,y}, Coef# = {0,1,2,3,4,5} for {c11,c12,c22,c1,c2,c0} */
    
    double C = 0;
    if(arg[0] == 0){
        C = 1;
    }else if(arg[0] == 2){
        C = 1;
    }else{
        C = 0;
    }
    return C;
}

void Lcbc::getCoef(Grid& G){
    coef = (double **) malloc(maxCoefNum*sizeof(double *)); 
    for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
        coef[coefNum] = (double *) malloc(G.Ng*sizeof(double));
    }
    int i[3];
    for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
        for(i[2] = 0; i[2]<G.Ngx[2]; i[2]++){
            for(i[1] = 0; i[1]<G.Ngx[1]; i[1]++){
                for(i[0] = 0; i[0]<G.Ngx[0]; i[0]++){
                    double arg[] = {((double) coefNum),G.x[0][i[0]],G.x[1][i[1]],G.x[2][i[2]]};
                    coef[coefNum][ind(i,G.Ngx)] = Cn(arg);
                }// end of i3 loop
            }// end of i2 loop
        }// end of i0 loop
    }// end of coefNum loop
}// end of getCoef function

void Lcbc::freeVariables(){
    free(LagrangeData); LagrangeData = NULL;
    for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
        free(coef[coefNum]); coef[coefNum] = NULL;
    }
    free(coef); coef = NULL;
}

void Lcbc::freeFaceVariables(){
    for(int axis = 0; axis<dim; axis++){
        int bdryNg = getBdryPointsNum(axis);
        for(int side = 0; side<2; side++){
            int face = ind2(side,axis,2,3);
            if(FaceMat[face].flag){
                for(int bdryPoint = 0; bdryPoint<bdryNg; bdryPoint++){
                    for(int vecNum = 0; vecNum<p; vecNum++){
                        free(FaceMat[face].CaVec[bdryPoint][vecNum]); FaceMat[face].CaVec[bdryPoint][vecNum] = NULL;
                        free(FaceMat[face].CbVec[bdryPoint][vecNum]); FaceMat[face].CbVec[bdryPoint][vecNum] = NULL;
                    }
                    free(FaceMat[face].CaVec[bdryPoint]); FaceMat[face].CaVec[bdryPoint] = NULL;
                    free(FaceMat[face].CbVec[bdryPoint]); FaceMat[face].CbVec[bdryPoint] = NULL;
                }
                
                delete [] FaceMat[face].eqNum;
                 free(FaceMat[face].CaVec);
                 free(FaceMat[face].CbVec);
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
                            free(CornerMat[corner].CaVec[bdryPt][gp]);
                            free(CornerMat[corner].CbVec[bdryPt][gp]);
                        }// end of gp loop
                        free(CornerMat[corner].CaVec[bdryPt]);
                        free(CornerMat[corner].CbVec[bdryPt]);
                    }// end of bdryPt loop
                    
                    free(CornerMat[corner].CaVec);
                    free(CornerMat[corner].CbVec);
                    
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
