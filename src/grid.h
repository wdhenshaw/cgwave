#ifndef grid_h
#define grid_h

class Grid{
public:
    static const int MaxDim = 3;
    double dx[MaxDim], *x[MaxDim];
    int indexRange[MaxDim][2], Ngx[MaxDim], Nx[MaxDim], dim, Ng;
    
    Grid(int *numGridPoints, int numGhost, int dimIn){
        dim = dimIn;
        Ng = 1;
        for(int axis = 0; axis<dim; axis++){
            Nx[axis]  = numGridPoints[axis];
            Ngx[axis] = 2*numGhost + numGridPoints[axis] + 1;
            dx[axis]  = (1.0/(double)numGridPoints[axis]);
            x[axis]   = new double[Ngx[axis]];
            Ng = Ng*Ngx[axis];
            
            /* save the indices where physical boundary begins and ends */
            for(int side = 0; side<=1; side++){
                indexRange[axis][side] = numGhost + side*numGridPoints[axis];
            }
            
            /* create the variables on the grid */
            for(int side = 0; side<Ngx[axis]; side++){
                x[axis][side] = ((double)(side - numGhost))*dx[axis];
            }// end of side loop
            
        }// end of axis loop
        
        if(MaxDim>dim){
            for(int k = dim; k<MaxDim; k++){
                for(int side = 0; side<=1; side++)
                    indexRange[k][side] = 0;
                
                x[k]    =  new double[1];
                x[k][0] = 0;
                dx[k]   = 0;
                Ngx[k]  = 1;
            }
        }
    }// end of constructor of Grid class
    
    Grid(){
    }
    
    void prepareGrid(int *numGridPoints, int numGhost, int dimIn){
        dim = dimIn;
        Ng = 1;
        for(int axis = 0; axis<dim; axis++){
            Nx[axis]  = numGridPoints[axis];
            Ngx[axis] = 2*numGhost + numGridPoints[axis] + 1;
            dx[axis]  = (1.0/(double)numGridPoints[axis]);
            x[axis]   = new double[Ngx[axis]];
            Ng = Ng*Ngx[axis];
            
            /* save the indices where physical boundary begins and ends */
            for(int side = 0; side<=1; side++){
                indexRange[axis][side] = numGhost + side*numGridPoints[axis];
            }
            
            /* create the variables on the grid */
            for(int side = 0; side<Ngx[axis]; side++){
                x[axis][side] = ((double)(side - numGhost))*dx[axis];
            }// end of side loop
            
        }// end of axis loop
        
        if(MaxDim>dim){
            for(int k = dim; k<MaxDim; k++){
                for(int side = 0; side<=1; side++)
                    indexRange[k][side] = 0;
                
                x[k]    = new double[1];
                x[k][0] = 0;
                dx[k]   = 0;
                Ngx[k]  = 1;
            }
        }
    }// end of prepareGrid
    
    ~Grid(){
        for(int axis = 0; axis<MaxDim; axis++){
            delete [] x[axis]; x[axis] = NULL;
        }
    }// end of destructor
};// end of class

#endif /* grid_h */
