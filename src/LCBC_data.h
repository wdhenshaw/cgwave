#ifndef LCBC_data_h
#define LCBC_data_h

#include <stdio.h>
#include <stdlib.h>

/// \brief In this file, we define the LcbcData class. Each LcbcData object holds known data values such as PDE coefficients, forcing and boundary data.

/// The LcbcData object carries known data values in a grid function called Fn. Fn can have components such as coefficient numbers, number of time derivatives, etc
/// Fn[component][i] where i = [i,j,k] (indices on the grid)

/// The LcbcData object also carries the following information
/// lth[3] is the length of Fn in each axis.
/// wth[3] is the number of ghost points in each axis.
/// indexLth is the number of components
///
/// Example coef[coefNum][i]
/// Example bdryData[nu][i] where nu = 0, 1, ..., p represents time derivatives

class LcbcData{
public:
    void initialize( int lthIn[3], int wthIn[3], int indexLthIn, double **FnIn = NULL, bool fixIn=false);
    void allocateSpace(int lth[3], int indexLth);
    void deleteLcbcData();
    
    int lth[3];   // the length of the grid function Fn in each axis.
    int wth[3];   // the number of ghost points in each axis.
    int indexLth; // the number of components
    
    bool allocatedSpace = false; // to tell if space has been allocated or not
    bool fix            = false; // to tell if the data needs to be extended to more grid points
    bool fill           = false; // to tell if user input data needs to be refilled in another grid function
    bool initialized    = false; // to tell if the LcbcData object is initialized
    double **Fn; // the grid function that holds boundary and forcing data on the grid
    
    /* Constructor */
    LcbcData(){
    }
    
    /* Destructor */
    ~LcbcData(){
        deleteLcbcData();
    }
};

void checkUserData(LcbcData *&classData, LcbcData *userData, int Nx[3], int Ngx[3], int faceEval[], int dim, int p, int extraDataGhost, int dataType, bool fix = false, bool overRide = false);
void checkData(LcbcData &classData, LcbcData &userData, int Nx[3], int Ngx[3], int axis, int p, int dim, int extraDataGhost, int FnType, bool overRide);
void fixData(LcbcData &data, int p, int dim, int axis);
void fillData(LcbcData &classData, LcbcData &userData);
int compareVectors(int vec[3], int cstVec[3], int minVec[3]);
void getBdryLcbcLthWth(int (&lcbcLth)[3], int (&lcbcWth)[3], int (&minLth)[3], int Ngx[3], int fixedAxis, int p, int dim);
void getForcingLthWth(int Nx[3], int (&lcbcLth)[3], int (&lcbcWth)[3], int (&minLth)[3], int fixedAxis, int dim, int p, int extraDataGhost);
void getVarAxis_data(int varAxis[2], int fixedAxis);
int compareVectors(int vec[3], int cstVec[3], int minVec[3]);
double extrapolate_data(double *gridFn, int index[3], int lth[3], int axis, int sign, int order);

#endif 
