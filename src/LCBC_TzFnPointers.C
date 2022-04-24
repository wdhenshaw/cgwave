#include "LCBC_TzFnPointers.h"

// =================== LCBC Function Pointers ===========================
extern OGFunction *ogFun;

double Cn_2(double *arg){
    /* This is a function that holds the coefficient functions of a 2D spatial operator Q */
    /* Qu = c11*uxx + c22*uyy + c1*ux + c2*uy + 2*c12*uxy + c0*u                          */
    /* The coefficients are numbers in the order of their appearance in above statement   */
    double C;
    double n = arg[0], x = arg[1], y = arg[2]; // n is the number of coefficient function and (x,y) are the spatial coordinates.
    
    if(arg[0] == 0){
        C = 1;                    // c11
    }else if(arg[0] == 1){
        C = 1;                    // c22
    }else if(arg[0] == 2){
        C = 0.0;                    // c1
    }else if(arg[0] == 3){
        C = 0.0;                    // c2
    }else if(arg[0] == 4){
        C = 0.0;                    // c12
    }else if(arg[0] == 5){
        C = 0.0;                    // c0
    }
    
    return C;
}

double Fn_2(double *arg){
    /* This is the forcing function take arg[3] = {0,x,y,t} (note that arg[0] is irrelevant to the computation) */
    double F;
    OGFunction & ue = *ogFun;
    double x = arg[1], y = arg[2], t = arg[3];
    
    double Arg[3] = {0,x,y};
    double c11 = Cn_2(Arg);
    Arg[0] = 1; double c22 = Cn_2(Arg);
    Arg[0] = 2; double c1  = Cn_2(Arg);
    Arg[0] = 3; double c2  = Cn_2(Arg);
    Arg[0] = 4; double c12 = Cn_2(Arg);
    Arg[0] = 5; double c0  = Cn_2(Arg);
    
    F = ue.t(x,y,0,0,t) - (c11*ue.xx(x,y,0,0,t) + c22*ue.yy(x,y,0,0,t) + 2*c12*ue.xy(x,y,0,0,t) + c1*ue.x(x,y,0,0,t) + c2*ue.y(x,y,0,0,t) + c0*ue(x,y,0,0,t));
    return F;
}

double Gn_2(double *arg){ // arranged (side,axis)
    /* This function computes the boundary functions on the 4 sides of a square domain */
    /* arg[0] = (side + 2*axis); holds the face number */
    /* arg[1] = x, arg[2] = y and arg[3] = t */
    
    OGFunction & ue = *ogFun;
    double G = 0;
    double x = arg[1], y = arg[2], t = arg[3];
    
    if(arg[0] == 0){
        G = ue(0,y,0,0,t);
    }else if(arg[0] == 2){
        G = ue(x,0,0,0,t);
    }else if(arg[0] == 1){
        G = ue(1,y,0,0,t);
    }else if(arg[0] == 3){
        G = ue(x,1,0,0,t);
    }
    return G;
}

double Cn_3(double *arg){
    /* This is a function that holds the coefficient functions of a 2D spatial operator Q */
    /* Qu = c11*uxx + c22*uyy + c1*ux + c2*uy + 2*c12*uxy + c0*u                          */
    /* The coefficients are numbers in the order of their appearance in above statement   */
    double C;
    
    if(arg[0] == 0){
        C = 1.0;                    // c11
    }else if(arg[0] == 1){
        C = 1.0;                    // c22
    }else if(arg[0] == 2){
        C = 1.0;                    // c33
    }else if(arg[0] == 3){
        C = 0.0;                    // c1
    }else if(arg[0] == 4){
        C = 0.0;                    // c2
    }else if(arg[0] == 5){
        C = 0.0;                    // c3
    }else if(arg[0] == 6){
        C = 0.0;                    // c12
    }else if(arg[0] == 7){
        C = 0.0;                    // c13
    }else if(arg[0] == 8){
        C = 0.0;                    // c23
    }else if(arg[0] == 9){
        C = 0.0;                    // c0
    }
    return C;
}

double Fn_3(double *arg){
    /* This is the forcing function take arg[3] = {0,x,y,t} (note that arg[0] is irrelevant to the computation) */
    double F;
    OGFunction & ue = *ogFun;
    double x = arg[1], y = arg[2], z = arg[3], t = arg[4];
    
    double Arg[] = {0,x,y,z};
                double c11 = Cn_3(Arg);
    Arg[0] = 1; double c22 = Cn_3(Arg);
    Arg[0] = 2; double c33 = Cn_3(Arg);
    Arg[0] = 3; double c1  = Cn_3(Arg);
    Arg[0] = 4; double c2  = Cn_3(Arg);
    Arg[0] = 5; double c3  = Cn_3(Arg);
    Arg[0] = 6; double c12 = Cn_3(Arg);
    Arg[0] = 7; double c13 = Cn_3(Arg);
    Arg[0] = 8; double c23 = Cn_3(Arg);
    Arg[0] = 9; double c0  = Cn_3(Arg);
    
    F = ue.t(x,y,z,0,t) - (c11*ue.xx(x,y,z,0,t)
                         + c22*ue.yy(x,y,z,0,t)
                         + c33*ue.zz(x,y,z,0,t)
                       + 2*c12*ue.xy(x,y,z,0,t)
                       + 2*c13*ue.xz(x,y,z,0,t)
                       + 2*c23*ue.yz(x,y,z,0,t)
                         + c1*ue.x(x,y,z,0,t)
                         + c2*ue.y(x,y,z,0,t)
                         + c3*ue.z(x,y,z,0,t)
                         + c0*ue(x,y,z,0,t));
    return F;
}

double Gn_3(double *arg){ // arranged (side,axis)
    /* This function computes the boundary functions on the 4 sides of a square domain */
    /* arg[0] = (side + 2*axis); holds the face number */
    /* arg[1] = x, arg[2] = y and arg[3] = t */
    
    OGFunction & ue = *ogFun;
    double G = 0;
    double x = arg[1], y = arg[2], z = arg[3], t = arg[4];
    
          if(arg[0] == 0){
        G = ue(0,y,z,0,t);
    }else if(arg[0] == 1){
        G = ue(1,y,z,0,t);
    }else if(arg[0] == 2){
        G = ue(x,0,z,0,t);
    }else if(arg[0] == 3){
        G = ue(x,1,z,0,t);
    }else if(arg[0] == 4){
        G = ue(x,y,0,0,t);
    }else if(arg[0] == 5){
        G = ue(x,y,1,0,t);
    }
    return G;
}
