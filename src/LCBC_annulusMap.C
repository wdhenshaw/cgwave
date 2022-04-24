#include "LCBC_annulusMap.h"
#include <math.h>
#define PI M_PI

extern OGFunction *ogFun;

double R0 = 0.5;
double R1 = 1.0;
double theta0 = 0.0;
double theta1 = 2*PI;
double gam = R1 - R0;
double psi = theta1 - theta0;

/* mapping information */
double G1(double r, double s){ double g1 = ( (gam*s + R0)*cos(psi*r + theta0)); return g1;}
double G2(double r, double s){ double g2 = ( (gam*s + R0)*sin(psi*r + theta0)); return g2;}
double sx(double r, double s){
    double D = ( cos(psi*r + theta0)/gam);
    return D;
}
double sy(double r, double s){
    double D = ( sin(psi*r + theta0)/gam);
    return D;
}
double rx(double r, double s){
    double D = (-sin(psi*r + theta0)/(psi*(gam*s + R0)));
    return D;
}
double ry(double r, double s){
    double D = ( cos(psi*r + theta0)/(psi*(gam*s + R0)));
    return D;
}

double sxx(double r, double s){
    double D = -(psi/gam)*sin(psi*r + theta0)*rx(r, s);
    return D;
}

double syy(double r, double s){
    double D = (psi/gam)*cos(psi*r + theta0)*ry(r,s);
    return D;
}// end of ryy

double sxy(double r, double s){
    double D = (-psi/gam)*sin(psi*r + theta0)*ry(r,s);
    return D;
}// end of rxy

double rxx(double r, double s){
    double D = (gam/psi)*sin(psi*r + theta0)*(1/((gam*s + R0)*(gam*s + R0)))*sx(r,s)
    - (1/(gam*s + R0))*cos(psi*r + theta0)*rx(r,s);
    return D;
}

double ryy(double r, double s){
    double D = (-gam/psi)*cos(psi*r + theta0)*(1/((gam*s + R0)*(gam*s + R0)))*sy(r,s)
    - (1/(gam*s + R0))*sin(psi*r + theta0)*ry(r,s);
    return D;
}

double rxy(double r, double s){
    double D = (gam/psi)*sin(psi*r + theta0)*(1/((gam*s + R0)*(gam*s + R0)))*sy(r,s)
    - (1/(gam*s + R0))*cos(psi*r + theta0)*ry(r,s);
    return D;
}

double annulus_Cn(double *arg){
    double r = arg[1], s = arg[2];
    double x = G1(r,s), y = G2(r,s);
    double C = 0;
    if(arg[0] == 0){
        C = rx(r,s)*rx(r,s) + ry(r,s)*ry(r,s);
    }else if(arg[0] == 1){
        C = sx(r,s)*sx(r,s) + sy(r,s)*sy(r,s);
    }else if(arg[0] == 2){
        C = rxx(r,s) + ryy(r,s); // c1, ux
    }else if(arg[0] == 3){
        C = sxx(r,s) + syy(r,s);
    }else if(arg[0] == 4){
        C = sx(r,s)*rx(r,s) + sy(r,s)*ry(r,s);
    }else if(arg[0] == 5){
        C = 0;
    }
    return C;
}

double annulus_Fn(double *arg){
    /* This is the forcing function take arg[3] = {0,x,y,t} (note that arg[0] is irrelevant to the computation) */
    double F;
    OGFunction & ue = *ogFun;
    double r = arg[1], s = arg[2], t = arg[3];
    double x = G1(r,s), y = G2(r,s);

    F = ue.t(x,y,0,0,t) - (ue.xx(x,y,0,0,t) + ue.yy(x,y,0,0,t));

    return F;
}

double annulus_Gn(double *arg){ // arranged (side,axis)
    /* This function computes the boundary functions on the 4 sides of a square domain */
    /* arg[0] = (side + 2*axis); holds the face number */
    /* arg[1] = x, arg[2] = y and arg[3] = t */

    OGFunction & ue = *ogFun;
    double r = arg[1], s = arg[2], t = arg[3], x, y;
    double G = 0;
    if(arg[0] == 0){
      x = G1(0,s), y = G2(0,s);
    }else if(arg[0] == 2){
        x = G1(r,0), y = G2(r,0);
    }else if(arg[0] == 1){
        x = G1(1,s), y = G2(1,s);
    }else if(arg[0] == 3){
        x = G1(r,1), y = G2(r,1);
    }else{
        printf("ERROR: incorrect arg[0] in annulus_Gn\n");
        exit(-1);
    }
    G = ue(x,y,0,0,t);
    return G;
}
