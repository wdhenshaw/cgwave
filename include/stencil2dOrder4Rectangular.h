! Stencil: nd=2, orderOfAccuracy=4, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
                                                           + scr(  3)*u(i1+0,i2-2,i3,m)                                                          \
                              + scr(  7)*u(i1-1,i2-1,i3,m) + scr(  8)*u(i1+0,i2-1,i3,m) + scr(  9)*u(i1+1,i2-1,i3,m)                             \
 + scr( 11)*u(i1-2,i2+0,i3,m) + scr( 12)*u(i1-1,i2+0,i3,m) + scr( 13)*u(i1+0,i2+0,i3,m) + scr( 14)*u(i1+1,i2+0,i3,m) + scr( 15)*u(i1+2,i2+0,i3,m)\
                              + scr( 17)*u(i1-1,i2+1,i3,m) + scr( 18)*u(i1+0,i2+1,i3,m) + scr( 19)*u(i1+1,i2+1,i3,m)                             \
                                                           + scr( 23)*u(i1+0,i2+2,i3,m)                                                          \
 FV(m)
