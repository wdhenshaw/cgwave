! Stencil: nd=2, orderOfAccuracy=6, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
                                                                                        + scr(  4)*u(i1+0,i2-3,i3,m)                                                                                       \
                                                           + scr( 10)*u(i1-1,i2-2,i3,m) + scr( 11)*u(i1+0,i2-2,i3,m) + scr( 12)*u(i1+1,i2-2,i3,m)                                                          \
                              + scr( 16)*u(i1-2,i2-1,i3,m) + scr( 17)*u(i1-1,i2-1,i3,m) + scr( 18)*u(i1+0,i2-1,i3,m) + scr( 19)*u(i1+1,i2-1,i3,m) + scr( 20)*u(i1+2,i2-1,i3,m)                             \
 + scr( 22)*u(i1-3,i2+0,i3,m) + scr( 23)*u(i1-2,i2+0,i3,m) + scr( 24)*u(i1-1,i2+0,i3,m) + scr( 25)*u(i1+0,i2+0,i3,m) + scr( 26)*u(i1+1,i2+0,i3,m) + scr( 27)*u(i1+2,i2+0,i3,m) + scr( 28)*u(i1+3,i2+0,i3,m)\
                              + scr( 30)*u(i1-2,i2+1,i3,m) + scr( 31)*u(i1-1,i2+1,i3,m) + scr( 32)*u(i1+0,i2+1,i3,m) + scr( 33)*u(i1+1,i2+1,i3,m) + scr( 34)*u(i1+2,i2+1,i3,m)                             \
                                                           + scr( 38)*u(i1-1,i2+2,i3,m) + scr( 39)*u(i1+0,i2+2,i3,m) + scr( 40)*u(i1+1,i2+2,i3,m)                                                          \
                                                                                        + scr( 46)*u(i1+0,i2+3,i3,m)                                                                                       \
 FV(m)
