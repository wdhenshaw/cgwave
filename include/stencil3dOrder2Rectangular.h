! Stencil: nd=3, orderOfAccuracy=2, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
                                                                                                \
                                 + scr(  5)*u(i1+0,i2+0,i3-1,m)                                \
                                                                                                \
                                 + scr( 11)*u(i1+0,i2-1,i3+0,m)                                \
 + scr( 13)*u(i1-1,i2+0,i3+0,m) + scr( 14)*u(i1+0,i2+0,i3+0,m) + scr( 15)*u(i1+1,i2+0,i3+0,m)\
                                 + scr( 17)*u(i1+0,i2+1,i3+0,m)                                \
                                                                                                \
                                 + scr( 23)*u(i1+0,i2+0,i3+1,m)                                \
                                                                                                \
 FV(m)
