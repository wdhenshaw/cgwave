! Stencil: nd=2, orderOfAccuracy=2, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
                              + scr(  2)*u(i1+0,i2-1,i3,m)                             \
 + scr(  4)*u(i1-1,i2+0,i3,m) + scr(  5)*u(i1+0,i2+0,i3,m) + scr(  6)*u(i1+1,i2+0,i3,m)\
                              + scr(  8)*u(i1+0,i2+1,i3,m)                             \
 FV(m)
