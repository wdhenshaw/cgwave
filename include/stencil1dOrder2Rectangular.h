! Stencil: nd=1, orderOfAccuracy=2, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
 + scr(  1)*u(i1-1,i2-1,i3+0,m) + scr(  2)*u(i1+0,i2-1,i3+0,m) + scr(  3)*u(i1+1,i2-1,i3+0,m)\
