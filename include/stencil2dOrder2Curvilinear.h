! Stencil: nd=2, orderOfAccuracy=2, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
 + sc(  1,i1,i2)*u(i1-1,i2-1,i3,m) + sc(  2,i1,i2)*u(i1+0,i2-1,i3,m) + sc(  3,i1,i2)*u(i1+1,i2-1,i3,m)\
 + sc(  4,i1,i2)*u(i1-1,i2+0,i3,m) + sc(  5,i1,i2)*u(i1+0,i2+0,i3,m) + sc(  6,i1,i2)*u(i1+1,i2+0,i3,m)\
 + sc(  7,i1,i2)*u(i1-1,i2+1,i3,m) + sc(  8,i1,i2)*u(i1+0,i2+1,i3,m) + sc(  9,i1,i2)*u(i1+1,i2+1,i3,m)\
 FV(m)
