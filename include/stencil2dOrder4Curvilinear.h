! Stencil: nd=2, orderOfAccuracy=4, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
 + sc(  1,i1,i2)*u(i1-2,i2-2,i3,m) + sc(  2,i1,i2)*u(i1-1,i2-2,i3,m) + sc(  3,i1,i2)*u(i1+0,i2-2,i3,m) + sc(  4,i1,i2)*u(i1+1,i2-2,i3,m) + sc(  5,i1,i2)*u(i1+2,i2-2,i3,m)\
 + sc(  6,i1,i2)*u(i1-2,i2-1,i3,m) + sc(  7,i1,i2)*u(i1-1,i2-1,i3,m) + sc(  8,i1,i2)*u(i1+0,i2-1,i3,m) + sc(  9,i1,i2)*u(i1+1,i2-1,i3,m) + sc( 10,i1,i2)*u(i1+2,i2-1,i3,m)\
 + sc( 11,i1,i2)*u(i1-2,i2+0,i3,m) + sc( 12,i1,i2)*u(i1-1,i2+0,i3,m) + sc( 13,i1,i2)*u(i1+0,i2+0,i3,m) + sc( 14,i1,i2)*u(i1+1,i2+0,i3,m) + sc( 15,i1,i2)*u(i1+2,i2+0,i3,m)\
 + sc( 16,i1,i2)*u(i1-2,i2+1,i3,m) + sc( 17,i1,i2)*u(i1-1,i2+1,i3,m) + sc( 18,i1,i2)*u(i1+0,i2+1,i3,m) + sc( 19,i1,i2)*u(i1+1,i2+1,i3,m) + sc( 20,i1,i2)*u(i1+2,i2+1,i3,m)\
 + sc( 21,i1,i2)*u(i1-2,i2+2,i3,m) + sc( 22,i1,i2)*u(i1-1,i2+2,i3,m) + sc( 23,i1,i2)*u(i1+0,i2+2,i3,m) + sc( 24,i1,i2)*u(i1+1,i2+2,i3,m) + sc( 25,i1,i2)*u(i1+2,i2+2,i3,m)\
 FV(m)
