! Stencil: nd=3, orderOfAccuracy=2, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
                                        + sc(  2,i1,i2,i3)*u(i1+0,i2-1,i3-1,m)                                       \
 + sc(  4,i1,i2,i3)*u(i1-1,i2+0,i3-1,m) + sc(  5,i1,i2,i3)*u(i1+0,i2+0,i3-1,m) + sc(  6,i1,i2,i3)*u(i1+1,i2+0,i3-1,m)\
                                        + sc(  8,i1,i2,i3)*u(i1+0,i2+1,i3-1,m)                                       \
 + sc( 10,i1,i2,i3)*u(i1-1,i2-1,i3+0,m) + sc( 11,i1,i2,i3)*u(i1+0,i2-1,i3+0,m) + sc( 12,i1,i2,i3)*u(i1+1,i2-1,i3+0,m)\
 + sc( 13,i1,i2,i3)*u(i1-1,i2+0,i3+0,m) + sc( 14,i1,i2,i3)*u(i1+0,i2+0,i3+0,m) + sc( 15,i1,i2,i3)*u(i1+1,i2+0,i3+0,m)\
 + sc( 16,i1,i2,i3)*u(i1-1,i2+1,i3+0,m) + sc( 17,i1,i2,i3)*u(i1+0,i2+1,i3+0,m) + sc( 18,i1,i2,i3)*u(i1+1,i2+1,i3+0,m)\
                                        + sc( 20,i1,i2,i3)*u(i1+0,i2-1,i3+1,m)                                       \
 + sc( 22,i1,i2,i3)*u(i1-1,i2+0,i3+1,m) + sc( 23,i1,i2,i3)*u(i1+0,i2+0,i3+1,m) + sc( 24,i1,i2,i3)*u(i1+1,i2+0,i3+1,m)\
                                        + sc( 26,i1,i2,i3)*u(i1+0,i2+1,i3+1,m)                                       \
 FV(m)
