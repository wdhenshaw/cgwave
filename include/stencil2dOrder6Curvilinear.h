! Stencil: nd=2, orderOfAccuracy=6, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
 + sc(  1,i1,i2)*u(i1-3,i2-3,i3,m) + sc(  2,i1,i2)*u(i1-2,i2-3,i3,m) + sc(  3,i1,i2)*u(i1-1,i2-3,i3,m) + sc(  4,i1,i2)*u(i1+0,i2-3,i3,m) + sc(  5,i1,i2)*u(i1+1,i2-3,i3,m) + sc(  6,i1,i2)*u(i1+2,i2-3,i3,m) + sc(  7,i1,i2)*u(i1+3,i2-3,i3,m)\
 + sc(  8,i1,i2)*u(i1-3,i2-2,i3,m) + sc(  9,i1,i2)*u(i1-2,i2-2,i3,m) + sc( 10,i1,i2)*u(i1-1,i2-2,i3,m) + sc( 11,i1,i2)*u(i1+0,i2-2,i3,m) + sc( 12,i1,i2)*u(i1+1,i2-2,i3,m) + sc( 13,i1,i2)*u(i1+2,i2-2,i3,m) + sc( 14,i1,i2)*u(i1+3,i2-2,i3,m)\
 + sc( 15,i1,i2)*u(i1-3,i2-1,i3,m) + sc( 16,i1,i2)*u(i1-2,i2-1,i3,m) + sc( 17,i1,i2)*u(i1-1,i2-1,i3,m) + sc( 18,i1,i2)*u(i1+0,i2-1,i3,m) + sc( 19,i1,i2)*u(i1+1,i2-1,i3,m) + sc( 20,i1,i2)*u(i1+2,i2-1,i3,m) + sc( 21,i1,i2)*u(i1+3,i2-1,i3,m)\
 + sc( 22,i1,i2)*u(i1-3,i2+0,i3,m) + sc( 23,i1,i2)*u(i1-2,i2+0,i3,m) + sc( 24,i1,i2)*u(i1-1,i2+0,i3,m) + sc( 25,i1,i2)*u(i1+0,i2+0,i3,m) + sc( 26,i1,i2)*u(i1+1,i2+0,i3,m) + sc( 27,i1,i2)*u(i1+2,i2+0,i3,m) + sc( 28,i1,i2)*u(i1+3,i2+0,i3,m)\
 + sc( 29,i1,i2)*u(i1-3,i2+1,i3,m) + sc( 30,i1,i2)*u(i1-2,i2+1,i3,m) + sc( 31,i1,i2)*u(i1-1,i2+1,i3,m) + sc( 32,i1,i2)*u(i1+0,i2+1,i3,m) + sc( 33,i1,i2)*u(i1+1,i2+1,i3,m) + sc( 34,i1,i2)*u(i1+2,i2+1,i3,m) + sc( 35,i1,i2)*u(i1+3,i2+1,i3,m)\
 + sc( 36,i1,i2)*u(i1-3,i2+2,i3,m) + sc( 37,i1,i2)*u(i1-2,i2+2,i3,m) + sc( 38,i1,i2)*u(i1-1,i2+2,i3,m) + sc( 39,i1,i2)*u(i1+0,i2+2,i3,m) + sc( 40,i1,i2)*u(i1+1,i2+2,i3,m) + sc( 41,i1,i2)*u(i1+2,i2+2,i3,m) + sc( 42,i1,i2)*u(i1+3,i2+2,i3,m)\
 + sc( 43,i1,i2)*u(i1-3,i2+3,i3,m) + sc( 44,i1,i2)*u(i1-2,i2+3,i3,m) + sc( 45,i1,i2)*u(i1-1,i2+3,i3,m) + sc( 46,i1,i2)*u(i1+0,i2+3,i3,m) + sc( 47,i1,i2)*u(i1+1,i2+3,i3,m) + sc( 48,i1,i2)*u(i1+2,i2+3,i3,m) + sc( 49,i1,i2)*u(i1+3,i2+3,i3,m)\
 FV(m)
