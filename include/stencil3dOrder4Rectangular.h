! Stencil: nd=3, orderOfAccuracy=4, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)\
                                                                                                                                                                \
                                                                                                                                                                \
                                                                 + scr( 13)*u(i1+0,i2+0,i3-2,m)                                                                \
                                                                                                                                                                \
                                                                                                                                                                \
                                                                                                                                                                \
                                                                 + scr( 33)*u(i1+0,i2-1,i3-1,m)                                                                \
                                 + scr( 37)*u(i1-1,i2+0,i3-1,m) + scr( 38)*u(i1+0,i2+0,i3-1,m) + scr( 39)*u(i1+1,i2+0,i3-1,m)                                \
                                                                 + scr( 43)*u(i1+0,i2+1,i3-1,m)                                                                \
                                                                                                                                                                \
                                                                 + scr( 53)*u(i1+0,i2-2,i3+0,m)                                                                \
                                 + scr( 57)*u(i1-1,i2-1,i3+0,m) + scr( 58)*u(i1+0,i2-1,i3+0,m) + scr( 59)*u(i1+1,i2-1,i3+0,m)                                \
 + scr( 61)*u(i1-2,i2+0,i3+0,m) + scr( 62)*u(i1-1,i2+0,i3+0,m) + scr( 63)*u(i1+0,i2+0,i3+0,m) + scr( 64)*u(i1+1,i2+0,i3+0,m) + scr( 65)*u(i1+2,i2+0,i3+0,m)\
                                 + scr( 67)*u(i1-1,i2+1,i3+0,m) + scr( 68)*u(i1+0,i2+1,i3+0,m) + scr( 69)*u(i1+1,i2+1,i3+0,m)                                \
                                                                 + scr( 73)*u(i1+0,i2+2,i3+0,m)                                                                \
                                                                                                                                                                \
                                                                 + scr( 83)*u(i1+0,i2-1,i3+1,m)                                                                \
                                 + scr( 87)*u(i1-1,i2+0,i3+1,m) + scr( 88)*u(i1+0,i2+0,i3+1,m) + scr( 89)*u(i1+1,i2+0,i3+1,m)                                \
                                                                 + scr( 93)*u(i1+0,i2+1,i3+1,m)                                                                \
                                                                                                                                                                \
                                                                                                                                                                \
                                                                                                                                                                \
                                                                 + scr(113)*u(i1+0,i2+0,i3+2,m)                                                                \
                                                                                                                                                                \
                                                                                                                                                                \
 FV(m)
