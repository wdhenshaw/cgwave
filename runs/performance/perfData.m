cpuTotal  =1; % location in array with cpu total
cpuSolve  =2; % location in array with cpu for solve
cpuARC    =3; % location in array with cpu adv rectangular/curvilinear
cpuUpw    =4; % location in array with cpu upwinding/dissipation
cpuBC     =5; % location in array with cpu boundary conditions
cpuInterp =6; % location in array with cpu interpolation
square1024 = 1; % grid name enumerator
solverName{1}="MEHA";
data(cpuTotal,  2 ,1,square1024)=     9.8 ; % data(:,order,solver,grid)= value
data(cpuSolve,  2 ,1,square1024)=     8.6 ; 
data(cpuARC,    2 ,1,square1024)=     5.8 ; 
data(cpuUpw,    2 ,1,square1024)=     0.0 ; 
data(cpuBC,     2 ,1,square1024)=     1.8 ; 
data(cpuInterp, 2 ,1,square1024)=     0.0; 
data(cpuTotal,  4 ,1,square1024)=    17.4 ; % data(:,order,solver,grid)= value
data(cpuSolve,  4 ,1,square1024)=    16.3 ; 
data(cpuARC,    4 ,1,square1024)=    13.4 ; 
data(cpuUpw,    4 ,1,square1024)=     0.0 ; 
data(cpuBC,     4 ,1,square1024)=     1.9 ; 
data(cpuInterp, 4 ,1,square1024)=     0.1; 
data(cpuTotal,  6 ,1,square1024)=    37.7 ; % data(:,order,solver,grid)= value
data(cpuSolve,  6 ,1,square1024)=    36.4 ; 
data(cpuARC,    6 ,1,square1024)=    32.8 ; 
data(cpuUpw,    6 ,1,square1024)=     0.0 ; 
data(cpuBC,     6 ,1,square1024)=     2.7 ; 
data(cpuInterp, 6 ,1,square1024)=     0.1; 
nonSquare1024 = 2; % grid name enumerator
data(cpuTotal,  2 ,1,nonSquare1024)=    14.0 ; % data(:,order,solver,grid)= value
data(cpuSolve,  2 ,1,nonSquare1024)=    12.2 ; 
data(cpuARC,    2 ,1,nonSquare1024)=     8.6 ; 
data(cpuUpw,    2 ,1,nonSquare1024)=     0.0 ; 
data(cpuBC,     2 ,1,nonSquare1024)=     2.3 ; 
data(cpuInterp, 2 ,1,nonSquare1024)=     0.1; 
data(cpuTotal,  4 ,1,nonSquare1024)=    32.1 ; % data(:,order,solver,grid)= value
data(cpuSolve,  4 ,1,nonSquare1024)=    29.8 ; 
data(cpuARC,    4 ,1,nonSquare1024)=    25.2 ; 
data(cpuUpw,    4 ,1,nonSquare1024)=     0.0 ; 
data(cpuBC,     4 ,1,nonSquare1024)=     3.3 ; 
data(cpuInterp, 4 ,1,nonSquare1024)=     0.1; 
data(cpuTotal,  6 ,1,nonSquare1024)=    81.5 ; % data(:,order,solver,grid)= value
data(cpuSolve,  6 ,1,nonSquare1024)=    77.9 ; 
data(cpuARC,    6 ,1,nonSquare1024)=    61.3 ; 
data(cpuUpw,    6 ,1,nonSquare1024)=     0.0 ; 
data(cpuBC,     6 ,1,nonSquare1024)=    15.4 ; 
data(cpuInterp, 6 ,1,nonSquare1024)=     0.1; 
