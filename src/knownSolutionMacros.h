#beginMacro knownSolutionMacroComments()
  // Macros for known solutions
  // This file is included in C++ and F90
#endMacro

#beginMacro getBoxHelmholtzParameters( freq, omega, kx,ky,kz, LANG )
  #If #LANG eq "C"  
    // This macro is used in: 
    //    userDefinedKnownSolution.bC
    //    userDefinedForcing.bC
    //    solveHelmholtz.bC 
    omega = frequencyArray(freq);
    kx  = rpar[1]*twoPi*(freq*.5+1.);
    ky  = rpar[2]*twoPi*(freq*.5+1.);
    kz  = rpar[3]*twoPi*(freq*.5+1.);  
    // kx  = (rpar[1]+freq)*twoPi;
    // ky  = (rpar[2]+freq)*twoPi;
    // kz  = (rpar[3]+freq)*twoPi; 
         
  #Else

    ! This macro is used in bcOptWave.bf90
    omega = frequencyArray(freq);
    kx = kxBoxHelmholtz*(freq*.5+1.)
    ky = kyBoxHelmholtz*(freq*.5+1.)
    kz = kzBoxHelmholtz*(freq*.5+1.)
    ! kx = kxBoxHelmholtz + twoPi*freq.
    ! ky = kyBoxHelmholtz + twoPi*freq.
    ! kz = kzBoxHelmholtz + twoPi*freq.    
  #End*0.
#endMacro  