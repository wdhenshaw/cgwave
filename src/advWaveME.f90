! This file automatically generated from advWaveME.bf90 with bpp.
!
! =======================================================================
! ============ Optimized advance routines for CgWave ====================
!              NEW MODIFIED EQUATION VERSIONS
!
! Feb 2022 : first version created from advWave.bf90 and hierDeriv.bf90
! =======================================================================
!
! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
! #Include "defineDiffOrder2f.h"
! #Include "defineDiffOrder4f.h"

! ! 6th-order
! #Include "defineDiffOrder6f.h"

! ! 8th order
! #Include "defineDiffOrder8f.h"

! define macros to evaluate derivatives for the 6th order method (from maple/makeGetDerivativesMacros.maple)
! #Include "../maple/defineGetSixthDerivativesMacros.h"


! ======================================================================================
!   Evaluate the TZ exact solution in 2D
! ======================================================================================

! ======================================================================================
!   Evaluate the TZ exact solution in 3D
! ======================================================================================

    


! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : beginLoops
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------




! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================


!   Include files that define macros to define the curvilinear derivative coefficients using the chain rule
! files created by chainRuleCoefficients.maple
! #Include "../maple/chainRuleCoefficientsMacro2d.h"
! #Include "../maple/chainRuleCoefficientsMacro3d.h"




! =========================================================================================
! Macro: Compute the forcing for the update of u
! =========================================================================================


! =========================================================================
! Compute errors in the derivatives and save the derivatives and errors
! =========================================================================
! #beginMacro computeErrors()
!   if( nd.eq.2 )then
!     call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex  )
!     call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey )
!     call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexx )
!     call ogDeriv(ep,0,1,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexy )
!     call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyy )
!   else
!     call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex  )
!     call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey  )
!     call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexx )
!     call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyy )
!   end if  

!   if( i1.eq.4 .and. i2.eq.6 )then  
!     write(*,'("hier: (i1,i2)=(",2i3,") uxx=",1pe12.4," uexx=",1pe12.4," err=",1pe8.2," uyy=",1pe12.4," err=",1pe8.2)') i1,i2,uxx,uexx,abs(uxx-uexx),uyy,abs(uyy-ueyy)
!   end if

!   maxErr(1) = max( maxErr(1),abs(ux-uex), abs(uy-uey) )
!   maxSol(1) = max( maxSol(1),abs(uex),abs(uey))

!   l2Err(1) = l2Err(1) + (ux-uex)**2 + (uy-uey)**2

!   maxErr(2) = max( maxErr(2),abs(uxx-uexx),abs(uyy-ueyy) ) 
!   maxErr(2) = max( maxErr(2),abs(uxy-uexy) )
!   maxSol(2) = max( maxSol(2),abs(uexx),abs(ueyy) )

!   l2Err(2) = l2Err(2) + (uxx-uexx)**2 + (uxy-uexy)**2 + (uyy-ueyy)**2

!   if( max( abs(uxx-uexx),abs(uyy-ueyy) )/maxSol(2) .gt. 1.e-3 )then
!     write(*,'("**hier: (i1,i2)=(",2i3,") uxx=",1pe12.4," uexx=",1pe12.4," err=",1pe8.2," uyy=",1pe12.4," err=",1pe8.2)') i1,i2,uxx,uexx,abs(uxx-uexx),uyy,abs(uyy-ueyy)
!   end if

!  ! Fill in the derivatives and errors:
!   d=0; 
!   ud(i1,i2,i3,d)=ux; ude(i1,i2,i3,d)=ux-uex; d=d+1;
!   ud(i1,i2,i3,d)=ux; ude(i1,i2,i3,d)=uy-uey; d=d+1;
!   if( maxDeriv.ge.2 )then
!     ud(i1,i2,i3,d)=uxx; ude(i1,i2,i3,d)=uxx-uexx; d=d+1;
!     ud(i1,i2,i3,d)=uxy; ude(i1,i2,i3,d)=uxy-uexy; d=d+1;
!     ud(i1,i2,i3,d)=uyy; ude(i1,i2,i3,d)=uyy-ueyy; d=d+1;
!   end if      

!   if( maxDeriv.ge.3 )then
!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxx  )
!       call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyy )
!       call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxy )
!       call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyy )
!     else
!       stop 3333
!     end if  
!     ! if( debug.gt.3 )then  
!     !   write(*,'(" (i1,i2)=(",2i3,") uxx=",1pe12.4," uexx=",1pe12.4," err=",1pe8.2," uyy=",1pe12.4," err=",1pe8.2)') i1,i2,uxx,uexx,abs(uxx-uexx),uyy,abs(uyy-ueyy)
!     ! end if
!     maxErr(3) = max( maxErr(3),abs(uxxx-uexxx), abs(uyyy-ueyyy) )
!     maxErr(3) = max( maxErr(3),abs(uxxy-uexxy), abs(uxyy-uexyy) )
!     maxSol(3) = max( maxSol(3),abs(uexxx),abs(uexxy),abs(uexyy),abs(ueyyy) )

!     ud(i1,i2,i3,d)=uxxx; ude(i1,i2,i3,d)=uxxx-uexxx; d=d+1;
!     ud(i1,i2,i3,d)=uxxy; ude(i1,i2,i3,d)=uxxy-uexxy; d=d+1;
!     ud(i1,i2,i3,d)=uxyy; ude(i1,i2,i3,d)=uxyy-uexyy; d=d+1;        
!     ud(i1,i2,i3,d)=uyyy; ude(i1,i2,i3,d)=uyyy-ueyyy; d=d+1; 

!     l2Err(3) = l2Err(3) + (uxxx-uexxx)**2 + (uxxy-uexxy)**2 + (uxyy-uexyy)**2 + (uyyy-ueyyy)**2

!   end if    

!   if( maxDeriv.ge.4 )then

!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,4,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxx  )
!       call ogDeriv(ep,0,3,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxy )
!       call ogDeriv(ep,0,2,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyy )
!       call ogDeriv(ep,0,1,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyyy )
!       call ogDeriv(ep,0,0,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyy )
!     else
!       stop 3333
!     end if  
!     ! if( debug.gt.3 )then  
!     !   write(*,'(" (i1,i2)=(",2i3,") uxx=",1pe12.4," uexx=",1pe12.4," err=",1pe8.2," uyy=",1pe12.4," err=",1pe8.2)') i1,i2,uxx,uexx,abs(uxx-uexx),uyy,abs(uyy-ueyy)
!     ! end if
!     maxErr(4) = max( maxErr(4),abs(uxxxx-uexxxx), abs(uyyyy-ueyyyy) )
!     maxErr(4) = max( maxErr(4),abs(uxxyy-uexxyy) )
!     maxErr(4) = max( maxErr(4),abs(uxyyy-uexyyy) )
!     maxErr(4) = max( maxErr(4),abs(uxxxy-uexxxy) )
!     maxSol(4) = max( maxSol(4),abs(uexxxx),abs(uexxyy),abs(ueyyyy) )

!     ud(i1,i2,i3,d)=uxxxx; ude(i1,i2,i3,d)=uxxxx-uexxxx; d=d+1;
!     ud(i1,i2,i3,d)=uxxxy; ude(i1,i2,i3,d)=uxxxy-uexxxy; d=d+1;
!     ud(i1,i2,i3,d)=uxxyy; ude(i1,i2,i3,d)=uxxyy-uexxyy; d=d+1;        
!     ud(i1,i2,i3,d)=uxyyy; ude(i1,i2,i3,d)=uxyyy-uexyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uyyyy; ude(i1,i2,i3,d)=uyyyy-ueyyyy; d=d+1; 

!    l2Err(4) = l2Err(4) + (uxxxx-uexxxx)**2 + (uxxxy-uexxxy)**2 + (uxxyy-uexxyy)**2 + (uxyyy-uexyyy)**2 + (uyyyy-ueyyyy)**2


!   end if 

!   if( maxDeriv.ge.5 )then
!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,5,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxx  )
!       call ogDeriv(ep,0,0,5,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyyy )
!       call ogDeriv(ep,0,4,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxy )
!       call ogDeriv(ep,0,3,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxyy )
!       call ogDeriv(ep,0,2,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyyy )
!       call ogDeriv(ep,0,1,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyyyy )
!     else
!       stop 3333
!     end if  
!     maxErr(5) = max( maxErr(5),abs(uxxxxx-uexxxxx ), abs(uyyyyy-ueyyyyy) )
!     maxErr(5) = max( maxErr(5),abs(uxxxxy-uexxxxy ), abs(uxxxyy-uexxxyy) )
!     maxErr(5) = max( maxErr(5),abs(uxxyyy-uexxyyy ), abs(uxyyyy-uexyyyy) )
!     maxSol(5) = max( maxSol(5),abs(uexxxxx),abs(uexxxxy),abs(uexxxyy),abs(uexxyyy),abs(uexyyyy),abs(ueyyyyy) )

!     ud(i1,i2,i3,d)=uxxxx; ude(i1,i2,i3,d)=uxxxx-uexxxx; d=d+1;
!     ud(i1,i2,i3,d)=uxxxy; ude(i1,i2,i3,d)=uxxxy-uexxxy; d=d+1;
!     ud(i1,i2,i3,d)=uxxyy; ude(i1,i2,i3,d)=uxxyy-uexxyy; d=d+1;        
!     ud(i1,i2,i3,d)=uxyyy; ude(i1,i2,i3,d)=uxyyy-uexyyy; d=d+1;   
!     ud(i1,i2,i3,d)=uyyyy; ude(i1,i2,i3,d)=uyyyy-ueyyyy; d=d+1;   

!     l2Err(5) = l2Err(5) + (uxxxxx-uexxxxx)**2 + (uxxxxy-uexxxxy)**2 + (uxxxyy-uexxxyy)**2 + (uxxyyy-uexxyyy)**2 + (uxyyyy-uexyyyy)**2 + (uyyyyy-ueyyyyy)**2

!   end if 

!   if( maxDeriv.ge.6)then
!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,6,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxx  )
!       call ogDeriv(ep,0,5,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxy )
!       call ogDeriv(ep,0,4,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxyy )
!       call ogDeriv(ep,0,3,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxyyy )
!       call ogDeriv(ep,0,2,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyyyy )
!       call ogDeriv(ep,0,1,5,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyyyyy )
!       call ogDeriv(ep,0,0,6,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyyyy )
!     else
!       stop 3333
!     end if  
!     if( i1.eq.5 .and. i2.eq.6 )then  
!       write(*,'(" (i1,i2)=(",2i3,") uxxxxxx=",1pe12.4," uexxxxxx=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxxx,uexxxxxx,abs(uxxxxxx-uexxxxxx)/max(1.,abs(uexxxxxx))
!       write(*,'(" (i1,i2)=(",2i3,") uxxxxxy=",1pe12.4," uexxxxxy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxxy,uexxxxxy,abs(uxxxxxy-uexxxxxy)/max(1.,abs(uexxxxxy))
!       write(*,'(" (i1,i2)=(",2i3,") uxxxxyy=",1pe12.4," uexxxxyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxyy,uexxxxyy,abs(uxxxxyy-uexxxxyy)/max(1.,abs(uexxxxyy))
!       write(*,'(" (i1,i2)=(",2i3,") uxxxyyy=",1pe12.4," uexxxyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxyyy,uexxxyyy,abs(uxxxyyy-uexxxyyy)/max(1.,abs(uexxxyyy))
!       write(*,'(" (i1,i2)=(",2i3,") uxxyyyy=",1pe12.4," uexxyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxyyyy,uexxyyyy,abs(uxxyyyy-uexxyyyy)/max(1.,abs(uexxyyyy))
!       write(*,'(" (i1,i2)=(",2i3,") uxyyyyy=",1pe12.4," uexyyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxyyyyy,uexyyyyy,abs(uxyyyyy-uexyyyyy)/max(1.,abs(uexyyyyy))
!       write(*,'(" (i1,i2)=(",2i3,") uyyyyyy=",1pe12.4," ueyyyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uyyyyyy,ueyyyyyy,abs(uyyyyyy-ueyyyyyy)/max(1.,abs(ueyyyyyy))
!     end if
!     maxErr(6) = max( maxErr(6),abs(uxxxxxx-uexxxxxx) )
!     maxErr(6) = max( maxErr(6),abs(uxxxxxy-uexxxxxy) )
!     maxErr(6) = max( maxErr(6),abs(uxxxxyy-uexxxxyy) )
!     maxErr(6) = max( maxErr(6),abs(uxxxyyy-uexxxyyy) )
!     maxErr(6) = max( maxErr(6),abs(uxxyyyy-uexxyyyy) )
!     maxErr(6) = max( maxErr(6),abs(uxyyyyy-uexyyyyy) )
!     maxErr(6) = max( maxErr(6),abs(uyyyyyy-ueyyyyyy) )
!     maxSol(6) = max( maxSol(6),abs(uexxxxxx),abs(uexxxxxy),abs(uexxxxyy),abs(uexxxyyy),abs(uexxyyyy),abs(uexyyyyy),abs(ueyyyyyy) )

!     ud(i1,i2,i3,d)=uxxxxx; ude(i1,i2,i3,d)=uxxxxx-uexxxxx; d=d+1;
!     ud(i1,i2,i3,d)=uxxxxy; ude(i1,i2,i3,d)=uxxxxy-uexxxxy; d=d+1;
!     ud(i1,i2,i3,d)=uxxxyy; ude(i1,i2,i3,d)=uxxxyy-uexxxyy; d=d+1;        
!     ud(i1,i2,i3,d)=uxxyyy; ude(i1,i2,i3,d)=uxxyyy-uexxyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uxyyyy; ude(i1,i2,i3,d)=uxyyyy-uexyyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uyyyyy; ude(i1,i2,i3,d)=uyyyyy-ueyyyyy; d=d+1; 

!     l2Err(6) = l2Err(6) + (uxxxxxx-uexxxxxx)**2 + (uxxxxxy-uexxxxxy)**2 + (uxxxxyy-uexxxxyy)**2 + (uxxxyyy-uexxxyyy)**2 + (uxxyyyy-uexxyyyy)**2 +!                           (uxyyyyy-uexyyyyy)**2 + (uyyyyyy-ueyyyyyy)**2

!   end if

!   if( maxDeriv.ge.7 )then
!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,7,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxxx )
!       call ogDeriv(ep,0,6,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxxy )
!       call ogDeriv(ep,0,5,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxyy )
!       call ogDeriv(ep,0,4,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxyyy )
!       call ogDeriv(ep,0,3,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxyyyy )
!       call ogDeriv(ep,0,2,5,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyyyyy )
!       call ogDeriv(ep,0,1,6,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyyyyyy )
!       call ogDeriv(ep,0,0,7,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyyyyy )
!     else
!       stop 3333
!     end if  
!     if( i1.eq.5 .and. i2.eq.6 )then  
!       write(*,'(" (i1,i2)=(",2i3,") uxxxxxxx=",1pe12.4," uexxxxxxx=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxxxx,uexxxxxxx,abs(uxxxxxxx-uexxxxxxx)/max(1.,abs(uexxxxxxx))
!     end if        
!     maxErr(7) = max( maxErr(7),abs(uxxxxxxx-uexxxxxxx  ) )
!     maxErr(7) = max( maxErr(7),abs(uxxxxxxy-uexxxxxxy  ) )
!     maxErr(7) = max( maxErr(7),abs(uxxxxxyy-uexxxxxyy  ) )
!     maxErr(7) = max( maxErr(7),abs(uxxxxyyy-uexxxxyyy  ) )
!     maxErr(7) = max( maxErr(7),abs(uxxxyyyy-uexxxyyyy  ) )
!     maxErr(7) = max( maxErr(7),abs(uxxyyyyy-uexxyyyyy  ) )
!     maxErr(7) = max( maxErr(7),abs(uxyyyyyy-uexyyyyyy  ) )
!     maxErr(7) = max( maxErr(7),abs(uyyyyyyy-ueyyyyyyy  ) )
!     maxSol(7) = max( maxSol(7),abs(uexxxxxxx),abs(uexxxxxxy),abs(uexxxxxyy),abs(uexxxxyyy),abs(uexxxyyyy),abs(uexxyyyyy),abs(uexyyyyyy),abs(ueyyyyyyy) )

!     ud(i1,i2,i3,d)=uxxxxxx; ude(i1,i2,i3,d)=uxxxxxx-uexxxxxx; d=d+1;
!     ud(i1,i2,i3,d)=uxxxxxy; ude(i1,i2,i3,d)=uxxxxxy-uexxxxxy; d=d+1;
!     ud(i1,i2,i3,d)=uxxxxyy; ude(i1,i2,i3,d)=uxxxxyy-uexxxxyy; d=d+1;        
!     ud(i1,i2,i3,d)=uxxxyyy; ude(i1,i2,i3,d)=uxxxyyy-uexxxyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uxxyyyy; ude(i1,i2,i3,d)=uxxyyyy-uexxyyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uxyyyyy; ude(i1,i2,i3,d)=uxyyyyy-uexyyyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uyyyyyy; ude(i1,i2,i3,d)=uyyyyyy-ueyyyyyy; d=d+1; 

!     l2Err(7) = l2Err(7) + (uxxxxxxx-uexxxxxxx)**2 + (uxxxxxxy-uexxxxxxy)**2 + (uxxxxxyy-uexxxxxyy)**2 + (uxxxxyyy-uexxxxyyy)**2 + (uxxxyyyy-uexxxyyyy)**2+!                           (uxxyyyyy-uexxyyyyy)**2 + (uxyyyyyy-uexyyyyyy)**2 + (uyyyyyyy-ueyyyyyyy)**2    

!   end if

!   if( maxDeriv.ge.8 )then
!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,8,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxxxx )
!       call ogDeriv(ep,0,7,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxxxy )
!       call ogDeriv(ep,0,6,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxxyy )
!       call ogDeriv(ep,0,5,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxxyyy )
!       call ogDeriv(ep,0,4,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxxyyyy )
!       call ogDeriv(ep,0,3,5,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxyyyyy )
!       call ogDeriv(ep,0,2,6,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyyyyyy )
!       call ogDeriv(ep,0,1,7,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyyyyyyy )
!       call ogDeriv(ep,0,0,8,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyyyyyy )
!     else
!       stop 3333
!     end if  
!     if( i1.eq.5 .and. i2.eq.6 )then  
!       write(*,'("D8 (i1,i2)=(",2i3,") uxxxxxxxx=",1pe12.4," uexxxxxxxx=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxxxxx,uexxxxxxxx,abs(uxxxxxxxx-uexxxxxxxx)/max(1.,abs(uexxxxxxxx))
!       write(*,'("D8 (i1,i2)=(",2i3,") uxxxxxxxy=",1pe12.4," uexxxxxxxy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxxxxy,uexxxxxxxy,abs(uxxxxxxxy-uexxxxxxxy)/max(1.,abs(uexxxxxxxy))
!       write(*,'("D8 (i1,i2)=(",2i3,") uxxxxxxyy=",1pe12.4," uexxxxxxyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxxxyy,uexxxxxxyy,abs(uxxxxxxyy-uexxxxxxyy)/max(1.,abs(uexxxxxxyy))
!       write(*,'("D8 (i1,i2)=(",2i3,") uxxxxxyyy=",1pe12.4," uexxxxxyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxxyyy,uexxxxxyyy,abs(uxxxxxyyy-uexxxxxyyy)/max(1.,abs(uexxxxxyyy))
!       write(*,'("D8 (i1,i2)=(",2i3,") uxxxxyyyy=",1pe12.4," uexxxxyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxxyyyy,uexxxxyyyy,abs(uxxxxyyyy-uexxxxyyyy)/max(1.,abs(uexxxxyyyy))
!       write(*,'("D8 (i1,i2)=(",2i3,") uxxxyyyyy=",1pe12.4," uexxxyyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxxyyyyy,uexxxyyyyy,abs(uxxxyyyyy-uexxxyyyyy)/max(1.,abs(uexxxyyyyy))
!       write(*,'("D8 (i1,i2)=(",2i3,") uxxyyyyyy=",1pe12.4," uexxyyyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxxyyyyyy,uexxyyyyyy,abs(uxxyyyyyy-uexxyyyyyy)/max(1.,abs(uexxyyyyyy))
!       write(*,'("D8 (i1,i2)=(",2i3,") uxyyyyyyy=",1pe12.4," uexyyyyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uxyyyyyyy,uexyyyyyyy,abs(uxyyyyyyy-uexyyyyyyy)/max(1.,abs(uexyyyyyyy))
!       write(*,'("D8 (i1,i2)=(",2i3,") uyyyyyyyy=",1pe12.4," ueyyyyyyyy=",1pe12.4," rel-err=",1pe8.2)') i1,i2,uyyyyyyyy,ueyyyyyyyy,abs(uyyyyyyyy-ueyyyyyyyy)/max(1.,abs(ueyyyyyyyy))
!     end if

!     maxErr(8) = max( maxErr(8),abs(uxxxxxxxx-uexxxxxxxx  ) )
!     maxErr(8) = max( maxErr(8),abs(uxxxxxxxy-uexxxxxxxy  ) )
!     maxErr(8) = max( maxErr(8),abs(uxxxxxxyy-uexxxxxxyy  ) )
!     maxErr(8) = max( maxErr(8),abs(uxxxxxyyy-uexxxxxyyy  ) )
!     maxErr(8) = max( maxErr(8),abs(uxxxxyyyy-uexxxxyyyy  ) )
!     maxErr(8) = max( maxErr(8),abs(uxxxyyyyy-uexxxyyyyy  ) )
!     maxErr(8) = max( maxErr(8),abs(uxxyyyyyy-uexxyyyyyy  ) )
!     maxErr(8) = max( maxErr(8),abs(uxyyyyyyy-uexyyyyyyy  ) )
!     maxErr(8) = max( maxErr(8),abs(uyyyyyyyy-ueyyyyyyyy  ) )
!     maxSol(8) = max( maxSol(8),abs(uexxxxxxxx),abs(uexxxxxxxy),abs(uexxxxxxyy),abs(uexxxxxyyy),abs(uexxxxyyyy),abs(uexxxyyyyy),abs(uexxyyyyyy),abs(uexyyyyyyy),abs(ueyyyyyyyy) )

!     ud(i1,i2,i3,d)=uxxxxxxx; ude(i1,i2,i3,d)=uxxxxxxx-uexxxxxxx; d=d+1;
!     ud(i1,i2,i3,d)=uxxxxxxy; ude(i1,i2,i3,d)=uxxxxxxy-uexxxxxxy; d=d+1;
!     ud(i1,i2,i3,d)=uxxxxxyy; ude(i1,i2,i3,d)=uxxxxxyy-uexxxxxyy; d=d+1;        
!     ud(i1,i2,i3,d)=uxxxxyyy; ude(i1,i2,i3,d)=uxxxxyyy-uexxxxyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uxxxyyyy; ude(i1,i2,i3,d)=uxxxyyyy-uexxxyyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uxxyyyyy; ude(i1,i2,i3,d)=uxxyyyyy-uexxyyyyy; d=d+1; 
!     ud(i1,i2,i3,d)=uxyyyyyy; ude(i1,i2,i3,d)=uxyyyyyy-uexyyyyyy; d=d+1;         
!     ud(i1,i2,i3,d)=uyyyyyyy; ude(i1,i2,i3,d)=uyyyyyyy-ueyyyyyyy; d=d+1;   

!     l2Err(8) = l2Err(8) + (uxxxxxxxx-uexxxxxxxx)**2 + (uxxxxxxxy-uexxxxxxxy)**2 + (uxxxxxxyy-uexxxxxxyy)**2 + (uxxxxxyyy-uexxxxxyyy)**2 + (uxxxxyyyy-uexxxxyyyy)**2+!                           (uxxxyyyyy-uexxxyyyyy)**2 + (uxxyyyyyy-uexxyyyyyy)**2 + (uxyyyyyyy-uexyyyyyyy)**2 + (uyyyyyyyy-ueyyyyyyyy)**2             
!   end if  

!   count = count + 1                  
! #endMacro







! --------- include macros ORDER=2 --------
! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=2, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=2, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=2, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=3, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=3, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================



! --------- include macros ORDER=4 --------
! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=2, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=2, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=2, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=3, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=3, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================


! --------- include macros ORDER=6 --------
! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! --------- include macros ORDER=8 --------
! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=2, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=2, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=2, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================

! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
! -----------------------------------------------------------------------
! Macro to assign neighbours of metric terms when comuputing derivatives:
!   rxi1g(-2), rxi1g(-1), rxi1g(1), rxi1g(2) (2 ghost)
!   rxi2g(-2), rxi2g(-1), rxi2g(1), rxi2g(2) (2 ghost)
!   rxi3g(-2), rxi3g(-1), rxi3g(1), rxi3g(2) (2 ghost, 3d)
! Extrapolate when necessary
! -----------------------------------------------------------------------

! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================






! =========================================================================================
!
! Macro: Advance the wave equation, EXPLICIT TIME-STEPPING
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   ORDERINTIME : 2 or 4 
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================



! Argument list

! **********************************************************************************
! Macro ADV_WAVE:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************



! ==============================================
! ============== CONSTRUCT A FILE ==============
! ==============================================









            subroutine advWaveME( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
!======================================================================
!   Advance a time step for Maxwells eqution
!     OPTIMIZED version for rectangular grids.
! nd : number of space dimensions
!
! ipar(0)  = option : option=0 - Maxwell+Artificial diffusion
!                           =1 - AD only
!======================================================================
            implicit none
            integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

            real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
            real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
            real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme

            real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
            real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

            integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
            integer bc(0:1,0:2),ierr

            integer ipar(0:*)
            real rpar(0:*)

            real frequencyArray(0:*)
            
!     ---- local variables -----
            integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime
            integer addForcing,orderOfDissipation,option
            integer useWhereMask,solveForE,solveForH,grid
            integer ex,ey,ez, hx,hy,hz

            integer rectangular,curvilinear
            parameter( rectangular=0, curvilinear=1 )
!...........end   statement functions


      ! write(*,*) 'Inside advWave...'

            gridType           =ipar(2)
            orderOfAccuracy    =ipar(3)

            if( orderOfAccuracy.eq.2 )then

        ! if( nd.eq.2 .and. gridType.eq.rectangular ) then
        !   call advWaveME2dOrder2r( ARGLIST() )
        ! else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
        !   call advWaveME2dOrder2c( ARGLIST() )
        ! else if( nd.eq.3 .and. gridType.eq.rectangular ) then
        !   call advWaveME3dOrder2r( ARGLIST() )
        ! else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
        !   call advWaveME3dOrder2c( ARGLIST() )
        ! else
        !   stop 2271
        ! end if
                stop 2222

            else if( orderOfAccuracy.eq.4 ) then
                if( nd.eq.2 .and. gridType.eq.rectangular )then
                    call advWaveME2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
       !  else if(nd.eq.2 .and. gridType.eq.curvilinear )then
       !    call advWaveME2dOrder4c( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       !    call advWaveME3dOrder4r( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       !    call advWaveME3dOrder4c( ARGLIST() )
       ! else
       !   stop 8843
              end if
              stop 4444

            else if( orderOfAccuracy.eq.6 ) then

                if( nd.eq.2 .and. gridType.eq.rectangular )then
                    call advWaveME2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
        ! else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        !   !call advWaveME2dOrder6c( ARGLIST() )
                else if(  nd.eq.3 .and. gridType.eq.rectangular )then
                    call advWaveME3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
        ! else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        !   !call advWaveME3dOrder6c( ARGLIST() )
              else
                  stop 8843
              end if

            else if( orderOfAccuracy.eq.8 ) then

                if( nd.eq.2 .and. gridType.eq.rectangular )then
                    call advWaveME2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
        !  else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        !    call advWaveME2dOrder8c( ARGLIST() )
                else if(  nd.eq.3 .and. gridType.eq.rectangular )then
                    call advWaveME3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
        !  else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        !    call advWaveME3dOrder8c( ARGLIST() )
                else
                    stop 8843
                end if

            else
                write(*,'(" advWaveME:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy

                stop 11122
            end if

            return
            end








