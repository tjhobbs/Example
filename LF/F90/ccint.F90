!*==CCINT.spg  processed by SPAG 6.72Dc at 21:29 on  2 Nov 2016
      FUNCTION CCINT(Q2,X,Kt,Nq,Lq,Mq,Ms,Typ)
!
!  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE THE
!      SCALAR DENSITY <cbar c>
!
!  WRITTEN: T. HOBBS (FEB. 2015)
!  MODIFIED:         (SEP. 2016)
! ***********************************************************************
      IMPLICIT NONE
!*--CCINT11
      INTEGER Typ
      REAL*8 Q2 , X , Kt , kt2 , sinv , wf
      REAL*8 CCINT , Nq , Lq , Mq , Ms
      REAL*8 pi , mn
      REAL*8 kpkm_sum , kpkm2 , mq4
 
      pi = 4*DATAN(1.D0)
      mn = 0.9382720813D0     ! KEEP ALL MASSES IN GeV!
!
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
      kt2 = Kt**2
      sinv = (kt2+(1.D0-X)*Mq**2+X*Ms**2+(1.D0-X)**2*Q2/4.D0)           &
           & /(X*(1.D0-X))
 
      kpkm_sum = 2.D0*(kt2+(1.D0-X)**2*Q2/4.D0)
      kpkm2 = (kt2-(1.D0-X)**2*Q2/4.D0)**2
 
      mq4 = kpkm2*(1/X**2+1/(1.D0-X)**2+2.D0/(X*(1.D0-X)))              &
          & + kpkm_sum*((Mq/X)**2+(Ms/(1.D0-X))**2+(Mq**2+Ms**2)        &
          & /(X*(1.D0-X))) + Mq**4/X**2 + Ms**4/(1.D0-X)                &
          & **2 + 2.D0*(Mq**2*Ms**2)/(X*(1.D0-X))
 
!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
      IF ( Typ==1 ) THEN
         wf = DEXP(-sinv/Lq**2)                                       ! GAUSSIAN
      ELSEIF ( Typ==2 ) THEN
         wf = 1.D0/(1.D0+sinv/Lq**2+mq4/(4.D0*Lq**4))                 ! MONOPOLE; DIV.!
      ELSEIF ( Typ==3 ) THEN
         wf = 1.D0/(1.D0+sinv/Lq**2+mq4/(4.D0*Lq**4))**2              ! DIPOLE
      ENDIF
 
      CCINT = (Nq/(Lq**4))*(1.D0/(16.D0*pi**2))                         &
            & *(kt2+(Mq+X*mn)**2-(1.D0-X)**2*Q2/4.D0)                   &
            & *(1.D0/(X**3*(1.D0-X)))*2.D0*Kt*wf       !A DIMENSIONLESS REDEF.!!!
      END
! ***************************************************************************
! ***************************************************************************
