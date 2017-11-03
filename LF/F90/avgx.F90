!*==AVGX.spg  processed by SPAG 6.72Dc at 00:17 on  3 Nov 2016
! **********************************************************************
      FUNCTION AVGX(Nq,Lq,Lqb,Ms,Msb,Typ)
!
!  A FUNCTION TO COMPUTE THE AVERAGE VALUE <x>_IC OF THE NONPERT. CHARM
!   DISTRIBUTION FUNCTION. THIS IS TO BE CALLED TO COMPUTE THE PARAMETER
!   ESTIMATIONS IN THE MCMC CODE. THIS HAS BEEN CONVERTED TO F90 USING
!   A BRUTE-FORCE CONVERTER. (MUST BE CHECKED!)
!. . . NB: CONVERTED WITH, http://fortran.uk/plusfortonline.php
!
!  WRITTEN: T.J.Hobbs (Nov.2.2016)
! **********************************************************************
! VARIABLE DECLARATIONS
      IMPLICIT NONE
!*--AVGX14
      INTEGER ix , NX , ikt , NKT , Typ
      PARAMETER (NX=100)
      PARAMETER (NKT=100)
      EXTERNAL F1INT , NQBF
      REAL*8 AVGX
      REAL*8 F1INT , NQBF
      REAL*8 x , xmin , xmax , xint , xpt(NX)
      REAL*8 kt , ktmin , ktmax , ktint
      REAL*8 pi , Lq , Lqb
      REAL*8 mq , mqb , Ms , Msb , Nq , nqb
      REAL*8 c(NX) , c_i
!***********************************************************************
      pi = 4.D0*DATAN(1.D0)
      mq = 1.3D0                  ! CHARM MASS
      mqb = mq                    ! STRUCK QUARK MASSES IDENTICAL!
 
! DETERMINE THE ANTI-QUARK NORM. CONSTANT FROM THE GPDs CONDITION.
      nqb = NQBF(Nq,mq,Lq,Lqb,Ms,Msb,Typ)
!__________________________________________________________________
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!.................PROCEED WITH THE <x> CALCULATION.............
      AVGX = 0.D0
! DEFINE THE BOUNDS IN x:
      xmin = 0.D0
      xmax = 0.999D0                               !INTEGRATE OVER x FULLY!
      xint = (xmax-xmin)/DBLE(NX)
 
      DO ix = 1 , NX
         x = xmin + xint*DBLE(ix)
 
         xpt(ix) = x
 
! DEFINE THE BOUNDS OF THE kT INTEGRAL --- i.e., kT \in [0, \infty):
         c(ix) = 0.D0
 
         ktmin = 0.D0
         ktmax = 10.D0
         ktint = (ktmax-ktmin)/DBLE(NKT)
         DO ikt = 1 , NKT
            kt = ktmin + ktint*DBLE(ikt)
!________________________________________________________________________________________
!.... FOR c+cbar(x):
            c_i = x*(F1INT(0.D0,x,kt,Nq,Lq,mq,Ms,Typ)                   &
                & +F1INT(0.D0,x,kt,nqb,Lqb,mqb,Msb,Typ))
!
!. . . . TJH: WE MULT. WITH x TO EVALUATE <x>_IC
!________________________________________________________________________________________
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
! WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
!----- FIRST FOR THE kT INTEGRAL:
            IF ( ikt==0 ) THEN
               c(ix) = c(ix) + c_i
            ELSEIF ( ikt/2*2/=ikt ) THEN
               c(ix) = c(ix) + 4.D0*c_i
            ELSEIF ( ikt/2*2==ikt ) THEN
               c(ix) = c(ix) + 2.D0*c_i
            ENDIF
 
         ENDDO
         c(ix) = (ktint/3.D0)*c(ix)
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
 
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!----- AND THEN DETERMINE NUMERICAL VALUES FOR THE x MOMENTS ALSO:
         IF ( ix==0 ) THEN
            AVGX = AVGX + c(ix)
         ELSEIF ( ix/2*2/=ix ) THEN
            AVGX = AVGX + 4.D0*c(ix)
         ELSEIF ( ix/2*2==ix ) THEN
            AVGX = AVGX + 2.D0*c(ix)
         ENDIF
 
      ENDDO
      AVGX = (xint/3.D0)*AVGX
!    . . . . . . . . . . . . . . . .
!________________________________________________________________________________________
      END
!*==F1INT.spg  processed by SPAG 6.72Dc at 00:17 on  3 Nov 2016
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!________________________________________________________________________________________
 
! ***********************************************************************
      FUNCTION F1INT(Q2,X,Kt,Nq,Lq,Mq,Ms,Typ)
!
!  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE F1(Q2)
!      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED
!      OVER IN THE CALLING PROGRAM
!
!  WRITTEN: T. HOBBS (OCT. 2014)
! ***********************************************************************
      IMPLICIT NONE
!*--F1INT106
      INTEGER Typ
      REAL*8 Q2 , X , Kt , kt2 , sinv , wf
      REAL*8 F1INT , Nq , Lq , Mq , Ms
      REAL*8 pi , mn
      REAL*8 kpkm_sum , x1f , x2f , x3f
      REAL*8 aterm , bterm
 
      pi = 4*DATAN(1.D0)
      mn = 0.9382720813D0     ! PROTON MASS; IN GeV!
!
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
      kt2 = Kt**2
      sinv = (kt2+(1.D0-X)*Mq**2+X*Ms**2+(1.D0-X)**2*Q2/4.D0)           &
           & /(X*(1.D0-X))
 
      kpkm_sum = 2.D0*(kt2+(1.D0-X)**2*Q2/4.D0)
 
      x1f = (1.D0/(4.D0*Lq**4))                                         &
          & *(1.D0/X**2+1.D0/(1.D0-X)**2+2.D0/(X*(1.D0-X)))
 
      x2f = (1.D0/(4.D0*Lq**4))                                         &
          & *((Mq/X)**2+(Ms/(1.D0-X))**2+(Mq**2+Ms**2)/(X*(1.D0-X)))
 
      x3f = (1.D0/(4.D0*Lq**4))                                         &
          & *(Mq**4/X**2+Ms**4/(1.D0-X)**2+2.D0*(Mq**2*Ms**2)           &
          & /(X*(1.D0-X)))
 
      aterm = 1.D0 + sinv/Lq**2 + x1f*(kpkm_sum/2.D0)                   &
            & **2 + x2f*kpkm_sum + x3f
 
      bterm = (1.D0-X)**2*kt2*Q2*x1f
 
!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
      IF ( Typ==1 ) THEN
         wf = DEXP(-sinv/Lq**2)                                       ! GAUSSIAN
      ELSEIF ( Typ==2 ) THEN
         wf = 0.D0                                                    ! MONOPOLE...
      ELSEIF ( Typ==3 ) THEN
         wf = (1.D0/2.D0)*(2.D0*aterm-bterm)                            &
            & /aterm**3*(DSQRT(aterm/(aterm-bterm)))**3               ! DIPOLE
      ENDIF
 
      F1INT = (Nq/(Lq**4))*(1.D0/(16.D0*pi**2))                         &
            & *(kt2+(Mq+X*mn)**2-(1.D0-X)**2*Q2/4.D0)                   &
            & *(1.D0/(X**2*(1.D0-X)))*2.D0*Kt*wf       !A DIMENSIONLESS REDEF.!!!
      END
!*==NQBF.spg  processed by SPAG 6.72Dc at 00:17 on  3 Nov 2016
! ***************************************************************************
 
! **********************************************************************
      FUNCTION NQBF(Nq,Mq,Lq,Lqb,Ms,Msb,Typ)
!
!   A SIMPLE FUNCTION TO DETERMINE THE NORMALIZATION PRE-FACTOR
!   ON THE ANTI-QUARK WAVEFUNCTION GIVEN A PARTICULAR SET OF
!   PARAMETERS (I.E., Nq, Lq, Lqb, m_S, m_Sb)
!
!    THE ANALOGUE OF `Nqb_Find.f', BUT NOW INVOLVING NUMERICAL
!    INTEGRATIONS OVER kT; ALSO, WAVE-FNCTS FLAGGED BY 'typ'
!
!  WRITTEN: T. Hobbs (DEC 2, 2014)
!  MODIFIED FOR CHARM: T. Hobbs (SEPT 27, 2016)
! **********************************************************************
! VARIABLE DECLARATIONS
      IMPLICIT NONE
!*--NQBF171
      INTEGER ikti , NKTI , iw , Typ
      PARAMETER (NKTI=100)
      EXTERNAL F1INT
      REAL*8 F1INT
      REAL*8 NQBF
      REAL*8 kti , ktimax , ktimin , ktiint
      REAL*8 w , wmin , wmax , wint , is , isb , is_i , isb_i , is_kti ,&
           & isb_kti
      REAL*8 pi , mn , Lq , Lqb
      REAL*8 Mq , mqb , Ms , Msb , Nq
!***********************************************************************
      pi = 4.D0*DATAN(1.D0)
      mn = 0.9382720813D0         ! PROTON MASS, GeV
      mqb = Mq                    ! STRUCK QUARK MASSES IDENTICAL!
!---------------------------------------------------------------------
!***********************************************************************
!.  .  .  .  DETERMINE THE RELATION BETWEEN WF NORMALIZATIONS .  .  .  .
      is = 0
      isb = 0
      wmin = 0.D0
      wmax = 0.999D0                         !INTEGRATE OVER x-->w FULLY!
      wint = (wmax-wmin)/100.D0
      DO iw = 1 , 100
         w = wmin + wint*DBLE(iw)
 
         is_kti = 0.D0                    !AS WELL AS INTERNALLY OVER kTI!
         isb_kti = 0.D0
         ktimin = 0.D0
         ktimax = 10.D0
         ktiint = (ktimax-ktimin)/DBLE(NKTI)
         DO ikti = 1 , NKTI
            kti = ktimin + ktiint*DBLE(ikti)
 
            is_i = F1INT(0.D0,w,kti,1.D0,Lq,Mq,Ms,Typ)
            isb_i = F1INT(0.D0,w,kti,1.D0,Lqb,Mq,Msb,Typ)
 
            IF ( ikti==0 ) THEN
               is_kti = is_kti + is_i
               isb_kti = isb_kti + isb_i
            ELSEIF ( ikti/2*2/=ikti ) THEN
               is_kti = is_kti + 4.D0*is_i
               isb_kti = isb_kti + 4.D0*isb_i
            ELSEIF ( ikti/2*2==ikti ) THEN
               is_kti = is_kti + 2.D0*is_i
               isb_kti = isb_kti + 2.D0*isb_i
            ENDIF
 
         ENDDO
         is_kti = (ktiint/3.D0)*is_kti
         isb_kti = (ktiint/3.D0)*isb_kti
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
 
         IF ( iw==0 ) THEN
            is = is + is_kti
            isb = isb + isb_kti
         ELSEIF ( iw/2*2/=iw ) THEN
            is = is + 4.D0*is_kti
            isb = isb + 4.D0*isb_kti
         ELSEIF ( iw/2*2==iw ) THEN
            is = is + 2.D0*is_kti
            isb = isb + 2.D0*isb_kti
         ENDIF
 
      ENDDO
      is = (wint/3.D0)*is
      isb = (wint/3.D0)*isb
 
!....THESE ENABLE US TO RELATE THE NORMALIZATION CONSTANTS; THEY ARE THEN:
      NQBF = Nq*(is/isb)
!________________________________________________________________________________________
      END
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
