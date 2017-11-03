C **********************************************************************
	 FUNCTION avgx(Nq,Lq,Lqb,mS,mSb,typ)
C
C  A FUNCTION TO COMPUTE THE AVERAGE VALUE <x>_IC OF THE NONPERT. CHARM
C   DISTRIBUTION FUNCTION. THIS IS TO BE CALLED TO COMPUTE THE PARAMETER
C   ESTIMATIONS IN THE MCMC CODE. THIS MUST BE CONVERTED TO F90; WILL
C   FIRST TRY TO DO SO USING A BRUTE-FORCE CONVERTER.
C
C  WRITTEN: T.J.Hobbs (Nov.2.2016)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nx,ikT,nkT,typ
	PARAMETER(nx=100)
	PARAMETER(nkT=100)
	EXTERNAL F1int, NqbF
	REAL*8  avgx
	REAL*8  F1int, NqbF
	REAL*8  x,xmin,xmax,xint,xpt(nx)
	REAL*8  kT,kTmin,kTmax,kTint
	REAL*8  pi,Lq,Lqb
	REAL*8  mq,mqb,mS,mSb,Nq,Nqb
	REAL*8  c(nx),c_I,xpc(nx)
C***********************************************************************
       pi = 4.D0*DATAN(1.D0)
       mq = 1.3D0                 ! CHARM MASS
       mqb = mq                   ! STRUCK QUARK MASSES IDENTICAL!

! DETERMINE THE ANTI-QUARK NORM. CONSTANT FROM THE GPDs CONDITION.
       Nqb = NqbF(Nq,mq,Lq,Lqb,mS,mSb,typ)
!__________________________________________________________________
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
!.................PROCEED WITH THE <x> CALCULATION.............
          avgx = 0.D0
! DEFINE THE BOUNDS IN x:
	  xmin = 0.D0
	     xmax = 0.999D0                        !INTEGRATE OVER x FULLY!
	     xint = (xmax-xmin)/DBLE(nx)

         DO ix = 1, nx
	      x = xmin + xint * DBLE(ix)
 
                 xpt(ix) = x

! DEFINE THE BOUNDS OF THE kT INTEGRAL --- i.e., kT \in [0, \infty):
	     c(ix) = 0.D0

	  kTmin = 0.D0
	     kTmax = 10.D0
	     kTint = (kTmax-kTmin)/DBLE(nkT)
	   DO ikT = 1, nkT
	      kT = kTmin + kTint * DBLE(ikT)
C________________________________________________________________________________________
!.... FOR c+cbar(x):
              c_I = x * ( F1int(0.D0,x,kT,Nq,Lq,mq,mS,typ)
     &            + F1int(0.D0,x,kT,Nqb,Lqb,mqb,mSb,typ) )
!
!. . . . TJH: WE MULT. WITH x TO EVALUATE <x>_IC
C________________________________________________________________________________________
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
! WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
!----- FIRST FOR THE kT INTEGRAL:
	     IF (iKT.EQ.0) THEN
	      c(ix) = c(ix) + c_I
	    ELSE IF (iKT/2*2.NE.iKT) THEN
	      c(ix)  = c(ix) + 4.D0*c_I
	    ELSE IF (iKT/2*2.EQ.iKT) THEN
	      c(ix)  = c(ix) + 2.D0*c_I
	     ENDIF

	   ENDDO
	    c(ix) = (kTint/3.D0) * c(ix)
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!----- AND THEN DETERMINE NUMERICAL VALUES FOR THE x MOMENTS ALSO:
             IF (ix.EQ.0) THEN
              avgx = avgx + c(ix)
            ELSE IF (ix/2*2.NE.ix) THEN
              avgx = avgx + 4.D0*c(ix)
            ELSE IF (ix/2*2.EQ.ix) THEN
              avgx = avgx + 2.D0*c(ix)
             ENDIF

          ENDDO
            avgx = (xint/3.D0) * avgx
!    . . . . . . . . . . . . . . . .
        RETURN
!________________________________________________________________________________________
	END
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!________________________________________________________________________________________

C ***********************************************************************
        FUNCTION F1int (Q2,x,kT,Nq,Lq,mq,mS,typ)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE F1(Q2)
C      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED
C      OVER IN THE CALLING PROGRAM
C
C  WRITTEN: T. HOBBS (OCT. 2014)
C ***********************************************************************
        IMPLICIT NONE
        INTEGER typ
        REAL*8  Q2,x,kT,kT2,SInv,WF
        REAL*8  F1int,Nq,Lq,mq,mS
        REAL*8  pi,mN
        REAL*8  kpkm_sum,X1f,X2f,X3f
        REAL*8  Aterm, Bterm

        pi = 4*DATAN(1.D0)
        mN = 0.9382720813D0   ! PROTON MASS; IN GeV!
C
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
           kT2 = kT**2
           SInv = ( kT2 + (1.D0-x)*mq**2 + x*mS**2
     &          + (1.D0-x)**2*Q2/4.D0 ) / (x * (1.D0-x) )

           kpkm_sum = 2.D0 * ( kT2 + (1.D0-x)**2*Q2/4.D0 )

           X1f = (1.D0/(4.D0*Lq**4)) * ( 1.D0/x**2 + 1.D0/(1.D0-x)**2
     &                                       + 2.D0/(x*(1.D0-x)) )

           X2f = (1.D0/(4.D0*Lq**4)) * ( (mq/x)**2 + (mS/(1.D0-x))**2
     &                                + (mq**2+mS**2)/(x*(1.D0-x)) )

           X3f = (1.D0/(4.D0*Lq**4)) * ( mq**4/x**2 + mS**4/(1.D0-x)**2
     &         + 2.D0*(mq**2*mS**2)/(x*(1.D0-x)) )

         Aterm = 1.D0 + SInv/Lq**2 + X1f*(kpkm_sum/2.D0)**2
     &         + X2f*kpkm_sum + X3f

         Bterm = (1.D0-x)**2 * kT2 * Q2 * X1f

!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.1) THEN
            WF = DEXP( -SInv / Lq**2 )                                ! GAUSSIAN
          ELSE IF (typ.EQ.2) THEN
            WF = 0.D0                                                 ! MONOPOLE...
          ELSE IF (typ.EQ.3) THEN
            WF = (1.D0/2.D0) * (2.D0*Aterm - Bterm)/Aterm**3          ! DIPOLE
     &         * ( DSQRT( Aterm / (Aterm - Bterm) ) )**3
          ENDIF

        F1int = (Nq/(Lq**4)) * ( 1.D0/(16.D0*pi**2) )  !A DIMENSIONLESS REDEF.!!!
     &             * ( kT2 + (mq + x*mN)**2 - (1.D0-x)**2*Q2/4.D0 )
     &             * ( 1.D0 / ( x**2 * (1.D0-x) ) ) * 2.D0*kT * WF
        RETURN
        END
C ***************************************************************************

C **********************************************************************
         FUNCTION NqbF(Nq,mq,Lq,Lqb,mS,mSb,typ)
C
C   A SIMPLE FUNCTION TO DETERMINE THE NORMALIZATION PRE-FACTOR
C   ON THE ANTI-QUARK WAVEFUNCTION GIVEN A PARTICULAR SET OF
C   PARAMETERS (I.E., Nq, Lq, Lqb, m_S, m_Sb)
C     
C    THE ANALOGUE OF `Nqb_Find.f', BUT NOW INVOLVING NUMERICAL
C    INTEGRATIONS OVER kT; ALSO, WAVE-FNCTS FLAGGED BY 'typ'
C
C  WRITTEN: T. Hobbs (DEC 2, 2014)
C  MODIFIED FOR CHARM: T. Hobbs (SEPT 27, 2016)
C **********************************************************************
C VARIABLE DECLARATIONS
        IMPLICIT NONE
        INTEGER iKTI,nKTI,iw,typ
        PARAMETER(nKTI=100)
        EXTERNAL F1int
        REAL*8  F1int
        REAL*8  NqbF
        REAL*8  kTI,kTImax,kTImin,kTIint
        REAL*8  w,wmin,wmax,wint, Is,Isb,Is_I,Isb_I,Is_kTI,Isb_kTI
        REAL*8  pi,mN,Lq,Lqb
        REAL*8  mq,mqb,mS,mSb,Nq
C***********************************************************************
       pi = 4.D0*DATAN(1.D0)
       mN = 0.9382720813D0        ! PROTON MASS, GeV
       mqb = mq                   ! STRUCK QUARK MASSES IDENTICAL!
!---------------------------------------------------------------------
C***********************************************************************
!.  .  .  .  DETERMINE THE RELATION BETWEEN WF NORMALIZATIONS .  .  .  .
            Is  = 0
            Isb = 0
          wmin = 0.D0
             wmax = 0.999D0                  !INTEGRATE OVER x-->w FULLY!
             wint = (wmax-wmin)/100.D0
          DO iw = 1, 100
              w = wmin + wint * DBLE(iw)

          Is_kTI = 0.D0                   !AS WELL AS INTERNALLY OVER kTI!
          Isb_kTI = 0.D0
          kTImin = 0.D0
             kTImax = 10.D0
             kTIint = (kTImax-kTImin) / DBLE(nKTI)
          DO iKTI = 1, nKTI
              kTI = kTImin + kTIint * DBLE(iKTI)

              Is_I = F1int(0.D0,w,kTI,1.D0,Lq,mq,mS,typ)
             Isb_I = F1int(0.D0,w,kTI,1.D0,Lqb,mq,mSb,typ)

             IF (iKTI.EQ.0) THEN
              Is_kTI = Is_kTI + Is_I
              Isb_kTI = Isb_kTI + Isb_I
            ELSE IF (iKTI/2*2.NE.iKTI) THEN
              Is_kTI = Is_kTI + 4.D0*Is_I
              Isb_kTI = Isb_kTI + 4.D0*Isb_I
            ELSE IF (iKTI/2*2.EQ.iKTI) THEN
              Is_kTI = Is_kTI + 2.D0*Is_I
              Isb_kTI = Isb_kTI + 2.D0*Isb_I
             ENDIF

          ENDDO
            Is_kTI = (kTIint/3.D0) * Is_kTI
            Isb_kTI = (kTIint/3.D0) * Isb_kTI
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

             IF (iw.EQ.0) THEN
              Is = Is + Is_kTI
              Isb = Isb + Isb_kTI
            ELSE IF (iw/2*2.NE.iw) THEN
              Is = Is + 4.D0*Is_kTI
              Isb = Isb + 4.D0*Isb_kTI
            ELSE IF (iw/2*2.EQ.iw) THEN
              Is = Is + 2.D0*Is_kTI
              Isb = Isb + 2.D0*Isb_kTI
             ENDIF

          ENDDO
            Is = (wint/3.D0) * Is
            Isb = (wint/3.D0) * Isb

!....THESE ENABLE US TO RELATE THE NORMALIZATION CONSTANTS; THEY ARE THEN:
          NqbF = Nq * (Is / Isb)
!________________________________________________________________________________________
        RETURN
        END
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
