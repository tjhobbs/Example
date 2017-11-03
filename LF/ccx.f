C **********************************************************************
	 PROGRAM ccx
C
C  A SIMPLE CODE TO COMPUTE/PLOT THE x DEPENDENCE OF THE FORM FACTOR
C  INTEGRANDS; ESSENTIALLY, AT Q2=0, THE UNINTEGRATED F1 = s(x) - sbar(x),
C  ETC. 
C
C  WE INTEGRATE NUMERICALLY OVER kT2; MOREOVER,
C  THE CODE CALLS THE EXTERNAL FUNCTION F1int CONTAINED in `GEGM_5P.f'
C
C  WRITTEN: T. Hobbs (SEPT 27, 2016)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nx,ikT,nkT,typ
	PARAMETER(nx=1000)
	PARAMETER(nkT=1000)
	EXTERNAL F1int, NqbF, CCint
	REAL*8  F1int, NqbF, CCint
	REAL*8  x,xmin,xmax,xint,xpt(nx)
	REAL*8  kT,kTmin,kTmax,kTint
	REAL*8  pi,Q20,Lq,Lqb,eq2
	REAL*8  mq,mqb,mS,mSb,Nq,Nqb
	REAL*8  c(nx),cb(nx),cc(nx),dc(nx),xdc(nx),c_I,cb_I,cc_I
	REAL*8  xpc(nx), C0, DC0, DXC0, PXC0, CC0, F2c(nx)
C***********************************************************************
       pi = 4.D0*DATAN(1.D0)
       eq2 = (2.D0/3.D0)**2       ! C QUARK CHARGE SQUARED
       Q20 = 0.D0                 ! COMPUTE AT Q2=0
       typ = 1                    ! USE THE STAND. GAUSSIAN WF
       mq = 1.3D0                 ! CHARM MASS
!       mq = 0.D0                 ! CHARM MASS: ZERO! JUST A TEST...
       mqb = mq                   ! STRUCK QUARK MASSES IDENTICAL!
!__________________________________________________________________
!_ _ _ _ VALUES OF THE FITTED PARAMETERS WE PLACE HERE _ _ _ _ _ _ _ _
       Nq = 0.5D-2 * (1.D0 / 1.705730D-3)
!       Nq = 2.D0
       Lq = 3.D0     ! PARAMETER GUESSES; INITIALLY, FROM 2014 PRD
       Lqb = 3.D0
       mS = 1.865D0  ! APPROX. D MESON MASS
       mSb = 3.D0    ! FROM EBERT ET AL., (2011)

!. . . .  FOR A TEST OF THE MCMC ALGORITHMS:

!       Nq = 186.42986252200060D0
!       Lq = 6.1930970789829551D0
!       Lqb = 3.0854578213313850D0
!       mS = 5.5090061969383379D0
!       mSb = 4.2241795705035745D0

!       Nq  = 268.19D0
!       Lq  = 2.3294D0
!       Lqb = 3.6545D0
!       mS  = 4.4555D0
!       mSb = 3.5210D0

!       Nq  = 270.D0        ! FROM THE REDUCED, 2D SCANS TO TEST
!       Lq  = 1.7514D0      ! STATISTICAL CONTROLS
!       Lqb = Lq
!       mS  = 4.2318D0
!       mSb = mSb

!       mS = 0.D0     ! AGAIN, ZERO MASS LIMIT AS A TEST... 
!       mSb = 0.D0    !

! DETERMINE THE ANTI-QUARK NORM. CONSTANT FROM THE GPDs CONDITION.
       Nqb = NqbF(Nq,mq,Lq,Lqb,mS,mSb,typ)
!__________________________________________________________________

!---------------------------------------------------------------------
C***********************************************************************
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
!.................PROCEED WITH THE PDF CALCULATION.............
          C0 = 0.D0
          DC0 = 0.D0
          DXC0 = 0.D0
          PXC0 = 0.D0
          CC0 = 0.D0
! DEFINE THE BOUNDS IN x:
	  xmin = 0.D0
	     xmax = 0.999D0                        !INTEGRATE OVER x FULLY!
	     xint = (xmax-xmin)/DBLE(nx)

         DO ix = 1, nx
	      x = xmin + xint * DBLE(ix)
 
                 xpt(ix) = x

! DEFINE THE BOUNDS OF THE kT INTEGRAL --- i.e., kT \in [0, \infty):
	     c(ix) = 0.D0
	     cb(ix) = 0.D0
	     cc(ix) = 0.D0

	  kTmin = 0.D0
	     kTmax = 10.D0
	     kTint = (kTmax-kTmin)/DBLE(nkT)
	   DO ikT = 1, nkT
	      kT = kTmin + kTint * DBLE(ikT)
C________________________________________________________________________________________
!.... FOR c(x):
              c_I = F1int(Q20,x,kT,Nq,Lq,mq,mS,typ)

!.... FOR cb(x):
              cb_I = F1int(Q20,x,kT,Nqb,Lqb,mqb,mSb,typ)

!.... FOR <cc>:
              cc_I = ( CCint(Q20,x,kT,Nq,Lq,mq,mS,typ)
     &             + CCint(Q20,x,kT,Nqb,Lqb,mqb,mSb,typ) )
C________________________________________________________________________________________
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
! WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
!----- FIRST FOR THE kT INTEGRAL:
	     IF (iKT.EQ.0) THEN
	      c(ix) = c(ix) + c_I
	      cb(ix) = cb(ix) + cb_I
	      cc(ix) = cc(ix) + cc_I
	    ELSE IF (iKT/2*2.NE.iKT) THEN
	      c(ix)  = c(ix) + 4.D0*c_I
	      cb(ix) = cb(ix) + 4.D0*cb_I
	      cc(ix) = cc(ix) + 4.D0*cc_I
	    ELSE IF (iKT/2*2.EQ.iKT) THEN
	      c(ix)  = c(ix) + 2.D0*c_I
	      cb(ix) = cb(ix) + 2.D0*cb_I
	      cc(ix) = cc(ix) + 2.D0*cc_I
	     ENDIF

	   ENDDO
	    c(ix) = (kTint/3.D0) * c(ix)
	    cb(ix) = (kTint/3.D0) * cb(ix)
	    cc(ix) = (kTint/3.D0) * cc(ix)
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

          dc(ix) = c(ix) - cb(ix)
          xdc(ix) = x * dc(ix)
          xpc(ix) = x * ( c(ix) + cb(ix) )
          F2c(ix) = eq2 * x * ( c(ix) + cb(ix) )
!          F2c(ix) = ( c(ix) + cb(ix) )
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!----- AND THEN DETERMINE NUMERICAL VALUES FOR THE x MOMENTS ALSO:
             IF (ix.EQ.0) THEN
              C0 = C0 + c(ix)
              DC0 = DC0 + dc(ix)
              DXC0 = DXC0 + xdc(ix)
              PXC0 = PXC0 + xpc(ix)
              CC0 = CC0 + cc(ix)
            ELSE IF (ix/2*2.NE.ix) THEN
              C0 = C0 + 4.D0*c(ix)
              DC0 = DC0 + 4.D0*dc(ix)
              DXC0 = DXC0 + 4.D0*xdc(ix)
              PXC0 = PXC0 + 4.D0*xpc(ix)
              CC0 = CC0 + 4.D0*cc(ix)
            ELSE IF (ix/2*2.EQ.ix) THEN
              C0 = C0 + 2.D0*c(ix)
              DC0 = DC0+ 2.D0*dc(ix)
              DXC0 = DXC0+ 2.D0*xdc(ix)
              PXC0 = PXC0+ 2.D0*xpc(ix)
              CC0 = CC0 + 2.D0*cc(ix)
             ENDIF

          ENDDO
            C0 = (xint/3.D0) * C0
            DC0 = (xint/3.D0) * DC0
            DXC0 = (xint/3.D0) * DXC0
            PXC0 = (xint/3.D0) * PXC0
            CC0 = (xint/3.D0) * CC0

         PRINT*, "THE TOTAL CHARM PROB.:", C0
         PRINT*, "THE FIRST MOMENT, x*{c+cbar}:", PXC0
         PRINT*, "THE ZEROTH MOMENT, c-cbar:", DC0
         PRINT*, "THE FIRST MOMENT, x*{c-cbar}:", DXC0
         PRINT*, "THE CHARM SIGMA TERM, m_c <cc>, MeV:", mq * CC0 * 1D3
!_   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _
C...WRITE DATA TO FILE
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
        OPEN (11,FILE='DATA/test_c_cb_0.5%.dat',
!        OPEN (11,FILE='DATA/c_cb_stat-x.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          DO ix=1,nx
          WRITE (11,*) xpt(ix),  c(ix),  cb(ix), F2c(ix)
!          WRITE (11,*) xpt(ix), F2c(ix)
        ENDDO
        CLOSE (11)

        OPEN (12,FILE='DATA/test_c-cb_0.5%.dat',
!        OPEN (12,FILE='DATA/c-cb_stat-x.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          DO ix=1,nx
          WRITE (12,*) xpt(ix),  xdc(ix)
        ENDDO
        CLOSE (12)

        OPEN (13,FILE='DATA/test_sigma_0.5%.dat',
!        OPEN (13,FILE='DATA/sigma_stat-x.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          DO ix=1,nx
          WRITE (13,*) xpt(ix),  cc(ix)
        ENDDO
        CLOSE (13)
!________________________________________________________________________________________
	END
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!________________________________________________________________________________________
