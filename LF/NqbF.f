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
	PARAMETER(nKTI=1000)
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
