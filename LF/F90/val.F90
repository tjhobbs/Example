!*==VAL.spg  processed by SPAG 6.72Dc at 00:21 on  3 Nov 2016
      PROGRAM VAL
      IMPLICIT NONE
!*--VAL4
 
      EXTERNAL AVGX
      REAL*8 AVGX
 
!        print*, avgx(0.0248D0,3.D0,3.D0,2.865D0,3.D0,1)
      PRINT * , AVGX(2.D0,3.D0,3.D0,1.865D0,3.D0,1)
      PRINT * , "THIS IS THE END."
 
      END
