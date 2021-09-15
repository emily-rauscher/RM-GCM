      SUBROUTINE PLNK(E,T1,D)
!
!     ******************************************************
!     *  Purpose             :  Calculate Planck Function  *
!     *  Subroutines Called  :  None                       *
!     *  Input               :  WAVE, NCOUNT               *
!     *  Output              :  PLANK                      *
!     * ****************************************************
!
!  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE PLANCK FUNCTION BETWEEN
!  ZERO AND THE SPECIFIED VALUE OF LAMBDA.  THUS (USING XL AS LAMBDA)
!  WE WANT TO INTEGRATE
!  R = INTEGRAL(XL=0 TO XL=XLSPEC) ( C1*XL**-5* / (EXP(C2/XL*T)-1) )*DXL
!  SUBSTITUTING U=C2/(XL*T), THE INTEGRAL BECOMES
!  R = A CONSTANT TIMES INTEGRAL (USPEC TO INFINITY) OF
!            ( U**3 / (EXP(U) - 1) )*DU
!  THE APPROXIMATIONS SHOWN HERE ARE ON PAGE 998 OF ABRAMOWITZ AND SEGUN
!  UNDER THE HEADING OF DEBYE FUNCTIONS.  C2 IS THE PRODUCT OF PLANCK'S
!  CONSTANT AND THE SPEED OF LIGHT DIVIDED BY BOLTZMANN'S CONSTANT.
!  C2 = 14390 WHEN LAMBDA IS IN MICRONS.
!  THE FACTOR 0.15399 IS THE RECIPROCAL OF SIX TIMES
!  THE SUM OF (1/N**2) FOR ALL N FROM ONE TO INFINITY.  IT IS CHOSEN TO
!  NORMALIZE THE INTEGRAL TO A MAXIMUM VALUE OF UNITY.
!  RADIATION IN REAL UNITS IS OBTAINED BY MULTIPLYING THE INTEGRAL BY
!  THE STEFAN-BOLTZMANN CONSTANT TIMES T**4.
!
!
!   Include implicit declarations
!
      include 'rprecision.h'
!
!   Local declarations
!
      DIMENSION AM(5)
!
      D            =   0.0
      V1           =   E/T1
!
      IF (V1 .LE. 1.) THEN
         D         =  1.0 - 0.15399*V1**3 *    &
                      (1./3.-V1/8. + V1**2/60. - V1**4/5040. +    &
                      V1**6/272160. - V1**8/13305600         )
      ENDIF
!
      IF ( V1 .GT. 1. .AND. V1 .LE. 50.) THEN
         DO 100 M   =  1,5
            A       =  FLOAT(M)*V1
            AM(M)   =  0.15399 * EXP(-A)/M**4 *    &
                       (((A+3.)*A+6.)*A+6.)
 100     CONTINUE
!
         D          =  AM(1)+AM(2)+AM(3)+AM(4)+AM(5)
      ENDIF
!
      D             =  D*T1**4
!
      RETURN
      END

