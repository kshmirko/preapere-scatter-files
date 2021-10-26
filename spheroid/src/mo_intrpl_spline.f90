module mo_intrpl_spline
! 	use iso_c_binding
	 contains
      SUBROUTINE intrpl_spline(MM,XX,YY,XARG,YFIT,key,		&
                                           KS,CS)

      INTEGER,          INTENT(IN)    :: MM
      INTEGER,          INTENT(INOUT) :: key
      DOUBLE PRECISION, INTENT(IN)    :: XX(MM),YY(MM),XARG
      DOUBLE PRECISION, INTENT(OUT)   :: YFIT
      DOUBLE PRECISION, INTENT(INOUT) :: KS(MM+4),CS(MM+4)

!c      INTEGER            ::  MM,key
!c      DOUBLE PRECISION   ::  XX(MM),YY(MM),XARG
!c      DOUBLE PRECISION   ::  YFIT

      INTEGER   :: LCK,LWRK,IFAIL
      DOUBLE PRECISION ::  WRK(6*MM+16) 

        IF(key.eq.0) THEN
        LCK=MM+4
        LWRK=6*MM+16
        IFAIL=0
        CALL E01BAF(MM,XX,YY,KS,CS,LCK,WRK,LWRK,IFAIL)
        key=1
        ENDIF
 
        IFAIL=0
        CALL E02BBF(MM+4,KS,CS,XARG,YFIT,IFAIL)
	RETURN
	END SUBROUTINE intrpl_spline
!C*******************

!C*******************SPLINE**INTEREPOLATION****ROUTINES**********************************************
      SUBROUTINE E01BAF(M,X,Y,K,C,LCK,WRK,LWRK,IFAIL)
!C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
!C     MARK 11.5(F77) REVISED. (SEPT 1985.)
!C
!C     ******************************************************
!C
!C     NPL ALGORITHMS LIBRARY ROUTINE SP3INT
!C
!C     CREATED 16/5/79.                        RELEASE 00/00
!C
!C     AUTHORS ... GERALD T. ANTHONY, MAURICE G.COX
!C                 J.GEOFFREY HAYES AND MICHAEL A. SINGER.
!C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
!C     MIDDLESEX TW11 OLW, ENGLAND
!C
!C     ******************************************************
!C
!C     E01BAF.  AN ALGORITHM, WITH CHECKS, TO DETERMINE THE
!C     COEFFICIENTS IN THE B-SPLINE REPRESENTATION OF A CUBIC
!C     SPLINE WHICH INTERPOLATES (PASSES EXACTLY THROUGH) A
!C     GIVEN SET OF POINTS.
!C
!C     INPUT PARAMETERS
!C        M        THE NUMBER OF DISTINCT POINTS WHICH THE
!C                    SPLINE IS TO INTERPOLATE.
!C                    (M MUST BE AT LEAST 4.)
!C        X        ARRAY CONTAINING THE DISTINCT VALUES OF THE
!C                    INDEPENDENT VARIABLE. NB X(I) MUST BE
!C                    STRICTLY GREATER THAN X(J) WHENEVER I IS
!C                    STRICTLY GREATER THAN J.
!C        Y        ARRAY CONTAINING THE VALUES OF THE DEPENDENT
!C                    VARIABLE.
!C        LCK      THE SMALLER OF THE ACTUALLY DECLARED DIMENSIONS
!C                    OF K AND C. MUST BE AT LEAST M + 4.
!C
!C     OUTPUT PARAMETERS
!C        K        ON SUCCESSFUL EXIT, K CONTAINS THE KNOTS
!C                    SET UP BY THE ROUTINE. IF THE SPLINE IS
!C                    TO BE EVALUATED (BY NPL ROUTINE E02BEF,
!C                    FOR EXAMPLE) THE ARRAY K MUST NOT BE
!C                    ALTERED BEFORE CALLING THAT ROUTINE.
!C        C        ON SUCCESSFUL EXIT, C CONTAINS THE B-SPLINE
!C                    COEFFICIENTS OF THE INTERPOLATING SPLINE.
!C                    THESE ARE ALSO REQUIRED BY THE EVALUATING
!C                    ROUTINE E02BEF.
!C        IFAIL    FAILURE INDICATOR
!C                    0 - SUCCESSFUL TERMINATION.
!C                    1 - ONE OF THE FOLLOWING CONDITIONS HAS
!C                        BEEN VIOLATED -
!C                        M AT LEAST 4
!C                        LK AT LEAST M + 4
!C                        LWORK AT LEAST 6 * M + 16
!C                    2 - THE VALUES OF THE INDEPENDENT VARIABLE
!C                        ARE DISORDERED. IN OTHER WORDS, THE
!C                        CONDITION MENTIONED UNDER X IS NOT
!C                        SATISFIED.
!C
!C     WORKSPACE (AND ASSOCIATED DIMENSION) PARAMETERS
!C        WRK     WORKSPACE ARRAY, OF LENGTH LWRK.
!C        LWRK    ACTUAL DECLARED DIMENSION OF WRK.
!C                    MUST BE AT LEAST 6 * M + 16.
!C
!C     .. Parameters ..
      CHARACTER(LEN=6)		::	       SRNAME
      PARAMETER         (SRNAME='E01BAF')
!C     .. Scalar Arguments ..
      INTEGER           IFAIL, LCK, LWRK, M
!C     .. Array Arguments ..
      DOUBLE PRECISION  C(LCK), K(LCK), WRK(LWRK), X(M), Y(M)
!C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, SS
      INTEGER           I, IERROR, M1, M2
!C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
!C     .. External Functions ..
!      INTEGER           P01ABF		!AH original
!AH     EXTERNAL          P01ABF
!C     .. External Subroutines ..
!AH     EXTERNAL          E02BAF
!C     .. Data statements ..
      DATA              ONE/1.0D+0/
!C     .. Executable Statements ..
      IERROR = 1
!C
!C     TESTS FOR ADEQUACY OF ARRAY LENGTHS AND THAT M IS GREATER
!C     THAN 4.
!C
      IF (LWRK.LT.6*M+16 .OR. M.LT.4) GO TO 80
      IF (LCK.LT.M+4) GO TO 80
!C
!C     TESTS FOR THE CORRECT ORDERING OF THE X(I)
!C
      IERROR = 2
      DO 20 I = 2, M
         IF (X(I).LE.X(I-1)) GO TO 80
   20 CONTINUE
!C
!C     INITIALISE THE ARRAY OF KNOTS AND THE ARRAY OF WEIGHTS
!C
      WRK(1) = ONE
      WRK(2) = ONE
      WRK(3) = ONE
      WRK(4) = ONE
      IF (M.EQ.4) GO TO 60
      DO 40 I = 5, M
         K(I) = X(I-2)
         WRK(I) = ONE
   40 CONTINUE
   60 M1 = M + 1
      M2 = M1 + M
!C
!C     CALL THE SPLINE FITTING ROUTINE
!C
      IERROR = 0
      CALL E02BAF(M,M+4,X,Y,WRK,K,WRK(M1),WRK(M2),C,SS,IERROR)
!C
!C     ALL THE TESTS PERFORMED BY E02BAF ARE REDUNDANT
!C     BECAUSE OF THE ABOVE TESTS AND ASSIGNMENTS, AND SO
!C     IERROR = 0 AFTER THIS CALL.
!C
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
!C
!C     END OF E01BAF.
!C
      END SUBROUTINE E01BAF
!C**********************************************************************************************************************
      SUBROUTINE E02BAF(M,NCAP7,X,Y,W,K,WORK1,WORK2,C,SS,IFAIL)
!C     NAG COPYRIGHT 1975
!C     MARK 5 RELEASE
!C     MARK 6 REVISED  IER-84
!C     MARK 8 RE-ISSUE. IER-224 (APR 1980).
!C     MARK 9A REVISED. IER-356 (NOV 1981)
!C     MARK 11.5(F77) REVISED. (SEPT 1985.)
!C
!C     NAG LIBRARY SUBROUTINE  E02BAF
!C
!C     E02BAF  COMPUTES A WEIGHTED LEAST-SQUARES APPROXIMATION
!C     TO AN ARBITRARY SET OF DATA POINTS BY A CUBIC SPLINE
!C     WITH KNOTS PRESCRIBED BY THE USER.  CUBIC SPLINE
!C     INTERPOLATION CAN ALSO BE CARRIED OUT.
!C
!C     COX-DE BOOR METHOD FOR EVALUATING B-SPLINES WITH
!C     ADAPTATION OF GENTLEMAN*S PLANE ROTATION SCHEME FOR
!C     SOLVING OVER-DETERMINED LINEAR SYSTEMS.
!C
!C     USES NAG LIBRARY ROUTINE  P01AAF.
!C
!C     STARTED - 1973.
!C     COMPLETED - 1976.
!C     AUTHOR - MGC AND JGH.
!C
!C     REDESIGNED TO USE CLASSICAL GIVENS ROTATIONS IN
!C     ORDER TO AVOID THE OCCASIONAL UNDERFLOW (AND HENCE
!C     OVERFLOW) PROBLEMS EXPERIENCED BY GENTLEMAN*S 3-
!C     MULTIPLICATION PLANE ROTATION SCHEME
!C
!C     WORK1  AND  WORK2  ARE WORKSPACE AREAS.
!C     WORK1(R)  CONTAINS THE VALUE OF THE  R TH  DISTINCT DATA
!C     ABSCISSA AND, SUBSEQUENTLY, FOR  R = 1, 2, 3, 4,  THE
!C     VALUES OF THE NON-ZERO B-SPLINES FOR EACH SUCCESSIVE
!C     ABSCISSA VALUE.
!C     WORK2(L, J)  CONTAINS, FOR  L = 1, 2, 3, 4,  THE VALUE OF
!C     THE  J TH  ELEMENT IN THE  L TH  DIAGONAL OF THE
!C     UPPER TRIANGULAR MATRIX OF BANDWIDTH  4  IN THE
!C     TRIANGULAR SYSTEM DEFINING THE B-SPLINE COEFFICIENTS.
!C
!C     .. Parameters ..
      CHARACTER(LEN=6)		::	       SRNAME
      PARAMETER         (SRNAME='E02BAF')
!C     .. Scalar Arguments ..
      DOUBLE PRECISION  SS
      INTEGER           IFAIL, M, NCAP7
!C     .. Array Arguments ..
      DOUBLE PRECISION  C(NCAP7), K(NCAP7), W(M), WORK1(M),		&
                       WORK2(4,NCAP7), X(M), Y(M)
!C     .. Local Scalars ..
      DOUBLE PRECISION  ACOL, AROW, CCOL, COSINE, CROW, D, D4, D5, D6,		&
                       D7, D8, D9, DPRIME, E2, E3, E4, E5, K0, K1, K2,		&
                       K3, K4, K5, K6, N1, N2, N3, RELEMT, S, SIGMA,		&
                       SINE, WI, XI
      INTEGER           I, IERROR, IPLUSJ, IU, J, JOLD, JPLUSL, JREV, L,		&
                       L4, LPLUS1, LPLUSU, NCAP, NCAP3, NCAPM1, R
!C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
!C     .. External Functions ..
!      INTEGER           P01ABF		!AH original
!AH     EXTERNAL          P01ABF
!C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
!C     .. Executable Statements ..
      IERROR = 4
!C     CHECK THAT THE VALUES OF  M  AND  NCAP7  ARE REASONABLE
      IF (NCAP7.LT.8 .OR. M.LT.NCAP7-4) GO TO 420
      NCAP = NCAP7 - 7
      NCAPM1 = NCAP - 1
      NCAP3 = NCAP + 3
!C
!C     IN ORDER TO DEFINE THE FULL B-SPLINE BASIS, AUGMENT THE
!C     PRESCRIBED INTERIOR KNOTS BY KNOTS OF MULTIPLICITY FOUR
!C     AT EACH END OF THE DATA RANGE.
!C
      DO 20 J = 1, 4
         I = NCAP3 + J
         K(J) = X(1)
         K(I) = X(M)
   20 CONTINUE
!C
!C     TEST THE VALIDITY OF THE DATA.
!C
!C     CHECK THAT THE KNOTS ARE ORDERED AND ARE INTERIOR
!C     TO THE DATA INTERVAL.
!C
      IERROR = 1
      IF (K(5).LE.X(1) .OR. K(NCAP3).GE.X(M)) GO TO 420
      DO 40 J = 4, NCAP3
         IF (K(J).GT.K(J+1)) GO TO 420
   40 CONTINUE
!C
!C     CHECK THAT THE WEIGHTS ARE STRICTLY POSITIVE.
!C
      IERROR = 2
      DO 60 I = 1, M
         IF (W(I).LE.0.0D0) GO TO 420
   60 CONTINUE
!C
!C     CHECK THAT THE DATA ABSCISSAE ARE ORDERED, THEN FORM THE
!C     ARRAY  WORK1  FROM THE ARRAY  X.  THE ARRAY  WORK1  CONTAINS
!C     THE
!C     SET OF DISTINCT DATA ABSCISSAE.
!C
      IERROR = 3
      WORK1(1) = X(1)
      J = 2
      DO 80 I = 2, M
         IF (X(I).LT.WORK1(J-1)) GO TO 420
         IF (X(I).EQ.WORK1(J-1)) GO TO 80
         WORK1(J) = X(I)
         J = J + 1
   80 CONTINUE
      R = J - 1
!C
!C     CHECK THAT THERE ARE SUFFICIENT DISTINCT DATA ABSCISSAE FOR
!C     THE PRESCRIBED NUMBER OF KNOTS.
!C
      IERROR = 4
      IF (R.LT.NCAP3) GO TO 420
!C
!C     CHECK THE FIRST  S  AND THE LAST  S  SCHOENBERG-WHITNEY
!C     CONDITIONS ( S = MIN(NCAP - 1, 4) ).
!C
      IERROR = 5
      DO 100 J = 1, 4
         IF (J.GE.NCAP) GO TO 160
         I = NCAP3 - J + 1
         L = R - J + 1
         IF (WORK1(J).GE.K(J+4) .OR. K(I).GE.WORK1(L)) GO TO 420
  100 CONTINUE
!C
!C     CHECK ALL THE REMAINING SCHOENBERG-WHITNEY CONDITIONS.
!C
      IF (NCAP.LE.5) GO TO 160
      R = R - 4
      I = 4
      DO 140 J = 5, NCAPM1
         K0 = K(J+4)
         K4 = K(J)
  120    I = I + 1
         IF (WORK1(I).LE.K4) GO TO 120
         IF (I.GT.R .OR. WORK1(I).GE.K0) GO TO 420
  140 CONTINUE
!C
!C     INITIALISE A BAND TRIANGULAR SYSTEM (I.E. A
!C     MATRIX AND A RIGHT HAND SIDE) TO ZERO. THE
!C     PROCESSING OF EACH DATA POINT IN TURN RESULTS
!C     IN AN UPDATING OF THIS SYSTEM. THE SUBSEQUENT
!C     SOLUTION OF THE RESULTING BAND TRIANGULAR SYSTEM
!C     YIELDS THE COEFFICIENTS OF THE B-SPLINES.
!C
  160 DO 200 I = 1, NCAP3
         DO 180 L = 1, 4
            WORK2(L,I) = 0.0D0
  180    CONTINUE
         C(I) = 0.0D0
  200 CONTINUE
      SIGMA = 0.0D0
      J = 0
      JOLD = 0
      DO 340 I = 1, M
!C
!C        FOR THE DATA POINT  (X(I), Y(I))  DETERMINE AN INTERVAL
!C        K(J + 3) .LE. X .LT. K(J + 4)  CONTAINING  X(I).  (IN THE
!C        CASE  J + 4 .EQ. NCAP  THE SECOND EQUALITY IS RELAXED TO
!C        INCLUDE
!C        EQUALITY).
!C
         WI = W(I)
         XI = X(I)
  220    IF (XI.LT.K(J+4) .OR. J.GT.NCAPM1) GO TO 240
         J = J + 1
         GO TO 220
  240    IF (J.EQ.JOLD) GO TO 260
!C
!C        SET CERTAIN CONSTANTS RELATING TO THE INTERVAL
!C        K(J + 3) .LE. X .LE. K(J + 4).
!C
         K1 = K(J+1)
         K2 = K(J+2)
         K3 = K(J+3)
         K4 = K(J+4)
         K5 = K(J+5)
         K6 = K(J+6)
         D4 = 1.0D0/(K4-K1)
         D5 = 1.0D0/(K5-K2)
         D6 = 1.0D0/(K6-K3)
         D7 = 1.0D0/(K4-K2)
         D8 = 1.0D0/(K5-K3)
         D9 = 1.0D0/(K4-K3)
         JOLD = J
!C
!C        COMPUTE AND STORE IN  WORK1(L) (L = 1, 2, 3, 4)  THE VALUES
!C        OF
!C        THE FOUR NORMALIZED CUBIC B-SPLINES WHICH ARE NON-ZERO AT
!C        X=X(I).
!C
  260    E5 = K5 - XI
         E4 = K4 - XI
         E3 = XI - K3
         E2 = XI - K2
         N1 = WI*D9
         N2 = E3*N1*D8
         N1 = E4*N1*D7
         N3 = E3*N2*D6
         N2 = (E2*N1+E5*N2)*D5
         N1 = E4*N1*D4
         WORK1(4) = E3*N3
         WORK1(3) = E2*N2 + (K6-XI)*N3
         WORK1(2) = (XI-K1)*N1 + E5*N2
         WORK1(1) = E4*N1
         CROW = Y(I)*WI
!C
!C        ROTATE THIS ROW INTO THE BAND TRIANGULAR SYSTEM USING PLANE
!C        ROTATIONS.
!C
         DO 320 LPLUS1 = 1, 4
            L = LPLUS1 - 1
            RELEMT = WORK1(LPLUS1)
            IF (RELEMT.EQ.0.0D0) GO TO 320
            JPLUSL = J + L
            L4 = 4 - L
            D = WORK2(1,JPLUSL)
            IF (ABS(RELEMT).GE.D) DPRIME = ABS(RELEMT)		&
               *SQRT(1.0D0+(D/RELEMT)**2)
            IF (ABS(RELEMT).LT.D) DPRIME = D*SQRT(1.0D0+(RELEMT/D)**2)
            WORK2(1,JPLUSL) = DPRIME
            COSINE = D/DPRIME
            SINE = RELEMT/DPRIME
            IF (L4.LT.2) GO TO 300
            DO 280 IU = 2, L4
               LPLUSU = L + IU
               ACOL = WORK2(IU,JPLUSL)
               AROW = WORK1(LPLUSU)
               WORK2(IU,JPLUSL) = COSINE*ACOL + SINE*AROW
               WORK1(LPLUSU) = COSINE*AROW - SINE*ACOL
  280       CONTINUE
  300       CCOL = C(JPLUSL)
            C(JPLUSL) = COSINE*CCOL + SINE*CROW
            CROW = COSINE*CROW - SINE*CCOL
  320    CONTINUE
         SIGMA = SIGMA + CROW**2
  340 CONTINUE
      SS = SIGMA
!C
!C     SOLVE THE BAND TRIANGULAR SYSTEM FOR THE B-SPLINE
!C     COEFFICIENTS. IF A DIAGONAL ELEMENT IS ZERO, AND HENCE
!C     THE TRIANGULAR SYSTEM IS SINGULAR, THE IMPLICATION IS
!C     THAT THE SCHOENBERG-WHITNEY CONDITIONS ARE ONLY JUST
!C     SATISFIED. THUS IT IS APPROPRIATE TO EXIT IN THIS
!C     CASE WITH THE SAME VALUE  (IFAIL=5)  OF THE ERROR
!C     INDICATOR.
!C
      L = -1
      DO 400 JREV = 1, NCAP3
         J = NCAP3 - JREV + 1
         D = WORK2(1,J)
         IF (D.EQ.0.0D0) GO TO 420
         IF (L.LT.3) L = L + 1
         S = C(J)
         IF (L.EQ.0) GO TO 380
         DO 360 I = 1, L
            IPLUSJ = I + J
            S = S - WORK2(I+1,J)*C(IPLUSJ)
  360    CONTINUE
  380    C(J) = S/D
  400 CONTINUE
      IERROR = 0
  420 IF (IERROR) 440, 460, 440
  440 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
  460 IFAIL = 0
      RETURN
      END SUBROUTINE E02BAF
!C***************************************************************************************************************************
      SUBROUTINE E02BBF(NCAP7,K,C,X,S,IFAIL)
!C     NAG LIBRARY SUBROUTINE  E02BBF
!C
!C     E02BBF  EVALUATES A CUBIC SPLINE FROM ITS
!C     B-SPLINE REPRESENTATION.
!C
!C     DE BOOR*S METHOD OF CONVEX COMBINATIONS.
!C
!C     USES NAG LIBRARY ROUTINE  P01AAF.
!C
!C     STARTED - 1973.
!C     COMPLETED - 1976.
!C     AUTHOR - MGC AND JGH.
!C
!C     NAG COPYRIGHT 1975
!C     MARK 5 RELEASE
!C     MARK 7 REVISED IER-141 (DEC 1978)
!C     MARK 11.5(F77) REVISED. (SEPT 1985.)
!C
!C     .. Parameters ..
      CHARACTER(LEN=6)		::	       SRNAME
      PARAMETER         (SRNAME='E02BBF')
!C     .. Scalar Arguments ..
      DOUBLE PRECISION  S, X
      INTEGER           IFAIL, NCAP7
!C     .. Array Arguments ..
      DOUBLE PRECISION  C(NCAP7), K(NCAP7)
!C     .. Local Scalars ..
      DOUBLE PRECISION  C1, C2, C3, E2, E3, E4, E5, K1, K2, K3, K4, K5,		&
                       K6
      INTEGER           IERROR, J, J1, L
!C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
!C     .. External Functions ..
!      INTEGER           P01ABF
!AH     EXTERNAL          P01ABF
!C     .. Executable Statements ..
      IERROR = 0
      IF (NCAP7.GE.8) GO TO 20
      IERROR = 2
      GO TO 120
   20 IF (X.GE.K(4) .AND. X.LE.K(NCAP7-3)) GO TO 40
      IERROR = 1
      S = 0.0D0
      GO TO 120
!C
!C     DETERMINE  J  SUCH THAT  K(J + 3) .LE. X .LE. K(J + 4).
!C
   40 J1 = 0
      J = NCAP7 - 7
   60 L = (J1+J)/2
      IF (J-J1.LE.1) GO TO 100
      IF (X.GE.K(L+4)) GO TO 80
      J = L
      GO TO 60
   80 J1 = L
      GO TO 60
!C
!C     USE THE METHOD OF CONVEX COMBINATIONS TO COMPUTE  S(X).
!C
  100 K1 = K(J+1)
      K2 = K(J+2)
      K3 = K(J+3)
      K4 = K(J+4)
      K5 = K(J+5)
      K6 = K(J+6)
      E2 = X - K2
      E3 = X - K3
      E4 = K4 - X
      E5 = K5 - X
      C2 = C(J+1)
      C3 = C(J+2)
      C1 = ((X-K1)*C2+E4*C(J))/(K4-K1)
      C2 = (E2*C3+E5*C2)/(K5-K2)
      C3 = (E3*C(J+3)+(K6-X)*C3)/(K6-K3)
      C1 = (E2*C2+E4*C1)/(K4-K2)
      C2 = (E3*C3+E5*C2)/(K5-K3)
      S = (E3*C2+E4*C1)/(K4-K3)
  120 IF (IERROR) 140, 160, 140
  140 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
  160 IFAIL = 0
      RETURN
      END SUBROUTINE E02BBF
!C***************************************************************************************************************
!C***************************************************************************************************************
      FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
!C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
!C     MARK 13 REVISED. IER-621 (APR 1988).
!C     MARK 13B REVISED. IER-668 (AUG 1988).
!C
!C     P01ABF is the error-handling routine for the NAG Library.
!C
!C     P01ABF either returns the value of IERROR through the routine
!C     name (soft failure), or terminates execution of the program
!C     (hard failure). Diagnostic messages may be output.
!C
!C     If IERROR = 0 (successful exit from the calling routine),
!C     the value 0 is returned through the routine name, and no
!C     message is output
!C
!C     If IERROR is non-zero (abnormal exit from the calling routine),
!C     the action taken depends on the value of IFAIL.
!C
!C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
!C                 output)
!C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
!C     IFAIL =-13: soft failure, noisy exit but standard messages from
!C                 P01ABF are suppressed
!C     IFAIL =  0: hard failure, noisy exit
!C
!C     For compatibility with certain routines included before Mark 12
!C     P01ABF also allows an alternative specification of IFAIL in which
!C     it is regarded as a decimal integer with least significant digits
!C     cba. Then
!C
!C     a = 0: hard failure  a = 1: soft failure
!C     b = 0: silent exit   b = 1: noisy exit
!C
!C     except that hard failure now always implies a noisy exit.
!C
!C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
!C
!C     .. Scalar Arguments ..
!      INTEGER :: 	P01ABF		!AH original
implicit none
      INTEGER                 IERROR, IFAIL, NREC, P01ABF
      CHARACTER*(*)           SRNAME
!C     .. Array Arguments ..
      CHARACTER(len=*)           REC(*)
!C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER(LEN=72)		::	            MESS
!C     .. External Subroutines ..
!AH      EXTERNAL                P01ABZ, X04AAF, X04BAF
!C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
!C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
!C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.		&
            (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
!C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
!C                 Hard failure
                  CALL X04BAF(NERR,		&
                          ' ** NAG hard failure - execution terminated'		&
                             )
                  CALL P01ABZ
               ELSE
!C                 Soft failure
                  CALL X04BAF(NERR,		&
                             ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
!C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',   &
       ' =',I6)
      END FUNCTION P01ABF
!C**************************************************************************************************************
      SUBROUTINE X04AAF(I,NERR)
!C     MARK 7 RELEASE. NAG COPYRIGHT 1978
!C     MARK 7C REVISED IER-190 (MAY 1979)
!C     MARK 11.5(F77) REVISED. (SEPT 1985.)
!C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
!C     (STORED IN NERR1).
!C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
!C     VALUE SPECIFIED BY NERR.
!C
!C     .. Scalar Arguments ..
      INTEGER           I, NERR
!C     .. Local Scalars ..
      INTEGER           NERR1
!C     .. Save statement ..
      SAVE              NERR1
!C     .. Data statements ..
      DATA              NERR1/6/
!C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END SUBROUTINE X04AAF
!C*********************************************************************************************************************
      SUBROUTINE X04BAF(NOUT,REC)
!C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
!C
!C     X04BAF writes the contents of REC to the unit defined by NOUT.
!C
!C     Trailing blanks are not output, except that if REC is entirely
!C     blank, a single blank character is output.
!C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
!C     then no output occurs.
!C
!C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
!C     .. Local Scalars ..
      INTEGER           I
!C     .. Intrinsic Functions ..
      INTRINSIC         LEN
!C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
!C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
!C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
!C
99999 FORMAT (A)
      END SUBROUTINE X04BAF
!C********************************************************************************************************************************
      SUBROUTINE P01ABZ
!C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
!C
!C     Terminates execution when a hard failure occurs.
!C
!C     ******************** IMPLEMENTATION NOTE ********************
!C     The following STOP statement may be replaced by a call to an
!C     implementation-dependent routine to display a message and/or
!C     to abort the program.
!C     *************************************************************
!C     .. Executable Statements ..
      STOP
      END SUBROUTINE P01ABZ
!C****************************************************************************************************************************
end module mo_intrpl_spline
