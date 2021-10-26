module mo_intrpl_linear
	contains
!C***********************************************************************
		REAL FUNCTION  LINEAR( X, Y, M, X1 )
!C ---------------------------------------------------------------------
!C       Linear interpolation function.                                
!C ---------------------------------------------------------------------
		IMPLICIT NONE
    REAL      X( * ), Y( * )
    INTEGER  :: M, N
    REAL     :: X1
!  ---------------------------------------------------------------------	  
    IF(X(2).GT.X(1)) THEN
      
    	IF ( X1.LT.X(1) )  THEN
      	LINEAR = Y(1) - ( X1-X(1) )*( Y(2)-Y(1) )/( X(2)-X(1) )
      ELSE IF ( X1.GT.X(M) )  THEN
        LINEAR = Y(M) + ( X1-X(M) )*( Y(M)-Y(M-1) )/( X(M)-X(M-1) )
      ELSE
      	DO N = 2, M
        	IF ( X1.GE.X(N-1) .AND. X1.LE.X(N) )  LINEAR = Y(N-1) +  &
              ( X1-X(N-1) )*( Y(N)-Y(N-1) )/( X(N)-X(N-1) )
!C*****************************
        	IF(X1.EQ.X(N-1)) LINEAR = Y(N-1)
        	IF(X1.EQ.X(N))   LINEAR = Y(N)
!C*****************************
				ENDDO
      END IF
    ELSE 
      IF ( X1.GT.X(1) )  THEN
        LINEAR = Y(1) - ( X1-X(1) )*( Y(2)-Y(1) )/( X(2)-X(1) )
      ELSE IF ( X1.LT.X(M) )  THEN
        LINEAR = Y(M) + ( X1-X(M) )*( Y(M)-Y(M-1) )/( X(M)-X(M-1) )
      ELSE
        DO  N = 2, M
        	IF ( X1.LE.X(N-1) .AND. X1.GE.X(N) ) LINEAR = Y(N-1) +  &
              ( X1-X(N-1) )*( Y(N)-Y(N-1) )/( X(N)-X(N-1) )
      
!C*****************************
        	IF(X1.EQ.X(N-1)) LINEAR = Y(N-1)
        	IF(X1.EQ.X(N))   LINEAR = Y(N)
!C*****************************
				ENDDO
     	END IF

    ENDIF 
    RETURN
  END FUNCTION  LINEAR
!C**************************  END OF LINEAR  ****************************
end module mo_intrpl_linear
