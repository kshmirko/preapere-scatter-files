SUBROUTINE SIZEDISDN(KN,IA,ID,NMD,CM,SM,RMM,RMIN,RMAX,RRR,AR,AC,KNpar,IDBG) 
!C*********************************************************************
!C**     Determining bi-modal LogNormal size distribution:           **
!C**         d(...)/dlnR in KN - grid radius points                  **
!C*********************************************************************
!C INPUT:
!C************
!C  KN   I(NSD) - number of radius points
!C        <0 logarithmic intervals
!C        >0 linear intervals
!C  IA   I  - defines quadrature 
!C     =0 - rectangular approximation
!C     =1 - trapezoidal approximation
!C  ID  I - dimension of d(...)/dlnR or d(...)/dR
!C         = 0 - number
!C         = 1 - radius
!C         = 2 - area
!C         = 3 - volume
!C  NSD I    - number of size distributions (up to 3)
!C  NMD I    - number of modes (up to 2)
!C  CM  R(NMD,NSD) - concentrations
!C  SM  R(NMD,NSD) - halfwidths 
!C  RM  R(NMD,NSD) - mean radii (mkm)
!C  RMIN  R(NSD) - minimum radius (mkm)
!C  RMAX  R(NSD) - maximum radius (mkm)
!C*****************************************************
!C  OUTPUT:
!C  RRR    R(NSD,KNpar) - Radii for Size Distribution
!C  AR     R(KNpar) - d()/dlnR  or d()/dR (in M)
!C  AC     R    - total concentration (M3/M3)
!C*****************************************************
!	USE mo_par_DLS, only: KNpar
  !PARAMETER (KNpar=KNpar)
!   use iso_c_binding
	implicit none
	INTEGER KNpar
  INTEGER KN,KNN,IDBG
  REAL CM(NMD),SM(NMD),RMM(NMD),RMIN,RMAX
  REAL, INTENT(OUT) :: AR(KNpar), AC, RRR(KNpar)
!cl      REAL CM0(NMD,NSD),SM0(NMD,NSD),RMM0(NMD,NSD)
  REAL CM0(NMD),RMM0(NMD)
  REAL CM1(NMD),RMM1(NMD)
  REAL CM2(NMD),RMM2(NMD)
  REAL CM3(NMD),RMM3(NMD)
	REAL CMR(NMD,0:3),RMMR(NMD,0:3),&
				&AAR(0:3),ACR(0:3),ARR(KNpar,0:3)
  REAL AA, AI, PI, RH, RR, SLOG
	integer ia, id, nmd, i, i1, im, ii, jj  
!C*****************************************************
!C  recalculation Size distribution d(...)/dlnR 
!C  parameters in different dimensions:
!C*****************************************************
 
	OPEN (7,FILE='LNPAR.dat',status='unknown')
	OPEN (9,FILE='SizeDis.dat',status='unknown')
	PI= ACOS( -1.0 )
  !DO ISD=1,NSD
  KNN=KN
  IF(KN.LT.0) KNN=-KN
  DO JJ=1,KNN
  	AR(JJ)=0.
		DO I1=0,3
      ARR(JJ,I1)=0.
		ENDDO
  ENDDO
  !ENDDO
	!print *, SM
  DO IM=1,NMD
!  	DO I=1,NSD
 	 	SM(IM)=EXP(SM(IM))
!    ENDDO
  ENDDO
  !print *, SM
	SELECT CASE (ID)
		CASE (0) 
  		DO IM=1,NMD
      	!DO I=1,NSD
        	CM0(IM)=CM(IM)
        	RMM0(IM)=RMM(IM)
        !ENDDO
      ENDDO
	  CASE (1) 
    	DO IM=1,NMD
      	!DO I=1,NSD
        	RMM0(IM)=EXP(LOG(RMM(IM))-LOG(SM(IM))*LOG(SM(IM)))
        	CM0(IM)=CM(IM)/RMM0(IM)*EXP(0.5*LOG(SM(IM))*LOG(SM(IM)))
        !ENDDO
      ENDDO
    CASE (2) 
      DO IM=1,NMD
      	!DO I=1,NSD
        	RMM0(IM)=EXP(LOG(RMM(IM))-2.0*LOG(SM(IM))*LOG(SM(IM)))
        	CM0(IM)=CM(IM)/(PI*RMM0(IM)*RMM0(IM))/  &   
		                    EXP(2.0*LOG(SM(IM))*LOG(SM(IM)))
        !ENDDO
      ENDDO
    CASE (3) 
    	DO IM=1,NMD
      	!DO I=1,NSD
       	  RMM0(IM)=EXP(LOG(RMM(IM))-3.0*LOG(SM(IM))*LOG(SM(IM)))
       	 	CM0(IM)=CM(IM)/(4.0/3.0)/   &
                      (PI*RMM0(IM)*RMM0(IM)*RMM0(IM))/&
											&EXP(9.0/2.0*LOG(SM(IM))*LOG(SM(IM)))
       	!ENDDO
      ENDDO
  END SELECT
	  
  DO IM=1,NMD
  	!DO I=1,NSD
	  	RMMR(IM,0)=RMM0(IM)
      CMR(IM,0) =CM0(IM)
    !ENDDO
  ENDDO

  I=0
  DO IM=1,NMD
  	!DO I=1,NSD
      CM1(IM)=CM0(IM)*RMM0(IM)*EXP(0.5*LOG(SM(IM))*LOG(SM(IM)))
      RMM1(IM)=EXP(LOG(RMM0(IM))+LOG(SM(IM))*LOG(SM(IM)))
      CM2(IM)=CM0(IM)*PI*RMM0(IM)*RMM0(IM)     &
	                        *EXP(2.0*LOG(SM(IM))*LOG(SM(IM)))
      RMM2(IM)=EXP(LOG(RMM0(IM))+2.0*LOG(SM(IM))*LOG(SM(IM)))
      CM3(IM)=CM0(IM)*4.0/3.0*PI*RMM0(IM)*RMM0(IM)*RMM0(IM)    &
	                        *EXP(9.0/2.0*LOG(SM(IM))*LOG(SM(IM)))
      RMM3(IM)=EXP(LOG(RMM0(IM))+3.0*LOG(SM(IM))*LOG(SM(IM))) 
      RMMR(IM,1)=RMM1(IM)      
      CMR(IM,1)=CM1(IM)
	  	RMMR(IM,2)=RMM2(IM) 
      CMR(IM,2)=CM2(IM)     
	  	RMMR(IM,3)=RMM3(IM)  
      CMR(IM,3)=CM3(IM)
	  !END DO
	END DO
!      OPEN (7,FILE='LNPAR.dat',status='unknown')
	IF (IDBG==1) THEN
	  !DO I=1,NSD 
	  	WRITE(7,*)'  halfwidth:'
	    WRITE(7,17) (LOG(SM(IM)),IM=1,NMD)
	    WRITE(7,*) ' number concentration, mean radius:'
	    WRITE(7,17) (CM0(IM),RMM0(IM),IM=1,NMD)
	    WRITE(7,*) ' radius concentration, mean radius:'
	    WRITE(7,17) (CM1(IM),RMM1(IM),IM=1,NMD)
	    WRITE(7,*) ' area concentration, mean radius:'
	    WRITE(7,17) (CM2(IM),RMM2(IM),IM=1,NMD)
	    WRITE(7,*) ' volume concentration, mean radius:'
	    WRITE(7,17) (CM3(IM),RMM3(IM),IM=1,NMD)
	  !ENDDO
	ENDIF
17    FORMAT (4E15.5) 
!C******************************************
  DO I1=0,3
  	!DO II=1,NSD
      AAR(I1)=0.0
      ACR(I1)=0.0
    !ENDDO
!cl      DO IM=1,NMD  ! delete after asking Oleg
!cl      AA(IM,)=0.0
!cl      AC(IM)=0.0
!cl      ENDDO
    !DO I=1,NSD
    	DO IM=1,NMD
      	AI=0.
      	CALL SDNORM(CMR(IM,I1),SM(IM),RMMR(IM,I1),   &
                                  RMIN,RMAX,AI)
      	AAR(I1)=AAR(I1)+AI
      	ACR(I1)=ACR(I1)+CMR(IM,I1)
      ENDDO
    !ENDDO
	ENDDO ! I1
!C***** KN<0 logarithmic intervals *******
	II = 1
	!DO II=1,NSD
		DO I=1,KNN
      AR(I)=0.0
    ENDDO
		IF (IDBG==1) THEN
    	WRITE(9,*) II,'  size distribution mode number'
		END IF
    IF(IA.EQ.1) THEN
    	IF(KN.LT.0) RH=(LOG(RMAX)-LOG(RMIN))/(KNN-1)
      IF(KN.GT.0) RH=((RMAX)-(RMIN))/(KNN-1)
    ENDIF
    IF(IA.EQ.0) THEN
      IF(KN.LT.0) RH=(LOG(RMAX)-LOG(RMIN))/(KNN-1)
      IF(KN.GT.0) RH=((RMAX)-(RMIN))/(KNN-1)
    ENDIF
    
		DO IM=1,NMD
      DO  I=1,KNN
      	IF(IA.EQ.1) THEN
      		IF(KN.LT.0) RR=EXP(LOG(RMIN)+(I-1)*RH)
      		IF(KN.GT.0) RR=RMIN+(I-1)*RH
      	ENDIF
      	IF(IA.EQ.0) THEN
      		IF(KN.LT.0) RR=EXP(LOG(RMIN)+(I-1)*RH)
      		IF(KN.GT.0) RR=RMIN+(I-1)*RH
      	ENDIF
      	RRR(I)=RR      
	 		 	DO I1=0,3
       		ARR(I,I1)=ARR(I,I1)+SLOG(CMR(IM,I1),SM(IM),RMMR(IM,I1),RR)
       	ENDDO ! I1
       	AR(I)=AR(I)+SLOG(CM(IM),SM(IM),RMM(IM),RR)
			enddo
		enddo
	!enddo
!c **  Write into 'SizeDis.dat'
	IF(IDBG==1) THEN
		!DO ISD=1,NSD
			DO I=1,KNN
	 			WRITE(9,*) RRR(I),AR(I)
	    ENDDO
	  !ENDDO ! ISD
	END IF
!c **  Write into 'LNPAR.dat'
	IF (IDBG==1) THEN
	  DO I1=0,3
			WRITE(7,*)
	    IF (I1.EQ.0) WRITE(7,*) ' ******* number distribution:'
	    IF (I1.EQ.1) WRITE(7,*) ' ******* radius distribution:'
			IF (I1.EQ.2) WRITE(7,*) ' ******* area   distribution:'
			IF (I1.EQ.3) WRITE(7,*) ' ******* volume distribution:'

	    !DO ISD=1,NSD
	      !IF (ISD.EQ.1) WRITE(7,*) 'first component:'
	      !IF (ISD.EQ.2) WRITE(7,*) 'second component:'
	      DO I=1,KNN
	      	WRITE(7,*) RRR(I),ARR(I,I1)
	      ENDDO
	      WRITE(7,*) (ACR(I1)-AAR(I1))/ACR(I1),   &
	      					'relative error in concentration due to limit: Rmin-Rmax'
	    !ENDDO ! ISD
	  ENDDO ! I1    
	END IF            

  !DO ISD=1,NSD
	  AC=ACR(3)
    DO I=1,KNN
			AR(I)=ARR(I,3)
    ENDDO
  !ENDDO ! ISD

  DO IM=1,NMD
  	!DO I=1,NSD
    	SM(IM)=LOG(SM(IM))
    !ENDDO
  ENDDO
  CLOSE(9)
  CLOSE(7)
  RETURN
END

SUBROUTINE SDNORM(C,S,RM,RMIN,RMAX,AC)
	implicit none
!C******************************************************************
!C**  Normalization of lognormal function d(...)/dlnR             ** 
!C**  for discrete interval of sizes: RMIN-RMAX                   **
!C******************************************************************
!C INPUT:
!C**********
!C  C    R  - concentration
!C  S    R  - halfwidth
!C  RM   R  - mean radius
!C  RMIN R  - minimum radius
!C  RMAX R  - maximum radius
!C******************************************************************* 
!C OUTPUT:
!C**********
!C  AC   R  - normalization constant
!C*******************************************************************
!	------------------------------------------------------------------------------------------------------
	real,intent(in)	::	C
	real,intent(in)	::	S
	real,intent(in)	::	RM
	real,intent(in)	::	RMIN
	real,intent(in)	::	RMAX
!	------------------------------------------------------------------------------
	real,intent(out)	::	AC
!	------------------------------------------------------------------------------
	real	::	AA
	real	::	ADD
	real	::	AHR
	real	::	PI 
	real	::	RR
	integer	::	I
	integer	::	KNSIM
!	------------------------------------------------------------------------------
  REAL SLOG
  KNSIM = 301
  PI    = ACOS( -1.0 )
  AHR   = (LOG(RMAX)-LOG(RMIN))/(KNSIM-1)
  AA    = 0.
  DO I=1,KNSIM
  	RR=EXP(LOG(RMIN)+(I-1)*AHR)
    ADD = SLOG(C,S,RM,RR)
    AA=AA+ADD*AHR
  ENDDO
  AC=AA
  RETURN
END
 
!cl      REAL FUNCTION SLOG(C,S,RM,RR)
!clC******************************************************************
!clC**     Lognormal function d(RR)/dlnR                            ** 
!clC******************************************************************     
!clC  C   R  - concentration
!clC  S   R  - halfwidth
!clC  RM  R  - mean radius
!clC  RR  R  - value
!clC****************************************************************** 
!cl      PI = ACOS( -1.0 )
!cl      IF(S.LE.1) WRITE(*,*) S,'TROUBLES in SLOG:S.LE.1'
!cl      SLOG=C/SQRT(2.0*PI)/LOG(S)*EXP((-0.5)*((DLOG(RR/RM)*
!cl     & DLOG(RR/RM))/(DLOG(S)*DLOG(S))))
!cl      RETURN
!cl      END
FUNCTION SLOG(C,S,RM,RR)
	implicit none
!C******************************************************************
!C**     Lognormal function d(RR)/dlnR                            ** 
!C******************************************************************     
!C  C   R  - concentration
!C  S   R  - halfwidth
!C  RM  R  - mean radius
!C  RR  R  - value
!C****************************************************************** 
	real	::	SLOG
!	------------------------------------------------------------------------------
	real, intent(in)	::	C
	real, intent(in)	::	S
	real, intent(in)	::	RM
	real, intent(in)	::	RR
!	------------------------------------------------------------------------------
	real	::	A1, A2, AA
	real	::	PI	
  real, parameter  ::  threshold = 60. !20.0*log(10.0)
!C****************************************************************** 
  PI = ACOS( -1.0 )
  IF(S.LE.1) THEN
  	WRITE(*,*) S,'TROUBLES in SLOG:S.LE.1'
    STOP 'STOP: in SLOG:S.LE.1'
  ENDIF
  A1=LOG(RR/RM)
  A2=LOG(S)
  AA = 0.5*(A1*A1)/(A2*A2)
  if(AA .gt. threshold) then
  	AA = 0.0
  else
    AA = EXP(-AA)
  end if
  SLOG=C/SQRT(2.0*PI)/A2*AA

  RETURN
END FUNCTION SLOG
 
SUBROUTINE SIZEDISDN1(KN,IA,ID,NMD,CM,SM,RMM,RMIN,RMAX,RRR,AR,AC,KNpar,IDBG) 
!C*********************************************************************
!C**     Determining bi-modal LogNormal size distribution:           **
!C**         d(...)/dlnR in KN - grid radius points                  **
!C*********************************************************************
!C INPUT:
!C************
!C  KN   I(NSD) - number of radius points
!C        <0 logarithmic intervals
!C        >0 linear intervals
!C  IA   I  - defines quadrature 
!C     =0 - rectangular approximation
!C     =1 - trapezoidal approximation
!C  ID  I - dimension of d(...)/dlnR or d(...)/dR
!C         = 0 - number
!C         = 1 - radius
!C         = 2 - area
!C         = 3 - volume
!C  NSD I    - number of size distributions (up to 3)
!C  NMD I    - number of modes (up to 2)
!C  CM  R(NMD,NSD) - concentrations
!C  SM  R(NMD,NSD) - halfwidths 
!C  RM  R(NMD,NSD) - mean radii (mkm)
!C  RMIN  R(NSD) - minimum radius (mkm)
!C  RMAX  R(NSD) - maximum radius (mkm)
!C*****************************************************
!C  OUTPUT:
!C  RRR    R(NSD,KNpar) - Radii for Size Distribution
!C  AR     R(KNpar) - d()/dlnR  or d()/dR (in M)
!C  AC     R    - total concentration (M3/M3)
!C*****************************************************
!	USE mo_par_DLS, only: KNpar
  !PARAMETER (KNpar=KNpar)
!use iso_c_binding
	implicit none
	INTEGER KNpar
  INTEGER KN,KNN,IDBG
  REAL CM(NMD),SM(NMD),RMM(NMD),RMIN,RMAX
  REAL, INTENT(INOUT) :: AR(KNpar), AC, RRR(KNpar)
!cl      REAL CM0(NMD,NSD),SM0(NMD,NSD),RMM0(NMD,NSD)
  
	integer ia, id, nmd
	
	call SIZEDISDN(KN,IA,ID,NMD,CM,SM,RMM,RMIN,RMAX,RRR,AR,AC,KNpar,IDBG)
end


SUBROUTINE SIZEDIS2(KN,NMD,CM,SM,RMM,RMIN,RMAX,RRR,AR,AC) 
!C*********************************************************************
!C**     Determining bi-modal LogNormal size distribution:           **
!C**         d(...)/dlnR in KN - grid radius points                  **
!C*********************************************************************
!C INPUT:
!C************
!C  KN   I(NSD) - number of radius points
!C        <0 logarithmic intervals
!C        >0 linear intervals
!C  IA   I  - defines quadrature 
!C     =0 - rectangular approximation
!C     =1 - trapezoidal approximation
!C  ID  I - dimension of d(...)/dlnR or d(...)/dR
!C         = 0 - number
!C         = 1 - radius
!C         = 2 - area
!C         = 3 - volume
!C  NSD I    - number of size distributions (up to 3)
!C  NMD I    - number of modes (up to 2)
!C  CM  R(NMD,NSD) - concentrations
!C  SM  R(NMD,NSD) - halfwidths 
!C  RM  R(NMD,NSD) - mean radii (mkm)
!C  RMIN  R(NSD) - minimum radius (mkm)
!C  RMAX  R(NSD) - maximum radius (mkm)
!C*****************************************************
!C  OUTPUT:
!C  RRR    R(NSD,KNpar) - Radii for Size Distribution
!C  AR     R(KNpar) - d()/dlnR  or d()/dR (in M)
!C  AC     R    - total concentration (M3/M3)
!C*****************************************************
 USE mo_par_DLS, only: KNpar
 USE mo_DLS, only: SD
	implicit none
  INTEGER, INTENT(IN) ::  KN, NMD
  REAL, INTENT(IN) :: CM(NMD),SM(NMD),RMM(NMD),RMIN,RMAX
  REAL, INTENT(OUT) :: AR(KNpar), AC, RRR(KNpar)
!cl      REAL CM0(NMD,NSD),SM0(NMD,NSD),RMM0(NMD,NSD)
  
	integer ia, id, IDBG
	IA=1
	ID=3
	IDBG=0
	call SIZEDISDN(KN,IA,ID,NMD,CM,SM,RMM,RMIN,RMAX,RRR,AR,AC,KNpar,IDBG)
	SD(1:(-KN)) = AR(1:(-KN))
	
end subroutine SIZEDIS2


SUBROUTINE SIZEDISN(KN,NMD,CM,SM,RMM,RMIN,RMAX,RRR,AR,AC) 
!C*********************************************************************
!C**     Determining bi-modal LogNormal size distribution:           **
!C**         d(...)/dlnR in KN - grid radius points                  **
!C*********************************************************************
!C INPUT:
!C************
!C  KN   I(NSD) - number of radius points
!C        <0 logarithmic intervals
!C        >0 linear intervals
!C  IA   I  - defines quadrature 
!C     =0 - rectangular approximation
!C     =1 - trapezoidal approximation
!C  ID  I - dimension of d(...)/dlnR or d(...)/dR
!C         = 0 - number
!C         = 1 - radius
!C         = 2 - area
!C         = 3 - volume
!C  NSD I    - number of size distributions (up to 3)
!C  NMD I    - number of modes (up to 2)
!C  CM  R(NMD,NSD) - concentrations
!C  SM  R(NMD,NSD) - halfwidths 
!C  RM  R(NMD,NSD) - mean radii (mkm)
!C  RMIN  R(NSD) - minimum radius (mkm)
!C  RMAX  R(NSD) - maximum radius (mkm)
!C*****************************************************
!C  OUTPUT:
!C  RRR    R(NSD,KNpar) - Radii for Size Distribution
!C  AR     R(KNpar) - d()/dlnR  or d()/dR (in M)
!C  AC     R    - total concentration (M3/M3)
!C*****************************************************
 USE mo_par_DLS, only: KNpar
 USE mo_DLS, only: SD
	implicit none
  INTEGER, INTENT(IN) ::  KN, NMD
  REAL, INTENT(IN) :: CM(NMD),SM(NMD),RMM(NMD),RMIN,RMAX
  REAL, INTENT(OUT) :: AR(KNpar), AC, RRR(KNpar)
!cl      REAL CM0(NMD,NSD),SM0(NMD,NSD),RMM0(NMD,NSD)
  
	integer ia, id, IDBG
	IA=1
	ID=0
	IDBG=0
	call SIZEDISDN(KN,IA,ID,NMD,CM,SM,RMM,RMIN,RMAX,RRR,AR,AC,KNpar,IDBG)
	SD(1:(-KN)) = AR(1:(-KN))
	
end subroutine SIZEDISN