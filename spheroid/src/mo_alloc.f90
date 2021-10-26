MODULE alloc
! 	use iso_c_binding
	REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: UF11,UF12,UF22,UF33,UF34,UF44,UFEA
END MODULE alloc	

subroutine alloc_DLS_array(key,keyEL,key_alloc) !bind(C)
! 	use iso_c_binding
	use alloc
	use alloc1
	use mo_par_DLS,     only: KN1par,KR1par,KMpar,KREpar,KIMpar
	  	  
  implicit none
  integer, intent(in)    :: key,keyEL,key_alloc
! key_alloc = 1 - allocate arrays
! key_alloc = 2 - deallocate arrays
  integer ierr

  if(key_alloc.eq.1) then
!c*********************************************************
!c ** FOR single scattering aerosol OPTICAL properties !!!*
!c *** ALLOCATE ARRAYS to be used in subroutine USMATRIX 
!c *** (in matrix_intrpl.f)
!c
	  IF(key.lt.4) THEN 
			ALLOCATE(UFEA(2,KN1par,KIMpar,KREpar),    stat=ierr)
	    if(ierr/=0) stop 'Can not allocate UFEA array'
	    UFEA=0.
	    ALLOCATE(UF11(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	    if(ierr/=0) stop 'Can not allocate UF11 array'
	    UF11=0.
	    if(keyEL.gt.1) then
	    	ALLOCATE(UF12(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate UF12 array'
	      UF12=0.
	    endif
	    if(keyEL.gt.2) then				
	    	ALLOCATE(UF22(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate UF22 array'
	    	UF22=0.		  
			endif
	    if(keyEL.gt.3) then
	    	ALLOCATE(UF33(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate UF33 array'
	      UF33=0.
	    endif
	    if(keyEL.gt.4) then
	    	ALLOCATE(UF34(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate UF34 array'
	      UF34=0.
	    endif
	    if(keyEL.gt.5) then		  
	    	ALLOCATE(UF44(KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate UF44 array'
	    	UF44=0.
		  endif 
	  ENDIF ! key<4

	  IF(key.gt.2) THEN 
	!c
	!c *** ALLOCATE ARRAYS (alloc1.mod)
	!c
	  	ALLOCATE(UEA(KR1par,2,KN1par,KIMpar,KREpar),    stat=ierr)
	    if(ierr/=0) stop 'Can not allocate UEA array'
	    UEA=0.
		  ALLOCATE(U11(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	    if(ierr/=0) stop 'Can not allocate U11 array'
	    U11=0. 
	    if(keyEL.gt.1) then
	    	ALLOCATE(U12(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate U12 array'
	      U12=0.
	    endif
	    if(keyEL.gt.2) then
	    	ALLOCATE(U22(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate U22 array'
	      U22=0.
	    endif
	    if(keyEL.gt.3) then
				ALLOCATE(U33(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate U33 array'
	      U33=0.
	    endif
		  if(keyEL.gt.4) then
	    	ALLOCATE(U34(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate U34 array'
	      U34=0.
	    endif
		  if(keyEL.gt.5) then
	    	ALLOCATE(U44(KR1par,KMpar,KN1par,KIMpar,KREpar),stat=ierr)
	      if(ierr/=0) stop 'Can not allocate U44 array'
	      U44=0.
		  endif 
	  ENDIF ! key.gt.2            
	else ! key_alloc = 2
!c
!c *** DEALLOCATE ARRAYS (alloc.mod) in subroutine USMATRIX 
!c *** (in matrix_intrpl.f)

  	IF(key.LT.4) THEN
    	DEALLOCATE(UFEA,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate UFEA array'
      DEALLOCATE(UF11,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate UF11 array'
      if(keyEL.gt.1) then
      	DEALLOCATE(UF12,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate UF12 array'
			endif
      if(keyEL.gt.2) then
      	DEALLOCATE(UF22,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate UF22 array'
      endif
      if(keyEL.gt.3) then
      	DEALLOCATE(UF33,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate UF33 array'
      endif
      if(keyEL.gt.4) then
        DEALLOCATE(UF34,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate UF34 array'
      endif
      if(keyEL.gt.5) then
        DEALLOCATE(UF44,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate UF44 array'
			endif 
    ENDIF ! key.ne.4
    IF(key.gt.2) THEN 
!c
!c *** DEALLOCATE ARRAYS (alloc1.mod)
!c
    	DEALLOCATE(UEA,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate UEA array'
      DEALLOCATE(U11,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate U11 array'
	  	if(keyEL.gt.1) then
        DEALLOCATE(U12,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U12 array'
	  	endif
	  	if(keyEL.gt.2) then   
        DEALLOCATE(U22,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U22 array'
	  	endif
      if(keyEL.gt.3) then
				DEALLOCATE(U33,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U33 array'
	  	endif
	  	if(keyEL.gt.4) then
        DEALLOCATE(U34,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U34 array'
	  	endif
	  	if(keyEL.gt.5) then
        DEALLOCATE(U44,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U44 array'
	  	endif 
    ENDIF ! key>2    
       

  endif ! key_alloc

  return
end subroutine alloc_DLS_array

   

