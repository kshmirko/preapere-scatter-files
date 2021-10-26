module cstr

contains
	
function c_to_f_string(s) result(str)
  use iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  integer i, nchars

  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)

  str = transfer(s(1:nchars), str)

end function c_to_f_string

end module 

subroutine set_comm_name(i, s) bind(C)
	use iso_c_binding
	use cstr
	use mo_DLS, only: comm_name
	implicit none
	character(kind=c_char, len=1), intent(in) :: s(*)
	integer, intent(in)	:: i
	comm_name(i+1) = c_to_f_string(s)
end subroutine set_comm_name

subroutine set_distname_O(s) bind(C)
	use iso_c_binding
	use cstr
	use mo_DLS, only: distname_O
	implicit none
	character(kind=c_char, len=1), intent(in) :: s(*)
	
	distname_O = c_to_f_string(s)
end subroutine set_distname_O

subroutine set_distname_F(s) bind(C)
	use iso_c_binding
	use cstr
	use mo_DLS, only: distname_F
	implicit none
	character(kind=c_char, len=1), intent(in) :: s(*)
	
	distname_F = c_to_f_string(s)
end subroutine set_distname_F

subroutine set_distname_N(s) bind(C)
	use iso_c_binding
	use cstr
	use mo_DLS, only: distname_N
	implicit none
	character(kind=c_char, len=1), intent(in) :: s(*)
	
	distname_N = c_to_f_string(s)
end subroutine set_distname_N

subroutine DLS_read_config(fname, len)
use iso_c_binding
use mo_par_DLS, only : STR_LEN
implicit none
character, dimension(len) :: fname
integer 				  :: len
character(len=STR_LEN) ::bname
integer I

I=1
! Initialize bname with empty char
bname = ''
! Copy fname to bname char by char
do while(I<=len)
    bname(i:i)=fname(i)
    i=i+1
end do

write(*,*) trim(bname)
call DLS_read_input(bname)
end subroutine DLS_read_config


subroutine DLS_read_input(fullpath_input)
! 	use iso_c_binding
	use mo_par_DLS
  use mo_DLS
     
  implicit none
	
	! number of aerosol component, need development for KSD>1
  ! integer,parameter :: KSD=1    
  
	!real, dimension(KMD,KSD)   :: CM, SM, RMM     
  !real, dimension(KSD,KNpar) :: RRR, AR, xgrid
  !real, dimension(KSD)       :: AC

  real			::	rgmin, rgmax, wlmin, wlmax, pi2
  !integer		::	key_SD, ID, NMD, NSD
  integer 	::	i, j, len_of_filename
  character(len=255) :: fullpath_input
	
!	namelist /Config/ key, keyEL, keySUB, keyLS, key_org, key_fx, key_RD1, &
!	&WL, RN, RK, rgmin, rgmax, wlmin, wlmax, key_SD, ID, NMD, KN, grid, SD, &
!	&distname_O,distname_F,distname_N,KR,R,RD,KM,ANGLE,NRATN,comm_name
!  ---------------------------------------
  NSD=KSD
!
! ** READ INPUT
! 
  !fullpath_input = 'input.dat'
  !fullpath_input = ''
  !call getarg(1,fullpath_input)
  if(fullpath_input .eq. '') then
    write(*,'(a)') 'Input settings file name is missing on command line.'
    stop 'stop in DLS_read_input'
  endif

  open(10, file=trim(fullpath_input),status='old')
  read(10,*) key, keyEL, keySUB, keyLS, key_org, key_fx, key_RD1
! **   key_RD =1 - volume mixture of spheroids                  ** c
! **           2 - surface area  mixture of spheroids           ** c
  key_RD=1

  !if(keySUB.eq.0) &
	!	write(*,*) key, keyEL, keySUB, keyLS, key_org, key_fx, key_RD1
	  
	read(10,*) WL,RN,RK,rgmin,rgmax,wlmin,wlmax
	
	!if(keySUB.eq.0) write(*,*) WL, RN, RK, rgmin, rgmax, wlmin, wlmax 
	read(10,*) key_SD,ID,NMD
	
	if(key_SD.eq.0) then
			read(10,*) 
      read(10,*) KN
			if(KN.gt.KNpar) then
      	write(*,*) 'in input.dat KN=',KN,' .gt. ', &
        	&'KNpar=',KNpar,' in mo_par_DLS'
				STOP 'STOP: KN should be < or = KNpar'
			endif ! KN&KNpar
      
			!if(keySUB.eq.0) write(*,*) 'KN=',KN
      
			do i=1,KN
      	read(10,*) grid(i), SD(i)
      	!write(*,*) i,grid(i), SD(i),'  - i,grid(i), SD(i)'
      enddo ! i
	else ! key_SD=1
		read(10,*) (CM(i),SM(i),RMM(i),i=1,NMD)
		read(10,*) KN
    !if(keySUB.eq.0) write(*,*) 'KN=',KN
    do i=1,KN
      read(10,*) grid(i)
    enddo ! i
	endif ! key_SD
  xgrid(1:KN)=grid(1:KN)	
!      write(*,*) 'after grid'
	read(10,'(a)') distname_O
	read(10,'(a)') distname_F
	read(10,'(a)') distname_N
	!write(*,'(a)') 'distname_F=',TRIM(distname_F)
	!write(*,'(a)') 'distname_O=',TRIM(distname_O)
	!write(*,'(a)') 'distname_N=',TRIM(distname_N)
  read(10,*) KR
	if(KR.gt.KRpar) then
    !write(*,*) 'in input.dat KR=',KR,' .gt. '   & 
    !     ,'KRpar=',KRpar,' in mo_par_DLS'
		STOP 'STOP: KR should be < or = KRpar'
	endif ! KR&KRpar
	do i=1,KR
		read(10,*) R(i), RD(i)
	enddo ! i
  read(10,*) KM
	if(KM.gt.KMpar) then
		write(*,*) 'in input.dat KM=',KM,' .gt. ', &
			&'Kmpar=',KMpar,' in mo_par_DLS'
		STOP 'STOP: KM should be < or = KMpar'
	endif ! KM&KMpar
  do j=1,KM
  	read(10,*) ANGLE(j)
  enddo ! j
! READ common names for Kernel files
  READ(10,*) NRATN
  IF(NRATN.gt.KR1par) STOP ' in DLS_read_input_bin: NRATN.gt.KR1par !!!'
  DO i=1,NRATN
    READ(10,*) comm_name(i)
  ENDDO ! i	  
	  
  close(10)
	  
  if(key.eq.4.and.key_org.eq.1) then
  	write(*,*) ' if key=4, key_org=',key_org   &
    	,' should be 0.'
		stop 'STOP in OPTCHAR1 (matrix_optchr.f)'
	endif 
	if(key_org.eq.1.and.key_fx.eq.1) then
    write(*,*) 'STOP: key_org=1 & key_fx=1'
    write(*,*) 'If you want to use key_fx=1'   &
              ,' change key_org to 0'
    write(*,*) 'If key_org=0 is your choice,'  &
              ,' check dir_name1 in matrix_fixget'  &
              , ' and run the code again.'
		stop 
	endif

! **   key_grid1 read grid radii and scat.angles which were used** c
! **             for kernel look up tables or fixed kernels     ** c
! **          =0 - 'grid1.dat'                                  ** c
! **           1 - 'grid1.dat.fix'                              ** c
	
	if(key.eq.2) then
		key_grid1=1
	else
		key_grid1=0
	endif
!
! ** rgrid min & max calculation for fixed or NEWORG kernels
!
  if(key.eq.1.or.key_org.eq.1) then
  	pi2=2.*ACOS(-1.)
		pomin=pi2*rgmin/wlmax
		pomax=pi2*rgmax/wlmin  
	endif ! key or key_org

!
! ** CALCULATE SD if key_SD=1
!      
 if(key_SD.eq.1) then
   	!OPEN (7,FILE='LNPAR.dat',status='unknown')
   	!OPEN (9,FILE='SizeDis.dat',status='unknown')
     	call SIZEDISDN ( -KN,1,ID, NMD,&
     			&CM (1:NMD),&
     			&SM (1:NMD),&
     			&RMM(1:NMD),&
     			&xgrid(1),&
     			&xgrid(KN),&
     			&RRR(:),AR(:),AC, KNpar,1 )
     	if(keySUB.eq.0) then
     	!write(*,*) 'CM=',CM(1:NMD),' AC=',AC
       		do i=1,KN
 	  		write(*,13) i,RRR(i),AR(i)
       		enddo ! i
 	endif
 	SD(1:KN)=AR(1:KN)
    ! 	CLOSE(9)
 	!CLOSE(7)
 endif ! key_SD


13 FORMAT(i4,4E12.4)
	!write(*, nml=Config)
	return
end subroutine DLS_read_input

