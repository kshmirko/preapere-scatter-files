MODULE mo_DLS
 	use iso_c_binding
	use mo_par_DLS, only  : KNpar, KRpar, KMpar, KR1par, KMD 
	implicit none
	integer(C_INT), BIND(C) :: 	key, key_RD, keyEL, keySUB, keyLS, key_org, key_fx, key_grid1, &
						 	&key_RD1, KN, KM, KR, NRATN, NDP          
	real(C_FLOAT), BIND(C)	::	WL, RN, RK, pomin, pomax, xext, xabs, xsca, albedo
	real(C_FLOAT), dimension (KRpar), BIND(C)	::	R

	real(C_FLOAT), dimension (KNpar), BIND(C) :: grid,SD
	real(C_FLOAT), dimension (KRpar), BIND(C) :: RD
	real(C_FLOAT), dimension (KMpar), BIND(C) :: f11, f12, f22, f33, f34, f44
	real(C_FLOAT), dimension (KMpar), BIND(C) :: ANGLE
	real(C_FLOAT), BIND(C)	::	XBLR,XLDR

	character (len=255)                    ::	distname_O, distname_F, distname_N
	character (len=255), dimension(KR1par) ::	comm_name   
	integer(C_INT), BIND(C) ::	key_SD, ID, NMD, NSD
	! number of aerosol component, need development for KSD>1
  	integer,parameter :: KSD=1 
	real(C_FLOAT), dimension(KMD), BIND(C)   :: CM, SM, RMM     
  	real(C_FLOAT), dimension(KNpar), BIND(C) :: RRR, AR, xgrid
  	real(C_FLOAT), BIND(C)	:: AC
END MODULE mo_DLS

      
