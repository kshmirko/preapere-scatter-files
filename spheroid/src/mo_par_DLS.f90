      MODULE mo_par_DLS
	use iso_c_binding
!c ********************************************************** c
!c **  Parameters for matrix_...f                          ** c
!c ********************************************************** c
!c **  KN1par - number of grid radii in original kernels   ** c
!c **           and fixed kernels                          ** c
!c **  KM1par - number of scattering angles in org.kernels ** c
!c **  KR1par - number of grid aspect ratios for original  ** c
!c **                                            kernels   ** c
!c **  KREpar - number of refractive index real parts      ** c
!c **  KIMpar - number of refractive index imaginary parts ** c
!c **                                                      ** c
!c **  KNpar  - number of grid radii for optical charact.  ** c
!c **  KRpar  - number of grid aspect ratios for axial     ** c
!c **           ratio distribution                         ** c
!c **  KMpar  - number of scattering angles in fixed       ** c
!c **           kernels and for opt.characteristics        ** c
!c **  KMD    - number of modes in size distribution       ** c
!c **           (up to 2)                                  ** c
!c **  rootdir - kernel directory name                     ** c
!c ********************************************************** c
                      
!        IMPLICIT NONE

! parameters for Original or Fixed Kernels

      INTEGER(C_INT), PARAMETER ::  KN1par = 41
      INTEGER, PARAMETER ::  KR1par = 25
      INTEGER, PARAMETER ::  KM1par = 181
      INTEGER, PARAMETER ::  KREpar = 35
      INTEGER, PARAMETER ::  KIMpar = 16 

! parameters for Optical Characteristics

      INTEGER, PARAMETER ::  KNpar  = 41
      INTEGER, PARAMETER ::  KRpar  = 25
      INTEGER, PARAMETER ::  KMpar  = 181 

! parameters for Size Distribution

      INTEGER, PARAMETER ::  KMD  = 2 
!      CHARACTER(len=255), PARAMETER :: rootdir='./'
!      CHARACTER(len=255), PARAMETER :: rootdir='/Users/lapionak/MyCodes/'
	  integer, parameter	::	STR_LEN = 255
      END MODULE mo_par_DLS

      
