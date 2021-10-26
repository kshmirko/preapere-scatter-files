!      SUBROUTINE OPTCHAR(key,key_RD,keyEL,keySUB,keyLS
!     &                  ,key_org,key_fx,key_grid1
!     &                  ,WL,RN,RK,KN,grid,SD
!     &                  ,KR,R,RD,KM,ANGLE,xext,xabs,xsca,albedo
!     &                  ,f11,f12,f22,f33,f34,f44,pomin,pomax
!     &                  ,distname_O,distname_F,distname_N,NDP
!     &                  ,key_RD1)

      SUBROUTINE OPTCHAR(NDP) !bind(C)

!c ** matrix_optchar_LS.f previous file name 30/06/2011
!c **************************************************************** c
!c **   08/26/02                                                 ** c
!c **   Subroutine calculates optical characteristics for given  ** c
!c **   size distribution, refractive index, axis ratio          ** c
!c **   distribution and wavelength.                             ** c
!c **                                                            ** c
!c **   In case RN or RK or R() is out of corresponding ranges   ** c
!c **   subroutine changes RN or RK or R() for edge value and    ** c
!c **   gives optical characteristics for new values             ** c
!c **************************************************************** c
!c **                                                            ** c
!c ** INPUT:                                                     ** c
!c **                                                            ** c
!c **   key  = 1 - create fixed kernels (for fixed axis          ** c
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **          2 - read fixed kernels from 'Rke...fix...'  files ** c
!c **          3 - create fixed kernels but don't save them      ** c
!c **          4 - don't create fixed kernels, calculate         ** c
!c **              opt.characteristics from original kernels     ** c
!c **   key_RD =1 - volume mixture of spheroids                  ** c
!c **           2 - surface area  mixture of spheroids           ** c
!c **   keyEL=  0 - do not calculate angular characreristics     ** c
!c **           1 - calculate F11                                ** c
!c **           2 - F11,F12                                      ** c
!c **           3 - F11,F12,F22                                  ** c
!c **           4 - F11,F12,F22,F33                              ** c
!c **           5 - F11,F12,F22,F33,F34                          ** c
!c **           6 - F11,F12,F22,F33,F34,F44                      ** c
!c **   key_org=0 - read original kernels simultaneously make    ** c
!c **               angle interpolation and calculate opt.char.  ** c
!c **           1 -  -"-, save new kernels in                    ** c
!c **                /NEWORG/ directory, STOP                    ** c
!c **        =0 -  create fixed kernels (for fixed axis          ** c
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **         1 - save fixed kernels with original kernel format ** c
!c **             in order to be used as input kernels;          ** c
!c **             'Rke...fix' kernels have to be renamed and moved* c
!c **             into directory 'dir_name'(see 'matrix_fixget.f')* c
!c **             The files can not be used if key=2.            ** c
!c **                                                            ** c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.dat'                                  ** c
!c **           1 - 'grid1.dat.new'                              ** c
!c **   WL   - wavelength                                        ** c
!c **   RN   - real part of the refractive index                 ** c
!c **   RK   - imaginary part of the refractive index            ** c
!c **   KN   - number of grid radii                              ** c
!c **   grid(KN) - grid radii                                    ** c
!c **   SD(KN)   - size distribution for grid radii              ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_F - .fix kernel directory name                  ** c
!c **   distname_N - new original kernel directory name          ** c
!c **   KR  - number of axis ratios                              ** c
!c **   R(KR)  - grid axis ratios                                ** c
!c **   RD(KR) - axis ratio distribution for grid axis ratios    ** c
!c **   KM   - number of scattering angles                       ** c
!c **   ANGLE(KM) - scattering angles                            ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   ext     - extinction                                     ** c
!c **   albedo  - absorption                                     ** c
!c **   f... - scattering matrix elements                        ** c
!c **************************************************************** c
! use iso_c_binding
         USE mo_par_DLS
         USE mo_DLS, only: keyEL, keySUB, keyLS, key_RD, key_RD1 &
                           &, WL, RN, RK, CM, SM, RMM &
                           &, KN, grid, SD &
                           &, KR, R, RD &
                           &, xext, xabs, xsca, albedo &
                           &, KM, ANGLE &
                           &, f11, f12, f22, f33, f34, f44, XLDR, XBLR

!      use mod_os
         use mo_intrpl_linear
         use phase_func
				 use mo_usea
         implicit none

         integer, intent(inout)    :: NDP

!     keySUB =0 original code, !=1 to be used as Subroutine
!     in inversion code

         real dlnr, xnorm

         !real USEA, US11, US12, US22, US33, US34, US44
         !COMMON/US1/US11(KMpar, KNpar)
         !COMMON/US2/US12(KMpar, KNpar)
         !COMMON/US3/US22(KMpar, KNpar)
         !COMMON/US4/US33(KMpar, KNpar)
         !COMMON/US5/US34(KMpar, KNpar)
         !COMMON/US6/US44(KMpar, KNpar)
         !COMMON/US0/USEA(2, KNpar)

         real, dimension(2) :: X, Y
         real PI
         real coeff, g
         integer i, j, j1, j2
         real(4) :: ang_f33, ang_f44

         ang_f33 = 20.0
         ang_f44 = 20.0

         PI = ACOS(-1.)
!c
!c ** GET MATRICES
!c

         CALL USMATRIX(dlnr, NDP)

! coeff=1e+3 because of at the time of kernel calculation cross section (um**2) was multiplied by coeffitient 1e-3 (???)
!
         coeff = 1e+3                  ! for AERONET dV/dlnR
         SD(1:KN) = SD(1:KN)*coeff*dlnr

         xext = DOT_PRODUCT(USEA(1, 1:KN), SD(1:KN))
         xabs = DOT_PRODUCT(USEA(2, 1:KN), SD(1:KN))

         xsca = xext - xabs
         albedo = xsca/xext

         if (keyEL .gt. 0) then
            do j = 1, KM
               f11(j) = DOT_PRODUCT(US11(j, 1:KN), SD(1:KN))
            end do ! j
            if (keySUB .eq. 0) f11(1:KM) = f11(1:KM)/xsca
         end if

         if (keyEL .gt. 1) then
            do j = 1, KM
               f12(j) = DOT_PRODUCT(US12(j, 1:KN), SD(1:KN))
            end do ! j
            if (keySUB .eq. 0) f12(1:KM) = -f12(1:KM)/xsca/f11(1:KM)
         end if

         if (keyEL .gt. 2) then
            do j = 1, KM
               f22(j) = DOT_PRODUCT(US22(j, 1:KN), SD(1:KN))
            end do ! j
            if (keySUB .eq. 0) f22(1:KM) = f22(1:KM)/xsca/f11(1:KM)
         end if
         if (keyEL .gt. 3) then
            do j = 1, KM
               f33(j) = DOT_PRODUCT(US33(j, 1:KN), SD(1:KN))
            end do ! j
            if (keySUB .eq. 0) f33(1:KM) = f33(1:KM)/xsca/f11(1:KM)
         end if
         if (keyEL .gt. 4) then
            do j = 1, KM
               f34(j) = DOT_PRODUCT(US34(j, 1:KN), SD(1:KN))
            end do ! j
            if (keySUB .eq. 0) f34(1:KM) = f34(1:KM)/xsca/f11(1:KM)
         end if

         if (keyEL .gt. 5) then
            do j = 1, KM
               f44(j) = DOT_PRODUCT(US44(j, 1:KN), SD(1:KN))
            end do ! j
            if (keySUB .eq. 0) f44(1:KM) = f44(1:KM)/xsca/f11(1:KM)
         end if

         SD(1:KN) = SD(1:KN)/coeff/dlnr
!c
!c ** SMOOTH f33 and f44
!c
         j1 = 0
         j2 = 0
         if (keyEL .gt. 3) then
            j = 1
            do while (ANGLE(j) .lt. (ang_f33) .or. j .gt. KM)
               j1 = j
               j = j + 1
            end do
            j = 1
            do while (ANGLE(j) .lt. (ang_f33 + 10.0) .or. j .gt. KM)
               j2 = j
               j = j + 1
            end do

            if (j2 .gt. j1) then
               X(1) = ANGLE(j1 - 1)
               X(2) = ANGLE(j2 + 1)
               Y(1) = f33(j1 - 1)
               Y(2) = f33(j2 + 1)
               Y(1) = LOG(Y(1))
               Y(2) = LOG(Y(2))
               do j = j1, j2
                  f33(j) = EXP(LINEAR(X, Y, 2, ANGLE(j)))
               end do ! j
            end if ! j2&j1
         end if ! keyEL>3
!c
         j1 = 0
         j2 = 0
         if (keyEL .gt. 5) then
            j = 1
            do while (ANGLE(j) .lt. (ang_f44) .or. j .gt. KM)
               j1 = j
               j = j + 1
            end do
            j = 1
            do while (ANGLE(j) .lt. (ang_f44 + 10.0) .or. j .gt. KM)
               j2 = j
               j = j + 1
            end do

            if (j2 .gt. j1) then
               X(1) = ANGLE(j1 - 1)
               X(2) = ANGLE(j2 + 1)
               Y(1) = f44(j1 - 1)
               Y(2) = f44(j2 + 1)
               Y(1) = LOG(Y(1))
               Y(2) = LOG(Y(2))
               do j = j1, j2
                  f44(j) = EXP(LINEAR(X, Y, 2, ANGLE(j)))
               end do ! j
            end if ! j2&j1
         end if ! keyEL>5

!c ** Check f11 norma
!c
         IF (keyEL .gt. 0) THEN
            if ((KM .eq. 180 .or. KM .eq. 181) .and. &
                (ANGLE(1) .eq. (0.) .and. ANGLE(KM) .eq. (180.))) then
               if (keySUB .eq. 0) then
                  call SINT(ANGLE, f11, KM, xnorm)
                  call ASYPAR(ANGLE, f11, KM, g)
               else ! keySUB=1
                  call SINT(ANGLE, f11/xsca, KM, xnorm)
                  call ASYPAR(ANGLE, f11/xsca, KM, g)
               end if
            end if ! KM
         END IF ! keyEL>0

!c ** WRITE APPROXIMATED OPTICAL CHARACTERISTICS
!c
         IF (keySUB .eq. 0) THEN

!c
            IF (keyEL .gt. 0) THEN
               if ((KM .eq. 180 .or. KM .eq. 181) .and. &
                   (ANGLE(1) .eq. (0.) .and. ANGLE(KM) .eq. (180.))) then
                  XBLR = 4.*PI/(albedo*f11(KM))
                  if (keyEL .gt. 2) then
                     XLDR = (1.-f22(KM))/(1.+f22(KM))*100.
                  end if
               else if (KM .gt. 0 .and. ANGLE(KM) .eq. 180.) then
                  XBLR = 4.*PI/(albedo*f11(KM))
                  if (keyEL .gt. 2) then
                     XLDR = (1.-f22(KM))/(1.+f22(KM))*100.
                  end if
               end if
            END IF ! keyEL>0

         END IF ! key_SUB=0
13       format('check: f11 norma=', f7.3, '   one=', f7.3)
10       FORMAT('wavelength=', e11.4, '  n=', e11.4, '  k=', e11.4)
11       FORMAT('ext=', e13.5, '  abs=', e13.5, &
                '  sca=', e13.5, '  albedo=', e13.5)
12       FORMAT(F7.2, 6E14.5)
14       FORMAT(A8, 8A13, 2A8)
15       FORMAT(F8.3, 8E13.5, 2F8.2)
         RETURN
      END SUBROUTINE OPTCHAR

!c ***************************************************************** c
