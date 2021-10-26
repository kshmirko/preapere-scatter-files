!      SUBROUTINE MATRIX_FIX(key,key_RD,keyEL,keySUB
!     &                     ,key_org,key_fx,key_grid1
!     &                     ,KR,R,RD,KN1,grid1,KM,ANGLE
!     &                     ,WAVEL,KRE,KIM,ARE,AIM,pomin,pomax
!     &                     ,NRATN,RATIO
!     &                     ,distname_O,distname_F,distname_N, NDP)
      SUBROUTINE MATRIX_FIX(KN1, grid1, WAVEL, KRE, KIM, ARE, AIM,&
				& RATIO, NDP)
!c ** matrix_fixget_S.f previous file name 30/06/2011
!c ** ANGLE SPLINE interpolation version
!c ** rootdir is defined in "mo_par_DLS.f90"
!c ** 12/04/03 f22 interpolation is logarithmic
!c ** 05/05/03 this version can be used to retrieve an aspect ratio
!c ** 13/10/03 IF(ANGLE<=40.)ln(f33)-interpol.
!c **          IF(ANGLE<=50.)ln(f44)-interpol.
!c ** 16/07/18 IF(ANGLE<=20.)ln(f33)-interpol.
!c **          IF(ANGLE<=20.)ln(f44)-interpol.

!c **************************************************************** c
!c **   Subroutine gets original and calculates fixed kernel     ** c
!c **   matrices for given                                       ** c
!c **   aspect ratio distribution and scattering angles          ** c
!c **                                                            ** c
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
!c **   keyEL=  1 - calculate F11                                ** c
!c **           2 - F11,F12                                      ** c
!c **           3 - F11,F12,F22                                  ** c
!c **           4 - F11,F12,F22,F33                              ** c
!c **           5 - F11,F12,F22,F33,F34                          ** c
!c **           6 - F11,F12,F22,F33,F34,F44                      ** c
!c **   key_org=0 - read original kernels simultaneously make    ** c
!c **               angle interpolation and calculate opt.char.  ** c
!c **           1 -  -"-, save new kernels in                    ** c
!c **                /NEWORG/ directory, STOP                    ** c
!c **   key_fx works when key=1                                  ** c
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
!c **          =0 - 'grid1.txt'                                  ** c
!c **           1 - 'grid1_fix.txt'                              ** c
!c **   KR  - number of aspect ratios                            ** c
!c **   R(KR)  - grid ratios                                     ** c
!c **   RD(KR) - aspect ratio distribution for grid aspect ratios** c
!c **   KM   - number of scattering angles for fixed kernels     ** c
!c **   ANGLE(KM)- scattering angles for fixed kernels           ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_F - .fix kernel directory name                  ** c
!c **   distname_N - new original kernel directory               ** c
!c **                                      name (key_org=1)      ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   UF...(KMpar,KN1par,KIMpar,KREpar) - kernels for          ** c
!c **           given aspect ratio distribution                  ** c
!c **                                                            ** c
!c **   UFEA(2,KN1par,KIMpar,KREpar) - extinction and absorption ** c
!c **                                                    kernels ** c
!c **              1 - extinction                                ** c
!c **              2 - absorption                                ** c
!c **   KN1 - number of grid radii from original or fixed kernels** c
!c **   grid1(KN1) - grid radii                                  ** c
!c **   WAVEL - wavelength from original or fixed kernels        ** c
!c **   KRE   - number of real parts of refr.ind.                ** c
!c **   KIM   - number of imaginary parts of refr.ind.           ** c
!c **   ARE(KRE) - real parts of refr.ind.                       ** c
!c **   AIM(KIM) - imaginary parts of refr.ind.                  ** c
!c **************************************************************** c
!c **************************************************************** c
!c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.txt'                                  ** c
!c **           1 - 'grid1_fix.txt'                              ** c
!c

      use, intrinsic :: IEEE_ARITHMETIC ! fd
      use alloc
      use alloc1
      USE mo_par_DLS
      use mo_DLS, only: key, key_RD, keyEL, keySUB, key_org, key_fx,&
			& key_grid1, KR, R, RD, KM, ANGLE, pomin, pomax, distname_O,&
			& distname_F, distname_N, NRATN, comm_name
!      use mod_os
      use mo_intrpl_linear
      use mo_intrpl_spline

      integer          ::  NEL(0:6)
      character(len=2) ::  NELC(0:6)

      real rmin1, rmax1, rgrid1, rgrid2
      real RATIO(KR1par), RRATN
      dimension grid1(KN1par), ANGLE1(KM1par), ANGLE2(KM1par),&
			& RAB1(KM1par), RAB2(KM1par, KN1par, KIMpar, KREpar)
      dimension RDc(KR), X(2), Y(2), ARE(KREpar), AIM(KIMpar)
      real WAVEL
!      real LINEAR                !AH original
      integer NDP
!cl     &, LINEAR_LN
      CHARACTER(255) name, dir_name_O, dir_name_F, dir_name_N
      CHARACTER(255) full_name
      CHARACTER(255) comm_name1

!c      save NDP
      INTEGER                      :: key_spln, NNEL, KKEL
      double precision, dimension(KM1par)    :: XXS2, YYS2
      double precision                       :: XARG, YFIT
      double precision, dimension(KM1par + 4)  :: KS1, CS1
!c ----- Timer ------
      !real*4, external :: etime,dtime

      real ::  tarray(2), T_RFM, T_CFM, T_CFM0
      real   ::  UTEST, tiny = 1e-29
      real :: xx1, xx2
      real :: ang_f33, ang_f44

      ang_f33 = 20.0
      ang_f44 = 20.0

!cl    write(*,*)'distname_O=',TRIM(distname_O)
!cl          write(*,*)'distname_F=',TRIM(distname_F)
!cl          write(*,*)'distname_N=',TRIM(distname_N)

      !dir_name_O=TRIM(rootdir)//TRIM(distname_O)//'/'
      !dir_name_F=TRIM(rootdir)//TRIM(distname_F)//'/'
      !dir_name_N=TRIM(rootdir)//TRIM(distname_N)//'/'

      dir_name_O = TRIM(distname_O)//'/'
      dir_name_F = TRIM(distname_F)//'/'
      dir_name_N = TRIM(distname_N)//'/'

      PI2 = 2.*ACOS(-1.)
      NEL = (/0, 11, 12, 22, 33, 34, 44/)
      NELC = (/'00', '11', '12', '22', '33', '34', '44'/)
      T_RFM = 0.
      T_CFM = 0.
      !write (*, *) 'key_grid1=', key_grid1

      if (key_grid1 .eq. 0) then
         full_name = TRIM(dir_name_O)//'grid1.txt'
      else
         full_name = TRIM(dir_name_F)//'grid1_fix.txt'
      end if ! key_grid1
      !write (*, *) 'file=', TRIM(full_name)
      OPEN (10, file=full_name, status='unknown')
      READ (10, *) KN1, XWL

      if (KN1 .gt. KN1par) then
         write (*, *) 'in GET_MATRIX: KN1=', KN1, ' KN1par=', KN1par
         stop ' !!! KN1.ne.KN1par'
      end if
      DO I = 1, KN1
         READ (10, *) grid1(I)
      END DO ! I
      rgrid1 = grid1(1)
      rgrid2 = grid1(KN1)
      READ (10, *) KM1
      if (KM1 .gt. KM1par) then
         write (*, *) 'in GET_MATRIX: KM1=', KM1, ' KM1par=', KM1par
         stop ' !!! KM1.ne.KM1par'
      end if
      DO J = 1, KM1
         READ (10, *) ANGLE1(J)
      END DO ! J
      CLOSE (10)
      if ((key .eq. 4) .and. (KM .ne. KM1)) then
         write (*, *) 'if key=4, in input.dat KM=', KM,&
				 &' should be equal to KM1=', KM1, ' in grid1.txt'
         stop 'STOP in MATRIX_FIX (matrix_fixget.f)'
      end if ! key
!c
!c ** Redefine grids for fixed or NEWORG kernels
!c
      NDPP = 0
      if (key .eq. 1 .or. key_org .eq. 1) then
         rgmin = pomin*XWL/pi2
         rgmax = pomax*XWL/pi2
!tl        write(*,*) 'pomin=',pomin,' XWL=',XWL,' pi2=',pi2,rgmin
!tl        stop
!c
         xx1 = rgmin
         xx2 = grid1(1)
         if (xx1 .lt. xx2) then
            write (*, *) 'XWL=', XWL,&
             &' IS wave length in kernels equal to XWL?'
            write (*, *) 'check input.dat: rgmin=', rgmin, &
						&' < ','grid1(1)=', grid1(1)
            STOP 'STOP in MARTRIX_FIX (matrix_fixget.f)'
         end if
         xx1 = rgmax
         xx2 = grid1(KN1)
         if (xx1 .gt. xx2) then
            write (*, *) 'XWL=', XWL,&
             &' IS wave length in kernels equal to XWL?'
            write (*, *) 'check input.dat: rgmax=', rgmax, ' > ',&
						 &'grid1(KN1)=', grid1(KN1)
            STOP 'STOP in MARTRIX_FIX (matrix_fixget.f)'
         end if

         ind = 0
         do i = 1, KN1 - 1
!c        write(*,*) 'i=',i,' rgmin=',rgmin,' grid1(i+1)=',grid1(i+1)
            if (rgmin .le. grid1(i + 1)) then
!c        xr1=grid1(i)
               nn1 = i
               ind = ind + 1
            end if
            if (ind .eq. 1) EXIT
         end do ! i
!c
         ind = 0
         do i = KN1, 2, -1
!c        write(*,*) 'i=',i,' rgmax=',rgmax,' grid1(i-1)=',grid1(i-1)
            if (rgmax .ge. grid1(i - 1)) then
!c        xr2=grid1(i)
               nn2 = i
               ind = ind + 1
            end if
            if (ind .eq. 1) EXIT
         end do ! i
!c
         if (rgmin .eq. grid1(1)) then
!c        xr1=grid1(1)
            nn1 = 1
         end if ! rgmin
         if (rgmax .eq. grid1(KN1)) then
!c        xr2=grid1(KN1)
            nn2 = KN1
         end if ! rgmax
         nn3 = nn2 - nn1 + 1
!cl        write(*,*) 'nn1=',nn1,' nn2=',nn2,' nn3=',nn3
!cl        write(*,*) 'xr1=',xr1,' xr2=',xr2,' KN1=',KN1

!cl        STOP 'TEST STOP'
!c
!      write(*,*) 'before grid1_fix.txt'
         if (key_org .eq. 0) then
            full_name = TRIM(dir_name_F)//'grid1_fix.txt'
         else
            full_name = TRIM(dir_name_N)//'grid1.txt'
         end if
         write (*, *) 'file=', TRIM(full_name)
         open (77, file=full_name, status='replace')
         WRITE (77, '(i4,f8.3)') nn3, XWL
         DO I = nn1, nn2
            WRITE (77, '(e15.7)') grid1(I)
         END DO ! I
         WRITE (77, *) KM
         DO J = 1, KM
            WRITE (77, '(f6.2)') ANGLE(J)
         END DO ! J
         close (77)
         WRITE (*, *)
         WRITE (*, *) '  ATTENTION:'
         WRITE (*, *) '    New input file grid1.txt (key=1,key_org=1) or'
         WRITE (*, *) '    grid1_fix.txt (key=1,key_org=0) has been created'
         WRITE (*, *) '    for your further calculations !!!'
         WRITE (*, *)
      end if ! key&key_org
!			этот кусок оставляем
			
!c
!c ** Define directories
!c
!c ** for AERONET
!cl      NNEL=2
!cl      if(key_f11.eq.0) NNEL=7
      NNEL = keyEL + 1
      KKEL = 6 - keyEL

      IF (key .EQ. 2) THEN

         open (10, file=TRIM(TRIM(dir_name_F)//"Rkernel1_fix_00.txt"),&
         &status = "old")
         IF (keyEL .gt. 0) THEN
            open(11,file=TRIM(TRIM(dir_name_F)//"Rkernel1_fix_11.txt"),&
						&status = "old")
            if (keyEL .gt. 1) open (12, file=TRIM(dir_name_F)//"Rkernel1_fix_12.txt",&
            &status = "old")
            if (keyEL .gt. 2) open (13, file=TRIM(dir_name_F)//"Rkernel1_fix_22.txt",&
            &status = "old")
            if (keyEL .gt. 3) open (14, file=TRIM(dir_name_F)//"Rkernel1_fix_33.txt",&
            &status = "old")
            if (keyEL .gt. 4)open (15, file=TRIM(dir_name_F)//"Rkernel1_fix_34.txt",&
            &status = "old")
            if (keyEL .gt. 5)open (16, file=TRIM(dir_name_F)//"Rkernel1_fix_44.txt",&
            &status = "old")
         END IF ! keyEL
         open (20, file='CHECK.dat', status='unknown')

         DO II = 1, NNEL
            WRITE (20, *) 'INPUT CHECK: Nelement=', NEL(II - 1)
            READ (10 + II - 1, *) key_RD1
            IF (key_RD .ne. key_RD1) then
               WRITE (20, 14) II, key_RD, key_RD1
               WRITE (20, *) 'STOP: key_RD.ne.keyRD1 in Rke...fix'
               WRITE (*, 14) II, key_RD, key_RD1
               WRITE (*, *) 'STOP: key_RD.ne.keyRD1 in Rke...fix'
               STOP
            END IF ! key_RD

            READ (10 + II - 1, *) rmin, rmax
            rmin1 = rmin
            rmax1 = rmax
            IF (rmin1 .ne. rgrid1 .or. rmax1 .ne. rgrid2) THEN
               WRITE (20, *) 'grid1(1)=', rgrid1, ' rmin=', rmin1,&
               &' rmax=', rmax1, ' grid1(KN1)=', rgrid2
               WRITE (*, *) 'grid1(1)=', rgrid1, ' rmin=', rmin1,&
               &' rmax=', rmax1, ' grid1(KN1)=', rgrid2
               WRITE (*, *) 'Compare1: grid1(1) in grid1.dat',&
               &' with RMIN in Rke...fix'
               WRITE (*, *) 'Compare2: grid1(KN1) in grid1.dat',&
               &' with RMAX in Rke...fix'
               STOP 'STOP: grid1(1)/grid(KN1).ne.rmin/rmax'
            END IF
            READ (10 + II - 1, *) kk
            WRITE (20, *) 'KR=', KR, ' KR=', kk
            READ (10 + II - 1, *)
            DO I = 1, KR
               READ (10 + II - 1, *) xx, yy
               WRITE (20, *) R(I), RD(I), ' R,RD', xx, yy, ' R,RD'
            END DO ! I
            READ (10 + II - 1, *) kk
            WRITE (20, *) 'KN1=', KN1, ' KN1=', -kk
            IF (KN1 .ne. -kk) THEN
               WRITE (*, *) 'in grid1.dat KN1=', KN1,&
               &' .ne. KN1=', -kk, ' in Rke...fix'
               STOP 'STOP in matrix_fixget'
            END IF
         END DO ! II

         DO II = 2, NNEL
            WRITE (20, *) 'INPUT CHECK: Nelement=', NEL(II - 1)
            READ (10 + II - 1, *) kk
            WRITE (20, *) 'KM=', KM, ' KM=', kk
            READ (10 + II - 1, *) ANGLE2(:KM)
            WRITE (20, *) 'SCATTERING ANGLES2:'
            WRITE (20, 15) ANGLE2(:KM)
            WRITE (20, *) 'SCATTERING ANGLES:'
            WRITE (20, 15) ANGLE(:KM)
         END DO ! II

         DO II = 1, NNEL
            WRITE (20, *) 'INPUT CHECK: Nelement=', NEL(II - 1)
            READ (10 + II - 1, *) xx, yy
            WRITE (20, *) 'ARE(1),ARE(KRE)', xx, yy
            READ (10 + II - 1, *) xx, yy
            WRITE (20, *) 'AIM(1),AIM(KIM)', xx, yy
            READ (10 + II - 1, *) KRE, KIM
            WRITE (20, *) 'KRE,KIM', KRE, KIM
            IF (KIM .lt. 0) KIM = -KIM
         END DO ! II
         !WRITE (*, *) 'Fixed kernel matrices have been read :'
         DO IRE = 1, KRE
            DO IIM = 1, KIM
               DO II = 1, NNEL
                  if (key_fx .eq. 1) then
                     READ (10 + II - 1, *) kk, aa
                  else ! key_fx=0
                     READ (10 + II - 1, *) kk
                  end if ! key_fx
                  READ (10 + II - 1, *) WAVEL, ARE(IRE), AIM(IIM)
                  IF (AIM(IIM) .LT. 0) AIM(IIM) = -AIM(IIM)
                  !cl        WRITE(20,*) 'WAVEL,ARE,AIM,NEL',
                  !cl     &                WAVEL,ARE(IRE),-AIM(IIM),NEL(II-1)
               END DO ! II

               READ (10, *)
               READ (10, 11) UFEA(1, :KN1, IIM, IRE)
               READ (10, *)
               READ (10, 11) UFEA(2, :KN1, IIM, IRE)

               if (keyEL .gt. 0) then
                  DO I = 1, KN1
                     READ (11, 11) UF11(1:KM, I, IIM, IRE)
                  END DO ! I
               end if
               if (keyEL .gt. 1) then
                  DO I = 1, KN1
                     READ (12, 11) UF12(1:KM, I, IIM, IRE)
                  END DO ! I
               end if
               if (keyEL .gt. 2) then
                  DO I = 1, KN1
                     READ (13, 11) UF22(1:KM, I, IIM, IRE)
                  END DO ! I
               end if

               if (keyEL .gt. 3) then
                  DO I = 1, KN1
                     READ (14, 11) UF33(1:KM, I, IIM, IRE)
                  END DO ! I
               end if
               if (keyEL .gt. 4) then
                  DO I = 1, KN1
                     READ (15, 11) UF34(1:KM, I, IIM, IRE)
                     do j = 1, KM
                        UTEST = UF34(j, I, IIM, IRE)
                        if ((abs(UTEST) .lt. tiny) .and. (abs(UTEST) .gt. 0.0)) UF34(j, I, IIM, IRE) = 0.0
                     end do ! j

                  END DO ! I
               end if
               if (keyEL .gt. 5) then
                  DO I = 1, KN1
                     READ (16, 11) UF44(1:KM, I, IIM, IRE)

                  END DO ! I
               end if

            END DO ! IIM
            !WRITE (*, 16) WAVEL, ARE(IRE)
         END DO ! IRE
         close (10)
         if (keyEL .gt. 0) close (11)
         if (keyEL .gt. 1) close (12)
         if (keyEL .gt. 2) close (13)
         if (keyEL .gt. 3) close (14)
         if (keyEL .gt. 4) close (15)
         if (keyEL .gt. 5) close (16)
         close (20)

         IF (keyEL .gt. 0) THEN
            UF11(1:KM, 1:KN1, 1:KIM, 1:KRE) =LOG(UF11(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
            if (keyEL .gt. 1) UF12(1:KM, 1:KN1, 1:KIM, 1:KRE) = UF12(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
            if (keyEL .gt. 2) UF22(1:KM, 1:KN1, 1:KIM, 1:KRE) = LOG(UF22(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
            if (keyEL .gt. 3) then
               DO J = 1, KM
                  UF33(J, 1:KN1, 1:KIM, 1:KRE) = UF33(J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
                  IF (ANGLE(J) .LE. ang_f33)UF33(J, 1:KN1, 1:KIM, 1:KRE) = LOG(UF33(J, 1:KN1, 1:KIM, 1:KRE))
               END DO ! J
            end if
            if (keyEL .gt. 4) UF34(1:KM, 1:KN1, 1:KIM, 1:KRE) = UF34(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
            if (keyEL .gt. 5) then
               DO J = 1, KM
                  UF44(J, 1:KN1, 1:KIM, 1:KRE) = UF44(J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
                  IF (ANGLE(J) .LE. ang_f44) UF44(J, 1:KN1, 1:KIM, 1:KRE) = LOG(UF44(J, 1:KN1, 1:KIM, 1:KRE))
               END DO ! J
            end if
         END IF ! keyEL

13       FORMAT(3E12.4, I4)
14       FORMAT('II=', i2, ' key_RD=', i2, ' key_RD1=', i2)
16       FORMAT(12x, 'wl=', f5.2, 2x, 'n=', f8.5)
!cl      WRITE(*,*) 'Fixed kernel matrices have been read'
         IF (key_RD .EQ. 1) WRITE (*, *) 'Volume mixture of spheroids'
         IF (key_RD .EQ. 2) WRITE (*, *) 'Surface area mixture of spheroids'

61       format('  Read fixed matr. ........ ', f8.3, ' min.')
      END IF ! key=2




10    FORMAT(3E15.7, ' rmin, rmax, RATIO')
11    FORMAT(7E15.7)
12    FORMAT(3E15.7, I4, '    wavelength, n, k, NEL')
15    FORMAT(7F12.2)

      RETURN
   END SUBROUTINE MATRIX_FIX

!c ****************************************************************
