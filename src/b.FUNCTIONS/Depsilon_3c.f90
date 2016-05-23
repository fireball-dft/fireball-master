! copyright info:
!
!                             @Copyright 2012
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman
!
!
! Depsilon_3c.f90
! Program Description
! ===========================================================================
!       This subroutine sets up deps/dratm, deps/dr1, and deps/dr2 in the
! three-center molecular system defined by sighat=(r2-r1)/|r2-r1| and
! piprimehat = sighat-cross-dhat where dvec = ratm -0.5*(r1+r2).
! The eps matrix is obtained by calling epsiln(dhat,sighat).
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
!> @author Barry J. Haycock
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Depsilon_3c (r1, r2, r21, distance12, ratm, rnabc,   &
     &                          eps3, dera3, der13)

        implicit none
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent(in) :: distance12 ! distance between atom 1 and atom 2

        real, intent(in) :: eps3 (3, 3)! 3c epsilon matrix z = r2 - r1, x = dhat
        real, intent(in) :: r1 (3)     ! position of atom 1
        real, intent(in) :: r2 (3)     ! position of atom 2
        real, intent(in) :: r21 (3)    ! vector between atom 1 and atom 2
        real, intent(in) :: ratm (3)   ! position of the potential
        real, intent(in) :: rnabc (3)  ! vector from bond charge to potential

! Output
        real, intent(out) :: dera3 (3, 3, 3) ! deps/dratm in the 3-center system
        real, intent(out) :: der13 (3, 3, 3) ! deps/dr1 in the 3-center system

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer ix

        real crossmag
        real r1dotr2
        real r1dotratm
        real r2dotratm
        real r2mag2
        real r1mag2
        real ratmmag2
        real sum

        real crossa (3)

! Procedure
! ===========================================================================
! Initialize the necessary arrays.
        dera3 = 0.0d0
        der13 = 0.0d0

! Calculate (r2 - r1) - cross - d
        crossa(1) = r21(2)*rnabc(3) - r21(3)*rnabc(2)
        crossa(2) = r21(3)*rnabc(1) - r21(1)*rnabc(3)
        crossa(3) = r21(1)*rnabc(2) - r21(2)*rnabc(1)
        crossmag = sqrt(crossa(1)*crossa(1) + crossa(2)*crossa(2)   &
     &                  + crossa(3)*crossa(3))

! Crossmag cannot be zero for what follows. However, this should rarely happen
! because if crossmag = 0 we probably have an "on-top" instead of a true
! three-center system.
        if (abs(crossmag) .lt. 1.0d-3) then
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*) ' *********** WARNING in deps3cent ************ '
          write (11,*) ' Vectors sighat and dhat dangerously colinear '
          write (11,*) ' sigma - cross - r21 = ', crossmag
          write (11,*) ' setting all 3-center deps to zero '
          close (11)
         return
        end if

        r2mag2 = r2(1)**2 + r2(2)**2 + r2(3)**2
        r1mag2 = r1(1)**2 + r1(2)**2 + r1(3)**2
        ratmmag2 = ratm(1)**2 + ratm(2)**2 + ratm(3)**2
        r1dotr2 = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
        r2dotratm = r2(1)*ratm(1) + r2(2)*ratm(2) + r2(3)*ratm(3)
        r1dotratm = r1(1)*ratm(1) + r1(2)*ratm(2) + r1(3)*ratm(3)

! Now calculate deps/dratm
        do ix = 1, 3
         do imu = 1, 3
          dera3(ix,imu,3) = 0.0d0
          sum = xlevi(ix,imu,1)*(r2(1) - r1(1))  &
     &        + xlevi(ix,imu,2)*(r2(2) - r1(2))  &
     &        + xlevi(ix,imu,3)*(r2(3) - r1(3))
          dera3(ix,imu,2) =  &
     &     (1.0d0/crossmag)*(sum - (eps3(imu,2)/crossmag)  &
     &     *((r1mag2 + r2mag2 - 2.0d0*r1dotr2)*ratm(ix)  &
     &       + (r1dotr2 - r1mag2 - r2dotratm + r1dotratm)*r2(ix)  &
     &       + (r1dotr2 - r2mag2 + r2dotratm - r1dotratm)*r1(ix)))
         end do
        end do

        dera3(1,1,1) = eps3(3,3)*dera3(1,2,2) - eps3(2,3)*dera3(1,3,2)
        dera3(1,2,1) = eps3(1,3)*dera3(1,3,2) - eps3(3,3)*dera3(1,1,2)
        dera3(1,3,1) = eps3(2,3)*dera3(1,1,2) - eps3(1,3)*dera3(1,2,2)

        dera3(2,1,1) = eps3(3,3)*dera3(2,2,2) - eps3(2,3)*dera3(2,3,2)
        dera3(2,2,1) = eps3(1,3)*dera3(2,3,2) - eps3(3,3)*dera3(2,1,2)
        dera3(2,3,1) = eps3(2,3)*dera3(2,1,2) - eps3(1,3)*dera3(2,2,2)

        dera3(3,1,1) = eps3(3,3)*dera3(3,2,2) - eps3(2,3)*dera3(3,3,2)
        dera3(3,2,1) = eps3(1,3)*dera3(3,3,2) - eps3(3,3)*dera3(3,1,2)
        dera3(3,3,1) = eps3(2,3)*dera3(3,1,2) - eps3(1,3)*dera3(3,2,2)

! Now calculate deps/dr1
        do ix = 1, 3
         do imu = 1, 3
          der13(ix,imu,3) = (1.0d0/distance12)                               &
     &                      *(eps3(imu,3)*eps3(ix,3) - delk(imu,ix))
          sum = xlevi(ix,imu,1)*(ratm(1) - r2(1))                            &
     &          + xlevi(ix,imu,2)*(ratm(2) - r2(2))                          &
     &          + xlevi(ix,imu,3)*(ratm(3) - r2(3))
          der13(ix,imu,2) =                                                  &
     &     (1.0d0/crossmag)*(sum - (eps3(imu,2)/crossmag)                    &
     &     *((r2dotratm - r1dotratm + r1dotr2 - r2mag2)*ratm(ix)             &
     &       + (r2dotratm + r1dotratm - r1dotr2 - ratmmag2)*r2(ix)           &
     &       + (ratmmag2 - 2.0d0*r2dotratm + r2mag2)*r1(ix)))
         end do
        end do

        der13(1,1,1) = eps3(3,3)*der13(1,2,2) - eps3(3,2)*der13(1,2,3)       &
     &                 - eps3(2,3)*der13(1,3,2) + eps3(2,2)*der13(1,3,3)
        der13(1,2,1) = eps3(1,3)*der13(1,3,2) - eps3(1,2)*der13(1,3,3)       &
     &                 - eps3(3,3)*der13(1,1,2) + eps3(3,2)*der13(1,1,3)
        der13(1,3,1) = eps3(2,3)*der13(1,1,2) - eps3(2,2)*der13(1,1,3)       &
     &                 - eps3(1,3)*der13(1,2,2) + eps3(1,2)*der13(1,2,3)

        der13(2,1,1) = eps3(3,3)*der13(2,2,2) - eps3(3,2)*der13(2,2,3)       &
     &                 - eps3(2,3)*der13(2,3,2) + eps3(2,2)*der13(2,3,3)
        der13(2,2,1) = eps3(1,3)*der13(2,3,2) - eps3(1,2)*der13(2,3,3)       &
     &                 - eps3(3,3)*der13(2,1,2) + eps3(3,2)*der13(2,1,3)
        der13(2,3,1) = eps3(2,3)*der13(2,1,2) - eps3(2,2)*der13(2,1,3)       &
     &                 - eps3(1,3)*der13(2,2,2) + eps3(1,2)*der13(2,2,3)

        der13(3,1,1) = eps3(3,3)*der13(3,2,2) - eps3(3,2)*der13(3,2,3)       &
     &                 - eps3(2,3)*der13(3,3,2) + eps3(2,2)*der13(3,3,3)
        der13(3,2,1) = eps3(1,3)*der13(3,3,2) - eps3(1,2)*der13(3,3,3)       &
     &                 - eps3(3,3)*der13(3,1,2) + eps3(3,2)*der13(3,1,3)
        der13(3,3,1) = eps3(2,3)*der13(3,1,2) - eps3(2,2)*der13(3,1,3)       &
     &                 - eps3(1,3)*der13(3,2,2) + eps3(1,2)*der13(3,2,3)

! Format Statements
! ===========================================================================

        return
        end
