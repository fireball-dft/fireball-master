! copyright info:
!
!                             @Copyright 2016
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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! Dassemble_3c.f90
! Program Description
! ===========================================================================
!!      This routine assembles all of the true three-center interactions.
!! These are true three-center in that iatom .ne. jatom .ne. katom.
!! The matrix elements look like <psi(1)|V(3)|psi(2)>.

!! It contains the following subroutines within the module:
!!
!!       Dassemble_vna_3c.f90 - assemble neutral atom potential matrix derivatives
!!
!! For a complete list of the interactions see the files 3c.Z1.Z2.Z3.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================

! Module Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
!! A third party term is when we consider the NA (or etc.) to be at the
!! origin and take the matrix element between a pair of neighbors, neither
!! of which is the NA (or etc.), and the pair is truly a pair, and not an
!! atom.
!
!! The generic form of the matrix element is
!! v = <mu(r-(r1+ratm)) | v(r-ratm) | nu(r-(r2+ratm))>.
!! We want the derivative with respect to ratm keeping the position of
!! orbitals mu, nu fixed (constant r1 + ratm).
!
!! The derivative f3 looks like for the neutral atom case:
!! f3na = - sum (all neighbors of atom alpha at (li,bi) but bi.ne.balph)
!!    * sum (all neighbors m of (li,bi), and not having b value balph)
!!    * rho(mu,nu,i,m)* deriv wrt balpha <i! v(balph) !j>.
!! Note the minus sign to make it "force-like".
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms. Breifly fa=-d/dratm, fb=-d/d(bi), fc=-d/d(bj)
! where matrix elements are < psi(r-bi) ! v(r-ratm) ! psi(r-bj)>
!
! ===========================================================================
        module M_Dassemble_3c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_3c
        use M_neighbors
        use M_rotations
        use M_Drotations
        use M_assemble_3c
        use M_density_matrix

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! Dassemble_vna_3c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates derivative of the neutral atom potential
! matrix interactions.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
!> @author Barry Haycock
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_vna_3c (s)

        implicit none
        include '../include/constants.h'
        include '../include/interactions_3c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom     !< the three parties involved
        integer ibeta, jbeta             !< cells for three atoms
        integer ineigh, mneigh           !< counter over neighbors
        integer in1, in2, in3            !< species numbers
        integer interaction, isorp       !< which interaction and subtype
        integer issh
        integer norb_mu, norb_nu         !< size of the block for the pair

        real distance_13, distance_23  !< distance from 3rd atom
        real dQ                          !< net charge on atom
        real rcutoff1_min, rcutoff2_min, rcutoff3_min, rend1, rend2  !< for smoothing
        real stinky13                    !< smoothing value
        real stinky23                    !< smoothing value
        real dstinky13                   !< smoothing value
        real dstinky23                   !< smoothing value
        real stinky                      !< smoothing value
        real, dimension (3) :: dstnA     !< smoothing value
        real, dimension (3) :: dstnB     !< smoothing value
        real, dimension (3) :: dstnC     !< smoothing value
        real xsmooth                     !< for smoothing function

        real z                           !< distance between r1 and r2
        real x, cost                     !< dnabc and angle

        integer ix                       !< dimension - x y or z

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2, r3, r12 !< positions
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - r3
        real, dimension (3) :: rhatA1    !< unit vector along r3 - r1
        real, dimension (3) :: rhatA2    !< unit vector along r3 - r2

        real, dimension (3, 3, 3) :: depsA  !< the Depsilon matrix for the bc
        real, dimension (3, 3, 3) :: depsB  !< the Depsilon matrix for the na

        real, dimension (3) :: amt, bmt
        integer imu, inu, mvalue, iindex

! bcnam = Hartree matrix in molecular coordinates
! bcnax = Hartree matrix in crystal coordinates
! dbcnam = derivative of Hartree matrix in molecular coordinates
! vdbcnam = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcnax = vectorized derivative of Hartree matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcnam
        real, dimension (:, :), allocatable :: bcnax
        real, dimension (:, :), allocatable :: dterm
        real, dimension (:, :), allocatable :: sterm
        real, dimension (:, :), allocatable :: emnpl

        real, dimension (:, :), allocatable :: dpbcnam
        real, dimension (:, :), allocatable :: dxbcnam
        real, dimension (:, :), allocatable :: dybcnam

        real, dimension (:, :, :), allocatable :: f3naMa
        real, dimension (:, :, :), allocatable :: f3naMb

        real, dimension (:, :, :), allocatable :: f3naXa
        real, dimension (:, :, :), allocatable :: f3naXb
        real, dimension (:, :, :), allocatable :: f3naXc
        real, dimension (:, :, :), allocatable :: f3naXat !temporary storage
        real, dimension (:, :, :), allocatable :: f3naXbt !temporary storage
        real, dimension (:, :, :), allocatable :: f3naXct !temporary storage

        real, dimension (:, :, :), allocatable :: DemnplA
        real, dimension (:, :, :), allocatable :: DemnplB
        real, dimension (:, :, :), allocatable :: DemnplC

        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z
        type(T_Fdata_cell_3c), pointer :: pFdata_cell
        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        interface
          function smoother (z, rend, x)
            real smoother
            real, intent(in) :: z, rend, x
          end function smoother
        end interface

        interface
           real function Dsmoother (r, rend, x)
             real, intent (in) :: r
             real, intent (in) :: rend
             real, intent (in) :: x
           end function Dsmoother
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          in3 = s%atom(ialpha)%imass
          r3 = s%atom(ialpha)%ratom

          ! loop over the common neigbor pairs of ialp
          do ineigh = 1, s%neighbors(ialpha)%ncommon
            mneigh = s%neighbors(ialpha)%neigh_common(ineigh)
            if (mneigh .ne. 0) then
              iatom = s%neighbors(ialpha)%iatom_common_j(ineigh)
              ibeta = s%neighbors(ialpha)%iatom_common_b(ineigh)
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max

              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max

              poverlap=>s%overlap(iatom)
              pdipole_z=>s%dipole_z(iatom)
              pvna=>s%vna(iatom)

              pS_neighbors=>poverlap%neighbors(mneigh)
              pdip_neighbors=>pdipole_z%neighbors(mneigh)
              pvna_neighbors=>pvna%neighbors(mneigh)

! Allocade DFdata blocks
              allocate (pvna_neighbors%aDblock(3, norb_mu, norb_nu))
              pvna_neighbors%aDblock = 0.00d0
              allocate (pvna_neighbors%bDblock(3, norb_mu, norb_nu))
              pvna_neighbors%bDblock = 0.00d0
              allocate (pvna_neighbors%cDblock(3, norb_mu, norb_nu))
              pvna_neighbors%cDblock = 0.00d0

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
              z = distance (r1, r2)

              ! unit vector in sigma direction.
              if (z .lt. 1.0d-05) then
                sighat(1) = 0.0d0
                sighat(2) = 0.0d0
                sighat(3) = 1.0d0
              else
                sighat = (r2 - r1)/z
              end if

! ****************************************************************************
! Find rnabc = vector pointing from center of bondcharge to r3
! This gives us the distance dnabc (or x value in the 2D grid).
              r12 = 0.5d0*(r1 + r2)
              x = distance (r12, r3)

! Find other distances -
              distance_13 = distance (r3, r1)
              distance_23 = distance (r3, r2)

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = (r3 - 0.5d0*(r1 + r2))/x
              end if

! Find the unit vector in r3-r1 direction.
              if (distance_13 .gt. 1.0d-05) then
                rhatA1(:) = (r3(:) - r1(:))/distance_13
              end if

! Find the unit vector in r3-r2 direction.
              if (distance_23 .gt. 1.0d-05) then
                rhatA2(:) = (r3(:) - r2(:))/distance_23
              end if

              cost = dot_product(sighat, rhat)
              if (abs(cost) - 1.0d0 .lt. 0.010d0) then                
                cycle
              end if              
              call epsilon_function (rhat, sighat, eps)
              allocate (sterm (norb_mu, norb_nu)); sterm =0.0d0
              allocate (dterm (norb_mu, norb_nu)); dterm =0.0d0
              allocate (emnpl (norb_mu, norb_nu)); emnpl =0.0d0
              allocate (DemnplA (3, norb_mu, norb_nu)); DemnplA =0.0d0
              allocate (DemnplB (3, norb_mu, norb_nu)); DemnplB =0.0d0
              allocate (DemnplC (3, norb_mu, norb_nu)); DemnplC =0.0d0

! Find the smoothing quantity - here we calculate the long-range effective
! monopole.  This term is included so that we obtain no discontinuities when
! atoms leave or enter the rcutoff_1 + rcutoff_2 range criteria.
! Therefore, "close" two-center interactions are exact, while more distant
! two-center integrals go to effective monopoles.  The monopoles are effective
! in the sense that the two atoms in the matrix element, each has a different
! charge.  Since they are separated, this gives a monopole contribution at long
! range.

! The smoothing function is found by calling smoother(r,rbegin,rend).
! We define our final matrix element answer as
! smoother(r)*exact_piece + (1 - smoother(r))*longrange.  The distance r is the
! distance of the third center from the "effective" center of the bondcharge.
! The effective center of the bondcharge is (d + rc1 - rc2)/2 from r1 except in
! weird cases (see below). The distance rbegin is the distance at which we
! include only exact answers and do not smooth. The distance rend is the
! distance past which smooth(r) is zero, so that the result is long-range only.
! We skipped self-interaction terms.
              xsmooth = 0.8d0  ! parameter for smoothing

              rcutoff1_min = 99.0d0
              do issh = 1, species(in1)%nssh
                rcutoff1_min = min(rcutoff1_min, species(in1)%shell(issh)%rcutoffA)
              end do

              rcutoff2_min = 99.0d0
              do issh = 1, species(in2)%nssh
                rcutoff2_min = min(rcutoff2_min, species(in2)%shell(issh)%rcutoffA)
              end do

              rcutoff3_min = 99.0d0
              do issh = 1, species(in3)%nssh
                rcutoff3_min = min(rcutoff3_min, species(in3)%shell(issh)%rcutoffA)
              end do

              rend1 = rcutoff1_min + rcutoff3_min
              stinky13 = smoother (distance_13, rend1, xsmooth)
              dstinky13 = Dsmoother (distance_13, rend1, xsmooth)

              rend2 = rcutoff2_min + rcutoff3_min
              stinky23 = smoother (distance_23, rend2, xsmooth)
              dstinky23 = Dsmoother (distance_23, rend2, xsmooth)

              stinky = stinky13*stinky23

              dstnB = dstinky13*stinky23*rhatA1
              dstnC = stinky13*dstinky23*rhatA2
              dstnA = - dstnB - dstnC

! Set value for emnpl
              sterm = pS_neighbors%block/2.0d0
              dterm = pdip_neighbors%block/z
              emnpl = (sterm - dterm)/distance_13 + (sterm + dterm)/distance_23

              do inu = 1, norb_nu
                do imu= 1, norb_mu
                  DemnplB(:,imu,inu) =                                       &
                    pS_neighbors%Dblock(:,imu,inu)/(2.0d0*distance_13)       &
     &              + pS_neighbors%Dblock(:,imu,inu)/(2.0d0*distance_23)     &
     &              + pdip_neighbors%Dblock(:,imu,inu)/(z*distance_13)       &
     &              - pdip_neighbors%Dblock(:,imu,inu)/(z*distance_23)       &
     &              - pdip_neighbors%block(imu,inu)*sighat(:)/(z**2*distance_13)&
     &              + pdip_neighbors%block(imu,inu)*sighat(:)/(z**2*distance_23)&
     &              + rhatA1(:)*(sterm(imu,inu) - dterm(imu,inu))/distance_13**2

                  DemnplC(:,imu,inu)=                                        &
                    - pS_neighbors%Dblock(:,imu,inu)/(2.0d0*distance_13)     &
     &              - pS_neighbors%Dblock(:,imu,inu)/(2.0d0*distance_23)      &
     &              - pdip_neighbors%Dblock(:,imu,inu)/(z*distance_13)       &
     &              + pdip_neighbors%Dblock(:,imu,inu)/(z*distance_23)       &
     &              + pdip_neighbors%block(imu,inu)*sighat(:)/(z**2*distance_13)&
     &              - pdip_neighbors%block(imu,inu)*sighat(:)/(z**2*distance_23)&
     &              + rhatA2(:)*(sterm(imu,inu) - dterm(imu,inu))/distance_23**2

                end do !imu= 1, norb_mu
              end do !inu = 1, norb_nu
              DemnplA = -DemnplB - DemnplC

              allocate(bcnam(norb_mu, norb_nu)); bcnam = 0.00d0
              allocate(bcnax(norb_mu, norb_nu)); bcnax = 0.00d0
              allocate(dpbcnam(norb_mu, norb_nu)); dpbcnam = 0.00d0
              allocate(dxbcnam(norb_mu, norb_nu)); dxbcnam = 0.00d0
              allocate(dybcnam(norb_mu, norb_nu)); dybcnam = 0.00d0

              allocate (f3naMa(3, norb_mu, norb_nu)); f3naMa= 0.00d0
              allocate (f3naMb(3, norb_mu, norb_nu)); f3naMb= 0.00d0
              allocate (f3naXa(3, norb_mu, norb_nu)); f3naXa= 0.00d0
              allocate (f3naXb(3, norb_mu, norb_nu)); f3naXb= 0.00d0
              allocate (f3naXc(3, norb_mu, norb_nu)); f3naXc= 0.00d0

              allocate (f3naXat(3, norb_mu, norb_nu)); f3naXat= 0.00d0
              allocate (f3naXbt(3, norb_mu, norb_nu)); f3naXbt= 0.00d0
              allocate (f3naXct(3, norb_mu, norb_nu)); f3naXct= 0.00d0

! For now we just do the neutral atom interactions.
! Charged atom interactions are assembled in assemble_ca_3c.f
! So set isorp = 0 within this subroutine.
!
!              interaction    subtypes     index
!
!      bcna         1           0..9(max)   1..10

! ****************************************************************************
!
! Get the D-matrices from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_bcna
! neutral plart
              isorp = 0
              call getDMEs_Fdata_3c (in1, in2, in3, interaction, isorp, x,     &
     &                              z, norb_mu, norb_nu, cost, rhat, sighat,   &
     &                                   bcnam, dpbcnam, dxbcnam, dybcnam)

              ! Rotate into crystal coordinates
              call rotate (in1, in2, eps, norb_mu, norb_nu, bcnam, bcnax)

! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.
              do ix = 1, 3

! The first piece will be the force with respect to atom 3.
                if (x .gt. 1.0d-5) then
                  amt(ix) = (sighat(ix) - cost*rhat(ix))/x
                else
                  amt = 0.0d0
                end if

                pFdata_bundle => Fdata_bundle_3c(in1, in2, in3)
                pFdata_cell =>                                                &
     &          pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(interaction,isorp,1))

                do iindex = 1, pFdata_cell%nME
                  imu = pFdata_cell%mu_3c(iindex)
                  inu = pFdata_cell%nu_3c(iindex)
                  mvalue = pFdata_cell%mvalue_3c(iindex)
! Now recover f3naMa which is a two-dimensional array
                  f3naMa(ix,imu,inu) = rhat(ix)*Dxbcnam(imu,inu)               &
      &                               + amt(ix)*Dpbcnam(imu,inu)

! The second piece will be the force with respect to atom 1.
                  bmt(ix) = (cost*sighat(ix) - rhat(ix))/z

! Now recover f3naMb which is a two-dimensional array
                  f3naMb(ix,imu,inu) = - sighat(ix)*dybcnam(imu, inu)          &
     &              + bmt(ix)*dpbcnam(imu, inu) - f3naMa(ix,imu,inu)/2.0d0
                end do !iindex
              end do !ix


! ***************************************************************************
! Convert to Crystal Coordinates
! ***************************************************************************
! The call to rotated does the rotations to crystal coordinates of these
! force things.
!
! For example:
! Suppose we have f_A(3,mu,nu), which is d/dratm M(mu,nu) where M(mu,nu)
! is in molecular. To transform M(mu,nu) to crystal, we need Udag * M * U.
! Therefore, f_A(3,mu,nu)[CRYSTAL] = (d/dratm Udag) * M * U
!                                   + Udag * M * (d/dratm U)
!                                   + Udag * f_A * U.
!
! So, to use this baby, put in deps3c (deps/dr1, deps/dr2, deps/dratm),
! and f_A and M.
!
! NOTE: rotated works on the assumption that we are adding derivatives,
! NOT forces. So f3naMa,... etc. MUST not yet be forcelike.
! We do the - sign for forces at the end.
! ***************************************************************************
! Force on the neutral atom with respect to atom 3 (f3naMa).
              call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu,          &
     &                      bcnam, f3naMa, f3naXa)

! Force on the neutral atom with respect to atom 1 (f3naMb).
              call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu,          &
     &                      bcnam, f3naMb, f3naXb)

! Make things force-like and determine f3naXc, whcih is found from Newtons Laws:
              f3naXa(:,:,:) = - f3naXa(:,:,:)
              f3naXb(:,:,:) = - f3naXb(:,:,:)
              f3naXc(:,:,:) = - f3naXa(:,:,:) - f3naXb(:,:,:)

!until here calculated in the same way
              pvna_neighbors%aDblock = pvna_neighbors%aDblock + f3naXa*P_eq2
              pvna_neighbors%bDblock = pvna_neighbors%bDblock + f3naXb*P_eq2
              pvna_neighbors%cDblock = pvna_neighbors%cDblock + f3naXc*P_eq2

!*********************************************************************************************************
!*********************************************************************************************************
! Charged atom cases
              do isorp = 1, species(in3)%nssh
                bcnam = 0.00d0; bcnax = 0.00d0; dpbcnam = 0.00d0; dxbcnam = 0.00d0
                dybcnam = 0.00d0; f3naMa= 0.00d0; f3naMb= 0.00d0; f3naXa= 0.00d0
                f3naXb= 0.00d0; f3naXc= 0.00d0; f3naXat= 0.00d0; f3naXbt= 0.00d0
                f3naXct= 0.00d0
                call getDMEs_Fdata_3c (in1, in2, in3, interaction, isorp, x,   &
     &                                 z, norb_mu, norb_nu, cost, rhat, sighat,&
     &                                 bcnam, dpbcnam, dxbcnam, dybcnam)

              ! Rotate into crystal coordinates
                call rotate (in1, in2, eps, norb_mu, norb_nu, bcnam, bcnax)

                do ix = 1, 3

! The first piece will be the force with respect to atom 3.
                  if (x .gt. 1.0d-5) then
                    amt(ix) = (sighat(ix) - cost*rhat(ix))/x
                  else
                    amt = 0.0d0
                  end if

                  pFdata_bundle => Fdata_bundle_3c(in1, in2, in3)
                  pFdata_cell =>                                             &
     &            pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(interaction,isorp,1))

                  do iindex = 1, pFdata_cell%nME
                    imu = pFdata_cell%mu_3c(iindex)
                    inu = pFdata_cell%nu_3c(iindex)
                    mvalue = pFdata_cell%mvalue_3c(iindex)
! Now recover f3naMa which is a two-dimensional array
                    f3naMa(ix,imu,inu) = rhat(ix)*Dxbcnam(imu,inu)            &
      &                                 + amt(ix)*Dpbcnam(imu,inu)

! The second piece will be the force with respect to atom 1.
                    bmt(ix) = (cost*sighat(ix) - rhat(ix))/z

! Now recover f3naMb which is a two-dimensional array
                    f3naMb(ix,imu,inu) = - sighat(ix)*dybcnam(imu, inu)       &
     &                + bmt(ix)*dpbcnam(imu, inu) - f3naMa(ix,imu,inu)/2.0d0

                  end do !iindex
                end do !ix

                call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu,         &
     &                        bcnam, f3naMa, f3naXa)

! Force on the neutral atom with respect to atom 1 (f3naMb).
                call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu,         &
     &                        bcnam, f3naMb, f3naXb)

! Make things force-like and determine f3naXc, whcih is found from Newtons Laws:
                f3naXa(:,:,:) = - f3naXa(:,:,:)
                f3naXb(:,:,:) = - f3naXb(:,:,:)
                f3naXc(:,:,:) = - f3naXa(:,:,:) - f3naXb(:,:,:)

!Ensenble the exact part
                f3naXat = f3naXat + f3naXa
                f3naXbt = f3naXbt + f3naXb
                f3naXct = f3naXct + f3naXc

                dQ = s%atom(ialpha)%shell(isorp)%dQ
                do inu = 1, norb_nu
                  do imu= 1, norb_mu
                    pvna_neighbors%aDblock(:,imu,inu) =                            &
     &                pvna_neighbors%aDblock(:,imu,inu)                            &
     &                + dQ*P_eq2*(stinky*f3naXat(:,imu,inu)                        &
     &                            - dstnA(:)*bcnax(imu,inu) + dstnA(:)*emnpl(imu,inu)&
     &                            - (1.0d0 - stinky)*DemnplA(:,imu,inu))

                    pvna_neighbors%bDblock(:,imu,inu) =                            &
     &                pvna_neighbors%bDblock(:,imu,inu)                            &
     &                + dQ*P_eq2*(stinky*f3naXbt(:,imu,inu)                        &
     &                            - dstnB(:)*bcnax(imu,inu) + dstnB(:)*emnpl(imu,inu)&
     &                            - (1.0d0 - stinky)*DemnplB(:,imu,inu))

                    pvna_neighbors%cDblock(:,imu,inu) =                            &
     &                pvna_neighbors%cDblock(:,imu,inu)                            &
     &                + dQ*P_eq2*(stinky*f3naXct(:,imu,inu)                        &
     &                            - dstnC(:)*bcnax(imu,inu) + dstnC(:)*emnpl(imu,inu)&
     &                            - (1.0d0 - stinky)*DemnplC(:,imu,inu))
                  end do !imu= 1, norb_mu
                end do !inu = 1, norb_nu

              end do !isorp = 1, species(in3)%nssh

              deallocate (bcnam, bcnax)
              deallocate (dpbcnam, dxbcnam, dybcnam)
              deallocate (f3naMa, f3naMb)
              deallocate (f3naXa, f3naXb, f3naXc)
              deallocate (f3naXat, f3naXbt, f3naXct)
              deallocate (sterm, dterm, emnpl)
              deallocate (DemnplA, DemnplB, DemnplC)
            end if !if (mneigh .ne. 0)
          end do !ineigh
        end do !iatom

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vna_3c

! End Module
! ===========================================================================
        end module M_Dassemble_3c
