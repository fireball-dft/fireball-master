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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! M_Dassemble_rho_McWEDA
! Module Description
! ===========================================================================
!>       This is a module containing  the assembler programs required
!! to calculate the matrix elements for the densities rho_in and rho_at
!! used in the McWEDA-Harris interactions (e.g. Sankey-Niklewski).
!! It contains the following subroutines within the module:
!!
!!       Dassemble_rho_2c.f90 - assembles  two center part for rho_in and rho_at
!!       Dassemble_rho_3c.f90 - three center part for rho_in
!!       Dassemble_rho_average.f90 - calculates the final result for average
!!                                  densities
!!       Dassemble_rho_weighted_2c.f90 - assembles two center part
!!                                      for Wrho, Wrho_bond
!!       Dassemble_rho_weighted_3c.f90 - three center part for Wrho
!!       Dassemble_S_weighted.f90 - assembles overlap_weighted (and averaged)
!!                                 PRB 71, 235101 (2005):
!!                                 denominators in Eqs. (19), (22) and (25)
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
         module M_Dassemble_rho_McWEDA
         use M_assemble_blocks
         use M_configuraciones
         use M_Fdata_2c
         use M_Fdata_3c
         use M_rotations
         use M_Drotations 
        
! Type Declaration
! ===========================================================================
! None

! module procedures
         contains

! ===========================================================================
! Dassemble_rho_2c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates neutral atom 2-center input density matrix
!! interactions (for rho_in and rho_bond)
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
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_rho_2c (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer in1, in2, in3, inu, imu !< species numbers
        integer jatom                   !< neighbor of iatom
        integer interaction, isubtype   !< which interaction and subtype
        integer num_neigh               !< number of neighbors
        integer matom                   !< matom is the self-interaction atom
        integer mbeta                   !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu        !< size of the block for the pair

        real Qneutral                   ! charge
        real z                          !< distance between r1 and r2

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of density matrix in molecular coordinates
! vdbcxcx = vectorized derivative of density matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx
        real, dimension (:, :), allocatable :: dbcxcm
        real, dimension (:, :, :), allocatable :: vdbcxcm
        real, dimension (:, :, :), allocatable :: vdbcxcx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: prho_in_neighbors
        type(T_assemble_neighbors), pointer :: prho_in
        type(T_assemble_block), pointer :: prho_bond_neighbors
        type(T_assemble_neighbors), pointer :: prho_bond

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          prho_in=>s%rho_in(iatom)
          prho_bond=>s%rho_bond(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            prho_in_neighbors=>prho_in%neighbors(ineigh)
            prho_bond_neighbors=>prho_bond%neighbors(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (prho_in_neighbors%Dblock(3,norb_mu, norb_nu))
            allocate (prho_bond_neighbors%Dblock(3,norb_mu, norb_nu))
            prho_in_neighbors%Dblock = 0.0d0
            prho_bond_neighbors%Dblock = 0.0d0

! SET-UP STUFF
! ***************************************************************************
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
            call epsilon_function (r2, sighat, eps)
            call Depsilon_2c (r1, r2, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! CALL GetDMES AND GET rho_in FOR ONTOP CASE (OFF-SITE matrix elements)
! ***************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then
               
! Do nothing here - special case. Interaction already calculated in atm case.

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in bcxcx, where x means crytal
! coordinates.

! FORCES - ONTOP LEFT CASE
! ****************************************************************************
! For the rho_in_ontopL case, the potential is in the first atom - left (iatom):
! <mu|(rho_mu + rho_nu)|nu> -> (left) <mu|(rho_mu)|nu>
! dbcxcm is the "scalar" derivative of the matrix; vdbcxcm is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in cry
              interaction = P_rho_ontopL
              in3 = in2

! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of denstiy matrix in molecular coordinates
! vdbcxcx = vectorized derivative of density matrix in crystal coordinates
              allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (bcxcx (norb_mu, norb_nu)); bcxcx = 0.0d0
              allocate (dbcxcm (norb_mu, norb_nu)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, norb_mu, norb_nu)); vdbcxcm = 0.0d0
              allocate (vdbcxcx (3, norb_mu, norb_nu)); vdbcxcx = 0.0d0
             
              do isubtype = 1, species(in1)%nssh
                Qneutral = species(in1)%shell(isubtype)%Qneutral
                
                call getDMEs_Fdata_2c (in1, in3, interaction, isubtype, z,   &
     &                                 norb_mu, norb_nu, bcxcm, dbcxcm)
                            
! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do
                
                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,  &
     &                        vdbcxcm, vdbcxcx)

                prho_in_neighbors%Dblock =                                   &
     &            prho_in_neighbors%Dblock + vdbcxcx*Qneutral
                prho_bond_neighbors%Dblock =                                 &
     &            prho_bond_neighbors%Dblock + vdbcxcx*Qneutral
              
              end do

! For the rho_in_ontopR case, the potential is in the second atom
! - right (iatom): <mu|(rho_mu + rho_nu)|nu> -> (right) <mu|(rho_nu)|nu>
              interaction = P_rho_ontopR
              in3 = in2
              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,   &
     &                                 norb_mu, norb_nu, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,  &
     &                        vdbcxcm, vdbcxcx)

                prho_in_neighbors%Dblock =                                   &
     &            prho_in_neighbors%Dblock + vdbcxcx*Qneutral
                prho_bond_neighbors%Dblock =                                 &
     &            prho_bond_neighbors%Dblock + vdbcxcx*Qneutral
              end do
              deallocate (bcxcm, bcxcx, dbcxcm, vdbcxcm, vdbcxcx)
            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms

! CALL DOSCENTROS AND GET rho_in FOR ATOM CASE
! ***************************************************************************
! The rho_in two-center terms are: ontop (L), ontop (R), and atom.
! First, do rho_in_atom case. Here we compute <i | v(j) | i> matrix elements.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          prho_in=>s%rho_in(iatom)
          prho_bond=>s%rho_bond(iatom)
          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            prho_in_neighbors=>prho_in%neighbors(matom)
            prho_bond_neighbors=>prho_bond%neighbors(matom)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! SET-UP STUFF
! *************************************************************************
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
            call epsilon_function (r2, sighat, eps)
            call Depsilon_2c (r1, r2, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

            if (iatom .eq. jatom .and. mbeta .eq. 0) then
! one center case : calculate both rho_in and rho_bond
              interaction = P_rho_atom
              in3 = in1
              
! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of denstiy matrix in molecular coordinates
! vdbcxcx = vectorized derivative of density matrix in crystal coordinates
              norb_nu = species(in3)%norb_max
              allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (bcxcx (norb_mu, norb_nu)); bcxcx = 0.0d0
              allocate (dbcxcm (norb_mu, norb_nu)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, norb_mu, norb_nu)); vdbcxcm = 0.0d0
              allocate (vdbcxcx (3, norb_mu, norb_nu)); vdbcxcx = 0.0d0

              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,   &
     &                                 norb_mu, norb_nu, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,  &
     &                        vdbcxcm, vdbcxcx)

                prho_in_neighbors%Dblock =                                   &
     &            prho_in_neighbors%Dblock + vdbcxcx*Qneutral
                prho_bond_neighbors%Dblock =                                 &
     &            prho_bond_neighbors%Dblock + vdbcxcx*Qneutral
              end do
              deallocate (bcxcm, bcxcx, dbcxcm, vdbcxcm, vdbcxcx)
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_rho_atom
              in3 = in1
                
! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of denstiy matrix in molecular coordinates
! vdbcxcx = vectorized derivative of density matrix in crystal coordinates
              norb_nu = species(in3)%norb_max
              allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (bcxcx (norb_mu, norb_nu)); bcxcx = 0.0d0
              allocate (dbcxcm (norb_mu, norb_nu)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, norb_mu, norb_nu)); vdbcxcm = 0.0d0
              allocate (vdbcxcx (3, norb_mu, norb_nu)); vdbcxcx = 0.0d0

              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral

                call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,   &
     &                                 norb_mu, norb_nu, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,  &
     &                        vdbcxcm, vdbcxcx)

                prho_in_neighbors%Dblock =                                   &
     &            prho_in_neighbors%Dblock + vdbcxcx*Qneutral
              end do
              deallocate (bcxcm, bcxcx, dbcxcm, vdbcxcm, vdbcxcx)
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_rho_2c


! ===========================================================================
! Dassemble_rho_3c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates neutral atom 3-center matrix interactions
!! for rho_in - used to evaluate exchange-correlation interactions.
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
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_rho_3c (s)
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
        integer isubtype                 !< which subtype
!        integer interaction, isorp       !< which interaction and subtype
        integer ix, iindex

        integer norb_mu, norb_nu         !< size of the block for the pair
        integer imu, inu
        real Qneutral                    !< charge
        real z                           !< distances between r1 and r2
        real x, cost                     !< dnabc and angle
 
        real, dimension (3, 3) :: eps     !< the epsilon matrix
!        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2, r3, r12, r21   !< positions
        real, dimension (3) :: sighat     !< unit vector along r2 - r1
        real, dimension (3) :: rhat       !< unit vector along bc - r3
        real, dimension (3) :: amt
        real, dimension (3) :: bmt
        
        real, dimension (3, 3, 3) :: depsA  !< the Depsilon matrix for the bond-charge
        real, dimension (3, 3, 3) :: depsB  !< the Depsilon matrix for the potential


! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! d..bcxcm = derivative of density matrix in molecular coordinates
! vdxcM.. = vectorized derivative of density matrix in molecular coordinates
! vdxcX.. = vectorized derivative of density matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx
        real, dimension (:, :), allocatable :: dpbcxcm
        real, dimension (:, :), allocatable :: dxbcxcm
        real, dimension (:, :), allocatable :: dybcxcm
        
        real, dimension (:, :, :), allocatable :: vdxcMa
        real, dimension (:, :, :), allocatable :: vdxcMb
        real, dimension (:, :, :), allocatable :: vdxcXa
        real, dimension (:, :, :), allocatable :: vdxcXb
        real, dimension (:, :, :), allocatable :: vdxcXc
        
        type(T_Fdata_cell_3c), pointer :: pFdata_cell
        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface
        type(T_assemble_neighbors), pointer :: prho_in
        type(T_assemble_block), pointer :: prho_in_neighbors

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
              
              prho_in=>s%rho_in(iatom)
              prho_in_neighbors=>prho_in%neighbors(mneigh)
              
              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max
              
! Allocade DFdata blocks
              allocate (prho_in_neighbors%aDblock(3, norb_mu, norb_nu)); prho_in_neighbors%aDblock = 0.00
              allocate (prho_in_neighbors%bDblock(3, norb_mu, norb_nu)); prho_in_neighbors%bDblock = 0.00
              allocate (prho_in_neighbors%cDblock(3, norb_mu, norb_nu)); prho_in_neighbors%cDblock = 0.00

! SET-UP STUFF
! ***************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
              r21 = r2 - r1
              
              z = distance (r1, r2)

              ! unit vector in sigma direction.
              if (z .lt. 1.0d-05) then
                sighat(1) = 0.0d0
                sighat(2) = 0.0d0
                sighat(3) = 1.0d0
              else
                sighat = (r2 - r1)/z
              end if

! ***************************************************************************
! Find rnabc = vector pointing from center of bondcharge to r3
! This gives us the distance dnabc (or x value in the 2D grid).
              r12 = 0.5d0*(r1 + r2)
              x = distance (r12, r3)

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = (r3 - 0.5d0*(r1 + r2))/x
              end if
              cost = dot_product(sighat, rhat)
              if (abs(cost) - 1.0d0 .lt. 0.10d0) then
                cycle
              end if               
              call epsilon_function (rhat, sighat, eps)
              call Depsilon_3c (r1, r2, r21, z, r3, rhat, eps, depsA, depsB)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates, stored in bcxcm. where m means molecular coordinates.
! Rotate the matrix into crystal coordinates. The rotated  matrix elements
! are stored in bcxcx, where x means crytal coordinates.
              ! matrix element derivatives
              allocate (bcxcm(norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (bcxcx(norb_mu, norb_nu)); bcxcx = 0.0d0
              allocate (dpbcxcm(norb_mu, norb_nu)); dpbcxcm = 0.0d0
              allocate (dxbcxcm(norb_mu, norb_nu)); dxbcxcm = 0.0d0
              allocate (dybcxcm(norb_mu, norb_nu)); dybcxcm = 0.0d0
              
              ! vectorial representations
              allocate (vdxcMa(3, norb_mu, norb_nu)); vdxcMa = 0.0d0
              allocate (vdxcMb(3, norb_mu, norb_nu)); vdxcMb = 0.0d0
              allocate (vdxcXa(3, norb_mu, norb_nu)); vdxcXa = 0.0d0
              allocate (vdxcXb(3, norb_mu, norb_nu)); vdxcXb = 0.0d0
              allocate (vdxcXc(3, norb_mu, norb_nu)); vdxcXc = 0.0d0
              
              do isubtype = 1, species(in3)%nssh
                Qneutral = species(in3)%shell(isubtype)%Qneutral
                
                call getDMEs_Fdata_3c (in1, in2, in3, P_rho_3c, isubtype, x, &
     &                                 z, norb_mu, norb_nu, cost, rhat,      &
     &                                 sighat, bcxcm, dpbcxcm, dxbcxcm, dybcxcm)

                ! Rotate into crystal coordinates
                call rotate (in1, in2, eps, norb_mu, norb_nu, bcxcm, bcxcx)
                do ix = 1, 3

! The first piece will be the force with respect to atom 3.
                 if (x .gt. 1.0d-5) then
                  amt(ix) = (sighat(ix) - cost*rhat(ix))/x
                 else
                  amt = 0.0d0
                 end if

                 pFdata_bundle => Fdata_bundle_3c(in1, in2, in3)
                 pFdata_cell =>                                               &
     &             pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(P_rho_3c,isubtype,1))

                 do iindex = 1, pFdata_cell%nME
                   imu = pFdata_cell%mu_3c(iindex)
                   inu = pFdata_cell%nu_3c(iindex)

! Now recover f3naMa which is a two-dimensional array
                   vdxcMa(ix,imu,inu) = rhat(ix)*dxbcxcm(imu,inu)            &
     &                                 + amt(ix)*dpbcxcm(imu,inu)

! The second piece will be the force with respect to atom 1.
                   bmt(ix) = (cost*sighat(ix) - rhat(ix))/z

! Now recover f3naMb which is a two-dimensional array
                   vdxcMb(ix,imu,inu) = - sighat(ix)*dybcxcm(imu, inu)       &
     &                                   + bmt(ix)*dpbcxcm(imu, inu)         &
     &                                   - vdxcMa(ix,imu,inu)/2.0d0
                 end do ! iindex
                end do ! ix

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
                call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu, bcxcm, &
     &                        vdxcMa, vdxcXa)

! Force on the neutral atom with respect to atom 1 (f3naMb).
                call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu, bcxcm, &
     &                        vdxcMb, vdxcXb)

! Make things force-like and determine f3naXc, whcih is found from Newtons Laws:
                vdxcXa(:,:,:) = - vdxcXa(:,:,:)
                vdxcXb(:,:,:) = - vdxcXb(:,:,:)
                vdxcXc(:,:,:) = - vdxcXa(:,:,:) - vdxcXb(:,:,:)

                prho_in_neighbors%aDblock = prho_in_neighbors%aDblock        &
     &                                     + vdxcXa*Qneutral
                prho_in_neighbors%bDblock = prho_in_neighbors%bDblock        &
     &                                     + vdxcXb*Qneutral
                prho_in_neighbors%cDblock = prho_in_neighbors%cDblock        &
     &                                     + vdxcXc*Qneutral
              end do ! isubtype = 1, species(in3)%nssh

              deallocate (bcxcm, bcxcx, dpbcxcm, dxbcxcm, dybcxcm)
              deallocate (vdxcMa, vdxcMb, vdxcXa, vdxcXb, vdxcXc)
            end if ! if (mneigh .ne. 0)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_rho_3c


! ===========================================================================
! Dassemble_rho_weighted_2c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine assembles  two center part for rho_in_weighted and
!! rho_bond_weighted.
!!       rho_in_weighted: numerator in Eq. (19): PRB 71, 235101 (2005)
!!       rho_bond_weighted : numerator in Eqs. (22), (25):PRB 71, 235101 (2005)
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
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_rho_weighted_2c (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh            !< counter over atoms and neighbors
        integer in1, in2, in3, inu, imu  !< species numbers
        integer jatom                    !< neighbor of iatom
        integer interaction, isubtype    !< which interaction and subtype
        integer num_neigh                !< number of neighbors
        integer matom                    !< matom is the self-interaction atom
        integer mbeta                    !< cell containing neighbor of iatom
          
        integer norb_mu, norb_nu        !< size of the block for the pair
        integer nssh_i, nssh_j         !< size of the block for the pair

        real z                           !< distance between r1 and r2
        real Qneutral
        
        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sighat     !< unit vector along r2 - r1  

! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of density matrix in molecular coordinates
! vdbcxcx = vectorized derivative of density matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: dbcxcm
        real, dimension (:, :, :), allocatable :: vdbcxcm

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pWrho_in_neighbors
        type(T_assemble_neighbors), pointer :: prho_in_weighted
        type(T_assemble_block), pointer :: pWrho_bond_neighbors
        type(T_assemble_neighbors), pointer :: prho_bond_weighted

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          prho_in_weighted=>s%rho_in_weighted(iatom)
          prho_bond_weighted=>s%rho_bond_weighted(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          nssh_i = species(in1)%nssh
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pWrho_in_neighbors=>prho_in_weighted%neighbors(ineigh)
            pWrho_bond_neighbors=>prho_bond_weighted%neighbors(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
! Allocate block size
            nssh_j = species(in2)%nssh
            allocate (pWrho_in_neighbors%Dblock(3, nssh_i, nssh_j))
            allocate (pWrho_bond_neighbors%Dblock(3, nssh_i, nssh_j))
            pWrho_in_neighbors%Dblock = 0.0d0
            pWrho_bond_neighbors%Dblock = 0.0d0

! Calculate the distance between the two centers.
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! CALL DOSCENTROS AND GET rho_in_weighted FOR ONTOP CASE
! (i.e. OFF-SITE matrix elements)
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm).
! No rotation requiered for this case
! (weights are spherical, see definition Eq. (18) McWEDA paper)
!
! For the rho_in_weighted_ontopL case, the potential is in the first atom -
! left (iatom): <mu|(rho_mu + rho_nu)|nu> -> (left) <mu|(rho_mu)|nu>
              interaction = P_rhoS_ontopL
              in3 = in2

! bcxcm = density matrix in molecular coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of denstiy matrix in molecular coordinates
              allocate (bcxcm (nssh_i, nssh_j)); bcxcm = 0.0d0
              allocate (dbcxcm (nssh_i, nssh_j)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, nssh_i, nssh_j)); vdbcxcm = 0.0d0

              do isubtype = 1, species(in1)%nssh
                Qneutral = species(in1)%shell(isubtype)%Qneutral
                call getDMEs_Fdata_2c (in1, in3, interaction, isubtype, z,   &
     &                                 nssh_i, nssh_j, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                
                do inu = 1, nssh_j !norb_nu
                  do imu = 1, nssh_i !norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

                pWrho_in_neighbors%Dblock =                                   &
     &            pWrho_in_neighbors%Dblock + vdbcxcm*Qneutral
                pWrho_bond_neighbors%Dblock =                                 &
     &            pWrho_bond_neighbors%Dblock + vdbcxcm*Qneutral
              end do

! For the rho_in_weighted_ontopR case, the potential is in the second atom -
! right (iatom): <mu|(rho_mu + rho_nu)|nu> -> (right) <mu|(rho_nu)|nu>
              interaction = P_rhoS_ontopR
              in3 = in2

              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,   &
     &                                 nssh_i, nssh_j, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, nssh_j !norb_nu
                  do imu = 1, nssh_i !norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

                pWrho_in_neighbors%Dblock =                                  &
     &            pWrho_in_neighbors%Dblock + vdbcxcm*Qneutral
                pWrho_bond_neighbors%Dblock =                                &
     &            pWrho_bond_neighbors%Dblock + vdbcxcm*Qneutral
              end do
              deallocate (bcxcm, dbcxcm, vdbcxcm)

            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms

! CALL DOSCENTROS AND GET rho_in_weighted FOR ATOM CASE
! ****************************************************************************
! The rho_in_weighted two-center terms are: ontop (L), ontop (R), and atom.
! First, do rho_in_weighted_atom case.
! Here we compute <i|v(j)|i> matrix elements.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          prho_in_weighted=>s%rho_in_weighted(iatom)
          prho_bond_weighted=>s%rho_bond_weighted(iatom)
          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          nssh_i = species(in1)%nssh
          num_neigh = s%neighbors(iatom)%neighn

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pWrho_in_neighbors=>prho_in_weighted%neighbors(matom)
            pWrho_bond_neighbors=>prho_bond_weighted%neighbors(matom)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Calculate the distance between the two centers.
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

            if (iatom .eq. jatom .and. mbeta .eq. 0) then
! one center case
              interaction = P_rhoS_atom
              in3 = in1

! bcxcm = density matrix in molecular coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of denstiy matrix in molecular coordinates
              allocate (bcxcm (nssh_i, nssh_i)); bcxcm = 0.0d0
              allocate (dbcxcm (nssh_i, nssh_i)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, nssh_i, nssh_i)); vdbcxcm = 0.0d0

              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,   &
     &                                 nssh_i, nssh_i, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

                pWrho_in_neighbors%Dblock =                                  &
     &            pWrho_in_neighbors%Dblock + vdbcxcm*Qneutral
                pWrho_bond_neighbors%Dblock =                                &
     &            pWrho_bond_neighbors%Dblock + vdbcxcm*Qneutral
              end do
              deallocate (bcxcm, dbcxcm, vdbcxcm)
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). No rotations here
              interaction = P_rhoS_atom
              in3 = in1

! bcxcm = density matrix in molecular coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of denstiy matrix in molecular coordinates
              allocate (bcxcm (nssh_i, nssh_i)); bcxcm = 0.0d0
              allocate (dbcxcm (nssh_i, nssh_i)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, nssh_i, nssh_i)); vdbcxcm = 0.0d0

              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,   &
     &                                 nssh_i, nssh_i, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, nssh_j !norb_nu
                  do imu = 1, nssh_i !norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

                pWrho_in_neighbors%Dblock =                                  &
     &            pWrho_in_neighbors%Dblock + vdbcxcm*Qneutral
              end do
              deallocate (bcxcm, dbcxcm, vdbcxcm)
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None - 

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_rho_weighted_2c


! ===========================================================================
! Dassemble_rho_weighted_3c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine assembles  three center part for rho_weighted
!!       rho_weighted : numerator in Eq. (19): PRB 71, 235101 (2005)
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
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_rho_weighted_3c (s)
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
        integer in1, in2, in3  !, inu, imu
        integer isubtype                 !< which subtype
!        integer ix, iindex 
!        integer norb_mu, norb_nu
        integer nssh_i, nssh_j           !< size of the block for the pair
        integer issh, jssh               !< counter over shells
        
        real z                           !< distances between r1 and r2
        real x, cost                     !< dnabc and angle
        real Qneutral

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2, r3, r12, r21     !< positions
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - r3
        real, dimension (3) :: amt
        real, dimension (3) :: bmt
      
        real, dimension (3, 3, 3) :: depsA  !< the Depsilon matrix for the bond-charge
        real, dimension (3, 3, 3) :: depsB  !< the Depsilon matrix for the potential

! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! d..bcxcm = derivative of density matrix in molecular coordinates
! vdxcM.. = vectorized derivative of density matrix in molecular coordinates
! vdxcX.. = vectorized derivative of density matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx
        real, dimension (:, :), allocatable :: dpbcxcm
        real, dimension (:, :), allocatable :: dxbcxcm
        real, dimension (:, :), allocatable :: dybcxcm
        
        real, dimension (:, :, :), allocatable :: vdxcMa
        real, dimension (:, :, :), allocatable :: vdxcMb
        real, dimension (:, :, :), allocatable :: vdxcMc
        
!        real, dimension (:, :, :), allocatable :: vdbcxcm
!        real, dimension (:, :, :), allocatable :: vdbcxcx
                
!        type(T_Fdata_cell_3c), pointer :: pFdata_cell
!        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface
        type(T_assemble_neighbors), pointer :: prho_in_weighted
        type(T_assemble_block), pointer :: prho_in_neighbors_weighted
! Allocate Arrays
! ===========================================================================

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
              nssh_i = species(in1)%nssh
              
              prho_in_weighted=>s%rho_in_weighted(iatom)
              prho_in_neighbors_weighted=>prho_in_weighted%neighbors(mneigh)
              
              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              nssh_j = species(in2)%nssh
! Allocade DFdata blocks
              allocate (prho_in_neighbors_weighted%aDblock(3, nssh_i, nssh_j))
              prho_in_neighbors_weighted%aDblock = 0.00
              allocate (prho_in_neighbors_weighted%bDblock(3, nssh_i, nssh_j))
              prho_in_neighbors_weighted%bDblock = 0.00
              allocate (prho_in_neighbors_weighted%cDblock(3, nssh_i, nssh_j))
              prho_in_neighbors_weighted%cDblock = 0.00              

! Allocate DFdata blocks
              allocate (prho_in_neighbors_weighted%aDblock(3, nssh_i, nssh_j)); prho_in_neighbors_weighted%aDblock = 0.00
              allocate (prho_in_neighbors_weighted%bDblock(3, nssh_i, nssh_j)); prho_in_neighbors_weighted%bDblock = 0.00
              allocate (prho_in_neighbors_weighted%cDblock(3, nssh_i, nssh_j)); prho_in_neighbors_weighted%cDblock = 0.00

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
              r21 = r2 - r1
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

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = (r3 - 0.5d0*(r1 + r2))/x
              end if

              cost = dot_product(sighat, rhat)
              if (abs(cost) - 1.0d0 .lt. 0.10d0) then
                cycle
              end if 
              
              call epsilon_function (rhat, sighat, eps)
              call Depsilon_3c (r1, r2, r21, z, r3, rhat, eps, depsA, depsB)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates, stored in bcxcm. where m means molecular coordinates.
! Rotate the matrix into crystal coordinates. The rotated  matrix elements
! are stored in bcxcx, where x means crytal coordinates.
              do isubtype = 1, species(in3)%nssh

                ! matrix element derivatives
                allocate (bcxcm (nssh_i, nssh_j)); bcxcm = 0.0d0
                allocate (bcxcx (nssh_i, nssh_j)); bcxcx = 0.0d0
                allocate (dpbcxcm (nssh_i, nssh_j)); dpbcxcm = 0.0d0
                allocate (dxbcxcm (nssh_i, nssh_j)); dxbcxcm = 0.0d0
                allocate (dybcxcm (nssh_i, nssh_j)); dybcxcm = 0.0d0

                ! vectorized derivatives
                allocate (vdxcMa (3, nssh_i, nssh_j)); vdxcMa = 0.0d0
                allocate (vdxcMb (3, nssh_i, nssh_j)); vdxcMb = 0.0d0
                allocate (vdxcMc (3, nssh_i, nssh_j)); vdxcMc = 0.0d0
                
                Qneutral = species(in3)%shell(isubtype)%Qneutral
                
               
                call getDMEs_Fdata_3c (in1, in2, in3, P_rhoS_3c, isubtype, x,&
     &                                 z, nssh_i, nssh_j, cost, rhat, sighat,&
     &                                 bcxcm, dpbcxcm, dxbcxcm, dybcxcm)
              
! The first piece will be the force with respect to atom 3.
                if (x .gt. 1.0d-5) then
                  amt(:) = (sighat(:) - cost*rhat(:))/x
                else
                  amt = 0.0d0
                end if
                do issh = 1, species(in1)%nssh
                   do jssh = 1, species(in2)%nssh
                     vdxcMa(:,issh,jssh) = rhat(:)*dxbcxcm(issh,jssh) + amt(:)*dpbcxcm(issh,jssh)
                     bmt(:) = (cost*sighat(:) - rhat(:))/z 
                     
                     vdxcMb(:,issh,jssh) = - sighat(:)*dybcxcm(issh,jssh)                             &
     &                                      + bmt(:)*dpbcxcm(issh,jssh) - vdxcMa(:,issh,jssh)/2.0d0
                     vdxcMc(:,issh,jssh) = - vdxcMa(:,issh,jssh) - vdxcMb(:,issh,jssh)
                   end do ! jssh
                end do ! issh
                 
                prho_in_neighbors_weighted%aDblock =                         &
     &            prho_in_neighbors_weighted%aDblock + vdxcMa*Qneutral
                prho_in_neighbors_weighted%bDblock =                         &
     &            prho_in_neighbors_weighted%bDblock + vdxcMb*Qneutral
                prho_in_neighbors_weighted%cDblock =                         &
     &            prho_in_neighbors_weighted%cDblock + vdxcMc*Qneutral
                  
              deallocate (bcxcm, bcxcx, dpbcxcm, dxbcxcm, dybcxcm)
              deallocate (vdxcMa, vdxcMb, vdxcMc)
              end do ! isubtype = 1, species(in3)%nssh
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_rho_weighted_3c


! ===========================================================================
! destroy_Dassemble_rho
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the rho_McWEDA
!! information.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_Dassemble_rho (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                             !< counter over atoms
        integer ineigh 
        integer mneigh                           !< counter over neighbors

! Procedure
! ===========================================================================
! three-center interactions
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%ncommon
            mneigh = s%neighbors(iatom)%neigh_common(ineigh)
            deallocate (s%rho_in(iatom)%neighbors(mneigh)%aDblock)
            deallocate (s%rho_in(iatom)%neighbors(mneigh)%bDblock)
            deallocate (s%rho_in(iatom)%neighbors(mneigh)%cDblock)
            deallocate (s%rho_in_weighted(iatom)%neighbors(mneigh)%aDblock)
            deallocate (s%rho_in_weighted(iatom)%neighbors(mneigh)%bDblock)
            deallocate (s%rho_in_weighted(iatom)%neighbors(mneigh)%cDblock)
          end do  
        end do

! two-center interactions
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%rho_in_weighted(iatom)%neighbors(ineigh)%block)
            deallocate (s%rho_bond_weighted(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%rho_in_weighted(iatom)%neighbors)
          deallocate (s%rho_bond_weighted(iatom)%neighbors)
        end do
        deallocate (s%rho_in_weighted)
        deallocate (s%rho_bond_weighted)

        do iatom = 1, s%natoms
          do ineigh=1, s%neighbors(iatom)%neighn
            deallocate (s%rho_in(iatom)%neighbors(ineigh)%block)
            deallocate (s%rho_bond(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%rho_in(iatom)%neighbors)
          deallocate (s%rho_bond(iatom)%neighbors)
        end do
        deallocate (s%rho_in)
        deallocate (s%rho_bond)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_Dassemble_rho

! End Module
! ===========================================================================
        end module M_Dassemble_rho_McWEDA
