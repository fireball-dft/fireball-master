! copyright info:
!
!                             @Copyright 2007
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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
! M_case4.f90
! Program Description
! ============================================================================
!      This is a module calculating the integrals of two centers for the
!
! ============================================================================
! Code written by:
! Hong Wang
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ============================================================================
! Module Declaration
! ============================================================================
        module M_vnl
        use M_species
        use M_atom_functions
        use M_atomPP_functions
        use M_integrals_2c

        implicit none

! Type Declaration
! ============================================================================
! two-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_2C, containing
! all Fdata for a particular interaction/subinteraction type.
! ===========================================================================

! module procedures
        contains

! ===========================================================================
! initialize_vnl
! ===========================================================================
! Program Description
! ===========================================================================
!       We need to determine how many interactions belong to each nspecies
! bundle pair. This routine just counts how many total interactions contribute
! to that bundle.  Something like overlap is obviously only 1 interaction
! added, but something like vna needs number of interactions based on the
! number of shells.
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
!
!and
!
! Barry Haycock
! Department of Physics,
! Dublin Institute of Technology,
! Dublin 2.
! +353 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_vnl
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_vnl


! ============================================================================
! vnl
! ============================================================================
! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the site of one of the orbitals.  The potential V(1) is something like Vxc
! for the exchange correlation potential, Vna for the neutral atom potential,
! or 1 for the overlap term.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
!
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960

! Subroutine Declaration
! ===========================================================================
        subroutine vnl_2c
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        real dmax                           !< max distance between two centers
        real drr                            !< distance between mesh points
        real d                              !< distance between the two centers
        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 30) filename
        character (len = 25) interactions

        logical skip

! Procedure
! ============================================================================
        write (*,*)
        write (*,*) ' ******************************************************* '
        write (*,*) '              P P   I N T E R A C T I O N S              '
        write (*,*) ' ******************************************************* '

! Assign values to the unrequired variables for this specific interaction.
        isorp = 0
        ideriv = 999

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            call make_munu_PP (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

            ! Open ouput file for this species pair
            write (filename, '("/vnl.",i2.2,".",i2.2,".dat")')               &
     &             species(ispecies)%nZ, species(jspecies)%nZ
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

! Set up grid loop control constants
            dmax = species(ispecies)%rcutoffA_max + species(jspecies)%rcutoff_PP
            drr = dmax/ float(ndd_vnl-1)
            d = -drr
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoff_PP

            ! Set integration limits
            rhomin = 0.0d0
            rhomax = min(rcutoff1, rcutoff2)

            ! open directory file
            write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')             &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &            status = 'unknown', position = 'append')
            write (13,100) pFdata_bundle%nFdata_cell_2c, P_vnl, isorp,   &
     &                     filename(2:30), pFdata_cell%nME, ndd_vnl, dmax
            close (unit = 13)

! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2,"_munu_2c.",i2.2,".",i2.2,".dat")')  &
     &             P_vnl, species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c), index_2c = 1, nME2c_max)

! Loop over grid
            write (*,200) species(ispecies)%nZ, species(jspecies)%nZ
            do igrid = 1, ndd_vnl
              d = d + drr

              ! Set integration limits
              zmin = max(-rcutoff1, d - rcutoff2)
              zmax = min(rcutoff1, d + rcutoff2)

              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_vnl, nrho_vnl, rint_vnl,      &
     &                                   phifactor, zmin, zmax, rhomin,      &
     &                                   rhomax, pFdata_cell%fofx)

              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c), index_2c = 1, nME2c_max)
            end do ! igrid
            write (11,*)
          end do ! jpecies
        end do ! ispecies

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vnl integrals for nZ = ', i3,               &
     &              ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vnl_2c


! ============================================================================
! make_munu_PP
! ============================================================================
! Subroutine Description
! ============================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)):
!
! num_orb (in1): number of orbitals in atom-type in1
! mu (index,in1,in2): the mu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2
! nu (index,in1,in2): the nu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2

! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
!
! For the "atom" case, the number of non-zero matrix elements
! is dependent only on one atom. Both wavefunctions are located at
! the same site. This does that by "acknowledging" the jspecies call,
! But comopletely ignoring it, except for addressing the N_Nu L_Nu M_Nu's
!
! Subroutine Declaration
! ============================================================================
        subroutine make_munu_PP (itype, ispecies, jspecies)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Auguments Declaration and Description
! ============================================================================
! None

! Parameters and Data Declaration
! ============================================================================
! None

! Input
! ============================================================================

! Local Variable Declaration adn Description
! ============================================================================
        integer, intent(in) :: itype, ispecies, jspecies

        integer index_2c                ! counter for matrix location - mu, nu
        integer issh
        integer jssh                    ! index for looping over shells
        integer mvalue
        integer n1, l1, n2, l2

        integer nME2c_max

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).
        ! cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        nME2c_max = 0

! First, find the maximum number of matrix elements for each species pair.
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          do jssh = 1, species(jspecies)%nssh_PP
            l2 = species(jspecies)%shell_PP(jssh)%lssh
            do mvalue = -min(l1,l2), min(l1,l2)
              nME2c_max = nME2c_max + 1
            end do
          end do
        end do

! Now allocate the sizes for mu_2c, nu_2c, and the quantum numbers NLM
! for each mu and nu pair.
! Overlap Interactions:
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
        pFdata_cell%nME = nME2c_max

        allocate (pFdata_cell%mu_2c(nME2c_max))
        allocate (pFdata_cell%nu_2c(nME2c_max))
        allocate (pFdata_cell%mvalue_2c(nME2c_max))

        allocate (pFdata_cell%N_mu(nME2c_max))
        allocate (pFdata_cell%L_mu(nME2c_max))
        allocate (pFdata_cell%M_mu(nME2c_max))

        allocate (pFdata_cell%N_nu(nME2c_max))
        allocate (pFdata_cell%L_nu(nME2c_max))
        allocate (pFdata_cell%M_nu(nME2c_max))

! Set the values for NLM of each mu, nu pair.
        index_2c = 0
        n1 = 0
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh_PP
            l2 = species(jspecies)%shell_PP(jssh)%lssh
            n2 = n2 + l2 + 1
            do mvalue = -min(l1,l2), min(l1,l2)
              index_2c = index_2c + 1
              pFdata_cell%mu_2c(index_2c) = n1 + mvalue
              pFdata_cell%nu_2c(index_2c) = n2 + mvalue
              pFdata_cell%mvalue_2c(index_2c) = 0

              pFdata_cell%N_mu(index_2c) = issh
              pFdata_cell%L_mu(index_2c) = l1
              pFdata_cell%M_mu(index_2c) = mvalue

              pFdata_cell%N_nu(index_2c) = jssh
              pFdata_cell%L_nu(index_2c) = l2
              pFdata_cell%M_nu(index_2c) = mvalue

            end do
            n2 = n2 + l2
          end do
          n1 = n1 + l1
        end do

! End Subroutine
! =============================================================================
        return
        end subroutine make_munu_PP


! ============================================================================
! rint_vnl
! ============================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================
        function rint_vnl (itype, ispecies, jspecies, isorp, d, rho, z1, z2, &
     &                     ideriv, index_2c)
        implicit none

        real rint_vnl

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) :: d            !< distance between the two centers
        real, intent (in) :: rho
        real, intent (in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1, r2
        real rmax                         ! maximum value of r along grid

        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = isorp
        idummy = ideriv
        dummy = d

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

!Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

!Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)
        r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n1)%rcutoffA
        psi1val = psiofr (r1, rmax, ispecies, n1)

!find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        psi2val = vPPofr (r2, jspecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r2, z2, psi2val)
        vofr = 1.0d0

!Actual function (Ylm's are calculated above and multilplied after integration
        rint_vnl = psi1val*vofr*psi2val*rho

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_vnl


! End Module
! =============================================================================
        end module M_vnl
