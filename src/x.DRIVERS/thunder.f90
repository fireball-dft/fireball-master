! copyright info:
!
!                             @Copyright 2011
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

! thunder.f90
! Program Description
! ===========================================================================
!      Driver to run FIREBALL version which includes BEGIN and CREATE.
!
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        program thunder

! /SYSTEM
        use M_species
        use M_atom_functions
        use M_atomPP_functions

! /BEGIN
        use M_rcatms

        implicit none

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer istructure, iseparate

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
!       real timei, timef

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call cpu_time (time_begin)
        logfile = 1
        iseparate = 1
        open (unit = logfile, file = 'output.log', status = 'replace')
        call welcome

! ===========================================================================
! ---------------------------------------------------------------------------
!             R E A D   I N   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
! Call read_input to define variables from the begin.input file
        call read_Fdata_location
        allocate (species_PP (nspecies))
        call read_begin
        call write_species

! ===========================================================================
! ---------------------------------------------------------------------------
!               G E N E R A T E    W A V E F U N C T I O N S
!                      A N D    P O T E N T I A L S
! ---------------------------------------------------------------------------
! ===========================================================================
        call rcatms

        call cpu_time (time_end)
        logfile = 1
        write (logfile,*) ' THUNDER RUNTIME : ', time_end - time_begin, '[sec] '
        write (*,*) ' THUNDER RUNTIME : ', time_end - time_begin, '[sec] '
        close (logfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Working on structure - ', a25)

! End Program
! ===========================================================================
        stop
        end program thunder
