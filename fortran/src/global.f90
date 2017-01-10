!    SNEASY:  Simple Nonlinear EArth SYstem model, a composite of
!               DOECLIM, Carbon Cycle, and MOC Boxmodel
!
!    Copyright (C) 2009  Klaus Keller, Nathan Urban, and Brian Tuttle
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!    Klaus Keller, klaus@psu.edu
!    Nathan Urban, nurban@psu.edu
!    Brian Tuttle, btuttle@psu.edu
!
!===========================================================================

MODULE global

    implicit none

    integer, parameter :: i4b = Selected_Int_Kind(9)
    integer, parameter :: i8b = Selected_Int_Kind(18)
    integer, parameter :: SP = Selected_Real_Kind(6,37)
    integer, parameter :: DP = Selected_Real_Kind(15,307)
    
    real(DP), parameter :: Pi = 3.141592653589793d0

    integer(i4b) :: nsteps
    real(DP) :: deltat

! Random seed variables:
    integer(i4b) :: rand_seed
    logical :: seeded

! MPI variables
    integer, parameter :: MASTER = 0
    integer(i4b) :: ierr, n_procs, iam

CONTAINS

!---------------------------------------------------------------------------
FUNCTION stringcat(string1, string2) RESULT(catstrings)
!  ==========================================================================
! |  This function does a simple concatenation of two strings.  It is        |
! |  included here because the C preprocessor treats Fortran's concatenation |
! |  operator as a comment and truncates the line at that point.  This is    |
! |  a way to work around that behavior.                                     |
!  ==========================================================================

    character(len=*) :: string1, string2
    character(len=len(string1)+len(string2)) :: catstrings

    catstrings = string1 // string2

END FUNCTION stringcat
!---------------------------------------------------------------------------
END MODULE global
