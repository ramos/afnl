!
! A MODULE for MIXMAX Random Number Generator
!
! Copyright (C) 2005  Alberto Ramos <alberto.ramos@desy.de>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA
! 

! ***************************************************
! *
MODULE MixMax_tlb
! *
! ***************************************************

  USE ISO_FORTRAN_ENV, Only : ERROR_UNIT


  IMPLICIT NONE

  Integer (kind=8) :: skipmat256(2,2) = (/ (/1,1/), (/2,2/)/)




End MODULE MixMax_tlb
