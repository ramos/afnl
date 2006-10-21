!
! A MODULE for numerical Integration and to solve ODE's 
!
! Copyright (C) 2003  Alberto Ramos <alberto@martin.ft.uam.es>
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
MODULE NumTypes
! *
! ***************************************************
! *
! * Define the numerical Types Single precision, 
! * Double Precision, and Complex double precision
! * to make the code portable.
! *
! ***************************************************

  ! Here we define our Single Precision (SP)
  ! and Double Precision (DP) data types, for
  ! both Real and Complex variables.

  Integer, Parameter :: SP = Kind(1.0)
  Integer, Parameter :: DP = Kind(1.D0)

  Integer, Parameter :: SPC = Kind( (1.0_SP, 1.0_SP) )
  Integer, Parameter :: DPC = Kind( (1.0_DP, 1.0_DP) )

End MODULE NumTypes

