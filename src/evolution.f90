!
! MODULE with routines to perform GA optimizations
!
! Copyright (C) 2005  Alberto Ramos <alberto@martin.ft.uam.es>
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

! $ v. 1.0; Released: 13/05/2008; $

! ***************************************************
! *
MODULE Evolution
! *
! ***************************************************

  USE NumTypes
  USE Error
  USE Constants
  USE Genetics

  ! The external function (i.e. problem dependent) Feval
  ! gives the fitness of an organism
  Interface Feval
     Subroutine Feval(Or)
       USE Genetics
       
       Type (Organism), Intent (in) :: Or
     End Subroutine Feval
  End Interface


!CONTAINS

End MODULE Evolution

