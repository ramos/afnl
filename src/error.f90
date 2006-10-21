!
! A MODULE for Error I/O
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

! $ v. 1.0; Released: 16/09/2004; $

! ***************************************************
! *
MODULE Error
! *
! ***************************************************

  Integer, Parameter :: stderr = 0

CONTAINS

! ***************************************************
! *
  Subroutine Abort(Routine, Msg)
! *
! ***************************************************
  
    Character (len=*), Intent (in), Optional :: Routine
    Character (len=*), Intent (in) :: Msg
    
    If (Present(Routine)) Then
       Write(stderr, *)'  Abort: IN ', Trim(Routine),': ',Msg
    Else
       Write(stderr, *)'  Abort: ', Msg
    End If

    Stop
    
    Return
  End Subroutine Abort

! ***************************************************
! *
  Subroutine Perror(Routine, Msg)
! *
! ***************************************************
  
    Character (len=*), Intent (in), Optional :: Routine
    Character (len=*), Intent (in) :: Msg
    
    If (Present(Routine)) Then
       Write(stderr, *)'  Error: IN ', Trim(Routine),': ',Msg
    Else
       Write(stderr, *)'  Error: ', Msg
    End If

    Return
  End Subroutine Perror


End MODULE Error
