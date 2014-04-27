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
MODULE MixMax
! *
! ***************************************************

  USE ISO_FORTRAN_ENV, Only : ERROR_UNIT


  IMPLICIT NONE

  Type, Private :: State
!  Type :: State
     Integer (kind=4) :: N
     Integer (kind=8), Allocatable :: V(:)
     Integer (kind=4) :: cnt
     Integer (kind=8) :: sumtot = 0
  End type State

  Integer (kind=8), Parameter, Private :: &
       & MULT=1073217536_8, BITS = 61_8, MERSBASE = 2305843009213693951_8, &
       & A = 9_8
  Integer (kind=8), Private :: SPECIAL

  Type (State), Private        :: rnd
!  Type (State)       :: rnd
  Procedure (Mulspec), Pointer, Private :: Mod_Mulspec => noinit

  Real (kind=8) :: DINV_MERSBASE=0.433680868994201773791060216479542685926876E-18_8
  
  Interface 
     Function Mulspec(nin)
       Integer (kind=8), Intent (in) :: nin
       Integer (kind=8) :: Mulspec
     End Function Mulspec
  End Interface

  Interface mxmx
     Module Procedure mxmx_int, mxmx_f64
  End Interface mxmx

  Private :: noinit, modmulM61, mxmx_int, mxmx_f64, fill_rnd, Mod_Mersenne, &
       & Mod_Mulspec_direct, Mod_Mulspec_default, Mod_Mulspec_zero, &
       & Mod_Mulspec_one, mxmx_error

CONTAINS

! ***************************************************
! *
  Function modmulM61(ia,ib)
! *
! ***************************************************

    Integer (kind=8 ), Intent (in) :: ia, ib
    Integer (kind=8) :: modmulM61
    Integer (kind=16) :: tmp
    
    tmp = Int(ia,kind=16)*Int(ib,kind=16)
    modmulM61 = Int(Modulo(tmp,Int(MERSBASE,kind=16)),kind=8)

    Return
  End Function modmulM61

! ***************************************************
! *
  Subroutine mxmx_seed_lcg(nseed)
! *
! ***************************************************
    
    Integer (kind=8), Intent (in) :: nseed
    Integer :: I


    rnd%V(1) = iand(nseed,MERSBASE)
    Do I = 2, rnd%N
       rnd%V(I) = modmulM61(rnd%V(I-1), MULT)
       rnd%sumtot = Mod_Mersenne(rnd%V(I)+rnd%sumtot)
    End Do
    rnd%cnt = 1

    Return
  End Subroutine mxmx_seed_lcg

! ***************************************************
! *
  Subroutine mxmx_seed_vielbein(nidx)
! *
! ***************************************************
    
    Integer (kind=8), Intent (in), Optional :: nidx

    If (present(nidx)) Then
       If (nidx > rnd%N) CALL mxmx_error("mxmx_seed_vielbein", &
            & 'index out of range')
       rnd%V(1:nidx-1) = 0_8
       rnd%V(nidx+1:)  = 0_8
       rnd%V(nidx)     = 1_8
    Else
       rnd%V(2:) = 0_8
       rnd%V(1)  = 1_8
    End If
    rnd%cnt = rnd%N

    Return
  End Subroutine mxmx_seed_vielbein

! ***************************************************
! *
  Subroutine mxmx_int(n)
! *
! ***************************************************

    Integer (kind=8), Intent (out) :: n
    Integer (kind=4) :: ipos
    
    ipos = rnd%cnt
    If (ipos > rnd%N) Then
       ipos = 1
       CALL fill_rnd()
    End If
    rnd%cnt = rnd%cnt + 1
    n = rnd%V(ipos)

    Return
  End Subroutine mxmx_int

! ***************************************************
! *
  Subroutine mxmx_f64(d)
! *
! ***************************************************

    Real (kind=8), Intent (out) :: d
    Integer (kind=4) :: ipos
    
    ipos = rnd%cnt
    If (ipos > rnd%N) Then
       ipos = 1
       CALL fill_rnd()
    End If
    rnd%cnt = rnd%cnt + 1
    d = Real(rnd%V(ipos),kind=8)*DINV_MERSBASE

    Return
  End Subroutine mxmx_f64

! ***************************************************
! *
  Subroutine mxmx_init(nin)
! *
! ***************************************************

    Integer, Intent (in), Optional :: nin
    Integer :: nmat

    If (Present(nin)) Then
       nmat = nin
    Else
       nmat = 3150
    End If
    
    Allocate(rnd%V(nmat))
    rnd%N = nmat
    Select Case (nmat)
    Case (1260)
       SPECIAL = 15_8
    Case (3150)
       SPECIAL = -11_8
    Case (1000)
       SPECIAL = 0_8
    Case (720)
       SPECIAL = 1_8
    Case (508)
       SPECIAL = 5_8
    Case (256)
       SPECIAL = -1_8
    Case (88)
       SPECIAL = 1_8
    Case (64)
       SPECIAL = 6_8
    Case (44)
       SPECIAL = 0_8
    Case (40)
       SPECIAL = 1_8
    Case (30)
       SPECIAL = 3_8
    Case (16)
       SPECIAL = 6_8
    Case (10)
       SPECIAL = -1_8
    Case Default
       CALL mxmx_error('mxmx_init', 'Possible values for N are: '//&
            &'3150 1260 1000 720 508 256 88 64 44 40 30 16 10')
    End Select

    Select Case (SPECIAL)
    Case (2:3)
       Mod_Mulspec => Mod_MulSpec_direct
    Case (-3:-2)
       Mod_Mulspec => Mod_MulSpec_direct
    Case (1)
       Mod_Mulspec => Mod_MulSpec_one
    Case (-1)
       Mod_Mulspec => Mod_MulSpec_default
    Case (0)
       Mod_Mulspec => Mod_MulSpec_zero
    Case Default
       Mod_Mulspec => Mod_MulSpec_default
    End Select
    rnd%cnt=1

    Return
  End Subroutine mxmx_init

! ***************************************************
! *
  Subroutine fill_rnd()
! *
! ***************************************************

    Integer (kind=8) :: tmpV, tmpP, oldsumtot, I, tmp2
    Integer (kind=16) :: newsum
    
    rnd%cnt = 1_8
    If (SPECIAL /= 0_8) tmp2 = rnd%V(2)

    oldsumtot  = rnd%sumtot
    newsum = 0_16 
    tmpV = Mod_Mersenne(rnd%V(1) + oldsumtot)
    tmpP = 0
    rnd%V(1) = tmpV
    Do I = 2, rnd%N
       tmpP = Mod_Mersenne(rnd%V(I) + tmpP)
       tmpV = Mod_Mersenne(tmpV + tmpP)
       rnd%V(I) = tmpV
       newsum = newsum + Int(tmpV,kind=16)
    End Do

    If (SPECIAL /= 0_8) Then
       tmp2 = Mod_MulSpec(tmp2)
       rnd%V(3) = Mod_Mersenne(rnd%V(3) + tmp2)
       newsum = newsum + Int(tmp2,kind=16)
    End If
    rnd%sumtot = Int(Modulo(newsum,int(MERSBASE,kind=16)),kind=8)

    Return
  End Subroutine fill_rnd

! ***************************************************
! *
  Function Mod_Mersenne(nin)
! *
! ***************************************************

    Integer (kind=8), Intent (in) :: nin
    Integer (kind=8) :: Mod_Mersenne
    
    Mod_Mersenne = iand(nin,MERSBASE) + ishft(nin, -BITS)

    Return
  End Function Mod_Mersenne

! ***************************************************
! *
  Function Mod_Mulspec_direct(nin) Result (Mod_Mulspec)
! *
! ***************************************************

    Integer (kind=8), Intent (in) :: nin
    Integer (kind=8) :: Mod_Mulspec
    
    Mod_Mulspec = Mod_Mersenne(SPECIAL*nin)

    Return
  End Function Mod_Mulspec_direct

! ***************************************************
! *
  Function Mod_Mulspec_default(nin) Result (Mod_Mulspec)
! *
! ***************************************************

    Integer (kind=8), Intent (in) :: nin
    Integer (kind=8) :: Mod_Mulspec
    
    Mod_Mulspec = modmulM61(SPECIAL,nin)

    Return
  End Function Mod_Mulspec_default

! ***************************************************
! *
  Function Mod_Mulspec_zero(nin) Result (Mod_Mulspec)
! *
! ***************************************************

    Integer (kind=8), Intent (in) :: nin
    Integer (kind=8) :: Mod_Mulspec
    
    Mod_Mulspec = nin
    Mod_Mulspec = 0_8

    Return
  End Function Mod_Mulspec_zero

! ***************************************************
! *
  Function Mod_Mulspec_one(nin) Result (Mod_Mulspec)
! *
! ***************************************************

    Integer (kind=8), Intent (in) :: nin
    Integer (kind=8) :: Mod_Mulspec
    
    Mod_Mulspec = nin

    Return
  End Function Mod_Mulspec_one

! ********************************
! *
    Subroutine mxmx_error(routine, msg)
! *
! ********************************

      Character (len=*), Intent (in) :: routine, msg

      Write(error_unit,*)'In '//Trim(routine)// ' :'
      Write(error_unit,'(5X,1A)')Trim(msg)
      Write(error_unit,*)

      Stop
      
      Return
    End Subroutine Mxmx_error

! ***************************************************
! *
  Function noinit(nin) Result (Mod_Mulspec)
! *
! ***************************************************

    Integer (kind=8), Intent (in) :: nin
    Integer (kind=8) :: Mod_Mulspec
    
    Mod_Mulspec = nin
    CALL mxmx_error('rnd state', 'Need to initialize MIXMAX with mxmx_init')

    Return
  End Function Noinit

End MODULE MixMax
