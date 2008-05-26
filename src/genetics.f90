!
! MODULE with routines to perform GA optimizations
!
! Copyright (C) 2008  Alberto Ramos <alberto@martin.ft.uam.es>
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
MODULE Genetics
! *
! ***************************************************

  USE NumTypes
  USE Error
  USE Constants
  USE Statistics

  Type Organism
     Character (len=1), Allocatable :: Genotype(:)
     Integer :: Ngene
     Real (kind=DP) :: Fitness
  End Type Organism

  Type Generation
     Type (Organism), Allocatable :: Memeber(:)
  End Type Generation


  Interface Assignment (=)
     Module Procedure Igual
  End Interface


  Interface Mutate
     Module Procedure Mutate_G, Mutate_Rnd
  End Interface

  Interface Crossover
     Module Procedure Crossover_1pt, Crossover_2pt
  End Interface

  Interface Encode     
     Module Procedure Encode_DP, Encode_SP
  End Interface

  Interface Decode     
     Module Procedure Decode_DP, Decode_SP
  End Interface


  Private Mutate_G, Mutate_Rnd, Crossover_1pt, Crossover_2pt, Igual, &
       & Encode_DP, Encode_SP, Decode_DP, Decode_SP


  ! List of non-mutable Genes. NOTE: It must be sorted!!
  Integer, Allocatable :: Non_Mutable(:)

CONTAINS

! ********************************
! *
  Subroutine Igual(Or1, Or2)
! *
! ********************************
    
    Type (Organism), Intent (in) :: Or2
    Type (Organism), Intent (out) :: Or1

    CALL Init_Organism(Or1, Or2%Ngene)

    Or1%Genotype = Or2%Genotype
    Or1%Fitness  = Or2%Fitness

    Return
  End Subroutine Igual


! ********************************
! *
  Subroutine Init_Organism(Or, Ng, Irnd)
! *
! ********************************

    Type (Organism), Intent (out) :: Or
    Integer, Intent (in) :: Ng
    Logical, Intent (in), Optional :: Irnd

    Integer :: Ic(Ng), I
    Character (len=20) :: fmt

    Allocate(Or%Genotype(Ng))
    Or%Ngene = Ng
    Or%Fitness = -Huge(1.0_DP)

    If (Present(Irnd)) Then
       If (Irnd) Then
          CALL Irand(Ic, 0, 9)
          Write(fmt,*)'(', Ng, 'I1.1)'
          Write(Or%Genotype,fmt)(Ic(I), I=1, Ng)
       End If
    End If


    Return
  End Subroutine Init_Organism

! ********************************
! *
  Subroutine Mutate_G(Or, Ng)
! *
! ********************************

    Type (Organism), Intent (inout) :: Or
    Integer, Intent (in) :: Ng

    Character (len=2) :: ch
    Integer :: Ipos
    
    If (Allocated(Non_Mutable)) Then
       Do Ipos = 1, Size(Non_Mutable)
          If (Ng == Non_Mutable(Ipos)) Return
       End Do
    End If


    ch(1:1) = Or%Genotype(Ng)
    If (ch(1:1) == '-') Then
       Or%Genotype(Ng:Ng) = ' '
    Else If (ch(1:1) == ' ') Then
       Or%Genotype(Ng:Ng) = '-'
    Else If (ch(1:1) /= '.') Then
       Write(*,*)'         *A: ', Or%Genotype(Ng)
       Ipos = Irand(0,9)
       Write(*,*)Ipos
       Write(Or%Genotype(Ng),'(1I1.1)')Ipos
       Write(*,*)'         *D: ', Or%Genotype(Ng)
    End If

    Return
  End Subroutine Mutate_G
  
! ********************************
! *
  Subroutine Mutate_Rnd(Or)
! *
! ********************************

    Type (Organism), Intent (inout) :: Or

    Integer :: Ng
    
    Ng = Irand(1,Or%Ngene)
    CALL Mutate(Or, Ng)

    Return
  End Subroutine Mutate_Rnd

! ********************************
! *
  Subroutine Crossover_1pt(Or1, Or2, Np)
! *
! ********************************
  
    Type (Organism), Intent (inout) :: Or1, Or2
    Integer, Intent (in), Optional :: Np
    
    Character, Allocatable :: Swp(:)
    Integer :: Ncp

    If (Present(Np)) Then
       Ncp = Np
    Else
       Ncp = Irand(1, Min(Or1%Ngene,Or2%Ngene))
    End If
    
    Allocate(Swp(Ncp))

    Swp(1:Ncp) = Or1%Genotype(1:Ncp)
    Or1%Genotype(1:Ncp) = Or2%Genotype(1:Ncp)
    Or2%Genotype(1:Ncp) = Swp(1:Ncp)

    Deallocate(Swp)

    Or1%Fitness = -Huge(1.0_DP)
    Or2%Fitness = -Huge(1.0_DP)

    Return
  End Subroutine Crossover_1pt


! ********************************
! *
  Subroutine Crossover_2pt(Or1, Or2, Np1, Np2)
! *
! ********************************
  
    Type (Organism), Intent (inout) :: Or1, Or2
    Integer, Intent (in) :: Np1, Np2

    Integer :: Ncp1, Ncp2, Ntmp1!, Ntmp2
    Character, Allocatable :: Swp(:)

    Ncp1 = Min(Np1, Np2)
    Ncp2 = Max(Np1, Np2)
!!$    Else
!!$       Ntmp1 = Irand(1, Min(Or1%Ngene,Or2%Ngene))
!!$       Ntmp2 = Irand(1, Min(Or1%Ngene,Or2%Ngene))
!!$       Ncp1 = Min(Ntmp1, Ntmp2)
!!$       Ncp2 = Max(Ntmp1, Ntmp2)
!!$    End If

    Ntmp1 = Ncp2-Ncp1+1
    Allocate(Swp(Ntmp1))

    Swp(:) = Or1%Genotype(Ncp1:Ncp2)
    Or1%Genotype(Ncp1:Ncp2) = Or2%Genotype(Ncp1:Ncp2)
    Or2%Genotype(Ncp1:Ncp2) = Swp(:)

    Deallocate(Swp)    

    Or1%Fitness = -Huge(1.0_DP)
    Or2%Fitness = -Huge(1.0_DP)

    Return
  End Subroutine Crossover_2pt

! ********************************
! *
  Function Encode_DP(X)
! *
! ********************************

    Real (kind=DP), Intent (in) :: X
    
    Character :: Encode_DP(24)
    Character (len=24) :: Str
    Integer :: I

    Write(Str,'(1I1,1F18.15,1I5.4)')Radix(X),Fraction(X),&
         & Exponent(X)

    Do I = 1, 24
       Encode_DP(I) = Str(I:I)
    End Do

    Return
  End Function Encode_DP

! ********************************
! *
  Function Encode_SP(X)
! *
! ********************************

    Real (kind=SP), Intent (in) :: X
    
    Character :: Encode_SP(16)
    Character (len=16) :: Str
    Integer :: I

    Write(Str,'(1I1,1F10.7,1I5.4)')Radix(X),Fraction(X),&
         & Exponent(X)

    Do I = 1, 16
       Encode_SP(I) = Str(I:I)
    End Do

    Return
  End Function Encode_SP

! ********************************
! *
  Subroutine Decode_DP(Str, Decode)
! *
! ********************************

    Character, Intent (in) :: Str(24)
    Real (kind=DP), Intent (out) :: Decode
    
    Character (len=24) :: tmp
    Integer :: rd, ex, I
    Real (kind=DP) :: fr

    Do I = 1, 24
       Tmp(I:I) = Str(I)
    End Do

    Read(Tmp,'(1I1,1F18.15,1I5.4)')rd, fr, ex
    If (ex > MaxExponent(1.0_DP)) ex = MaxExponent(1.0_DP)-1
    If (ex < MinExponent(1.0_DP)) ex = MinExponent(1.0_DP)+1
    Decode = Real(rd,kind=DP)**ex * fr

    Return
  End Subroutine Decode_DP

! ********************************
! *
  Subroutine Decode_SP(Str, Decode)
! *
! ********************************

    Character, Intent (in) :: Str(16)
    Real (kind=SP), Intent (out) :: Decode
    
    Character (len=16) :: Tmp
    Integer :: rd, ex, I
    Real (kind=SP) :: fr

    Do I = 1, 16
       Tmp(I:I) = Str(I)
    End Do

    Read(Tmp,'(1I1,1F10.7,1I5.4)')rd, fr, ex
    If (ex > MaxExponent(1.0_SP)) ex = MaxExponent(1.0_SP)-1
    If (ex < MinExponent(1.0_SP)) ex = MinExponent(1.0_SP)+1
    Decode = Real(rd,kind=DP)**ex * fr
    
    

    Return
  End Subroutine Decode_SP
  

End MODULE Genetics

