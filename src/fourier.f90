
! MODULE with definition and operations of Fourier series
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

! **************************************************
! *
MODULE Fourier
! *
! ***************************************************

  USE Constants, ONLY: PI => PI_DP, DPI => TWOPI_DP, UnitImag_DPC
  USE NumTypes

  IMPLICIT NONE

  Type Fourier_Serie
     Complex (kind=DPC), Allocatable :: Coef(:)
     Integer :: Nterm
  End Type Fourier_Serie
  
  Type Fourier_Serie_2D
     Complex (kind=DPC), Allocatable :: Coef(:,:)
     Integer :: Nterm
  End Type Fourier_Serie_2D
  

  Interface Operator (+)
     Module Procedure Add, Add_2D, Pos, Pos_2D
  End Interface
  
  Interface Operator (-)
     Module Procedure Sub, Sub_2D, Neg, Neg_2D
  End Interface
  
  Interface Operator (*)
     Module Procedure Prod, Prod_2D, Prodcte, Prodcte_2D, Prodcte2,&
          & Prodcte2_2D, ProdCcte, ProdCcte_2D, ProdCcte2, &
          & ProdCcte2_2D 
  End Interface
  
  Interface Operator (**)
     Module Procedure ExpS, ExpS_2D
  End Interface
  
  Interface Assignment (=)
     Module Procedure Equal, Equal_2D
  End Interface
  

  Interface Init_Serie
     Module Procedure Init_Serie_1D, Init_Serie_2D
  End Interface

  Interface Equal_Func
     Module Procedure Equal_Func_1D, Equal_Func_2D
  End Interface

  Interface Eval_Serie
     Module Procedure Eval_Serie_1D, Eval_Serie_2D
  End Interface

  Interface Unit
     Module Procedure Unit_1D, Unit_2D
  End Interface

  Interface DFT
     Module Procedure DFT_1D, DFT_2D
  End Interface

  Interface FastFT
     Module Procedure FastFT_1D
  End Interface

  Interface ExpS
     Module Procedure ExpS, ExpS_2D
  End Interface

  Interface NewExp
     Module Procedure NewExp_1D, NewExp_2D
  End Interface

  Interface Save_Serie
     Module Procedure Save_1D, Save_2D
  End Interface

  Interface Read_Serie
     Module Procedure Read_1D, Read_2D
  End Interface

  Interface Conjg
     Module Procedure ConjgFS_1D, ConjgFS_2D
  End Interface

  Interface cfttree
     Module Procedure CfttreeDummy
  End Interface cfttree

  Interface FFT
     Module Procedure FFT_1D, FastFT_1D
  End Interface FFT

  Private PI, DPI, Add, Add_2D, Prod, Prod_2D, Equal, Equal_2D,&
       & Eval_Serie_1D, Eval_Serie_2D, Sub, Sub_2D, Prodcte, &
       & Prodcte_2D, ExpS, ExpS_2D, Neg, Neg_2D, Pos, Pos_2D, &
       & DFT_1D, DFT_2D, NewExp_1D, NewExp_2D, ConjgFS_1D, &
       & ConjgFS_2D

CONTAINS


! ********************************
! *
  Subroutine Unit_1D(Serie, Nterm)
! *
! ********************************
! * Allocate space for the serie
! * Serie, and returns it with
! * the unit serie (0 fourier 
! * coefficient = 1).
! ********************************

    Type (Fourier_Serie), Intent (out) :: Serie
    Integer, Intent (in) :: Nterm

    CALL Init_Serie(Serie, Nterm)

    Serie%Coef = (0.0_DP, 0.0_DP)
    Serie%Coef(0) = (1.0_DP, 0.0_DP)

    Return
  End Subroutine Unit_1D

! ********************************
! *
  Subroutine Unit_2D(Serie, Nterm)
! *
! ********************************
! * Allocate space for the serie
! * Serie, and returns it with
! * the unit serie (00 fourier 
! * coefficient = 1).
! ********************************

    Type (Fourier_Serie_2D), Intent (out) :: Serie
    Integer, Intent (in) :: Nterm

    CALL Init_Serie(Serie, Nterm)

    Serie%Coef = (0.0_DP, 0.0_DP)
    Serie%Coef(0,0) = (1.0_DP, 0.0_DP)

    Return
  End Subroutine Unit_2D

! ********************************
! *
  Type (Fourier_Serie_2D) Function ProdCcte_2D(D, Serie2) Result (Prod)
! *
! ********************************
! * Calculates the product of a 
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie2
    Complex (kind=DPC), Intent (in) :: D

    CALL Init_Serie(Prod, Serie2%Nterm)

    Prod%Coef = D * Serie2%Coef

    Return
  End Function ProdCcte_2D

! ********************************
! *
  Type (Fourier_Serie_2D) Function Prodcte_2D(D, Serie2) Result (Prod)
! *
! ********************************
! * Calculates the product of a 
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie2
    Real (kind=DP), Intent (in) :: D

    CALL Init_Serie(Prod, Serie2%Nterm)

    Prod%Coef = D * Serie2%Coef

    Return
  End Function Prodcte_2D

! ********************************
! *
  Type (Fourier_Serie_2D) Function ProdCcte2_2D(Serie2, D) Result (Prod)
! *
! ********************************
! * Calculates the product of a 
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie2
    Complex (kind=DPC), Intent (in) :: D

    CALL Init_Serie(Prod, Serie2%Nterm)

    Prod%Coef = D * Serie2%Coef

    Return
  End Function ProdCcte2_2D

! ********************************
! *
  Type (Fourier_Serie_2D) Function Prodcte2_2D(Serie2, D) Result (Prod)
! *
! ********************************
! * Calculates the product of a 
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie2
    Real (kind=DP), Intent (in) :: D

    CALL Init_Serie(Prod, Serie2%Nterm)

    Prod%Coef = D * Serie2%Coef

    Return
  End Function Prodcte2_2D

! ********************************
! *
  Type (Fourier_Serie_2D) Function Prod_2D(Serie1, Serie2) Result (Prod)
! *
! ********************************
! * Calculates the product of two
! * Fourier series.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie1, Serie2

    Integer :: I1, I2, N1, N2, Ntot

    If (Serie1%Nterm > Serie2%Nterm) Then
       Ntot = Serie2%Nterm
    Else
       Ntot = Serie1%Nterm
    End If
    

    CALL Init_Serie(Prod, Ntot)
    
    Do N2 = -Ntot, Ntot
       Do N1 = -Ntot, Ntot
          Prod%Coef(N1,N2) = (0.0_DP, 0.0_DP)
          Do I1 = Max(N1 - Ntot, -Ntot), Min(N1 + Ntot, Ntot)
             Do I2 = Max(N2 - Ntot, -Ntot), Min(N2 + Ntot, Ntot)
                Prod%Coef(N1, N2) = Prod%Coef(N1,N2) + &
                     & Serie1%Coef(I1,I2) * Serie2%Coef(N1-I1, N2-I2)
             End Do
          End Do
       End Do
    End Do
    Prod%Nterm = Ntot
    

    Return
  End Function Prod_2D

! ********************************
! *
  Type (Fourier_Serie) Function ProdCcte(D, Serie2)
! *
! ********************************
! * Calculates the product of a
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie2
    Complex (kind=DPC), Intent (in) :: D
    
    CALL Init_Serie(ProdCcte, Serie2%Nterm)

    ProdCcte%Coef = D * Serie2%Coef

    Return
  End Function ProdCcte

! ********************************
! *
  Type (Fourier_Serie) Function Prodcte(D, Serie2)
! *
! ********************************
! * Calculates the product of a
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie2
    Real (kind=DP), Intent (in) :: D
    
    CALL Init_Serie(Prodcte, Serie2%Nterm)

    Prodcte%Coef = D * Serie2%Coef

    Return
  End Function Prodcte

! ********************************
! *
  Type (Fourier_Serie) Function ProdCcte2(Serie2, D)
! *
! ********************************
! * Calculates the product of a
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie2
    Complex (kind=DPC), Intent (in) :: D
    
    CALL Init_Serie(ProdCcte2, Serie2%Nterm)

    ProdCcte2%Coef = D * Serie2%Coef

    Return
  End Function ProdCcte2

! ********************************
! *
  Type (Fourier_Serie) Function Prodcte2(Serie2, D)
! *
! ********************************
! * Calculates the product of a
! * Fourier series and a number.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie2
    Real (kind=DP), Intent (in) :: D
    
    CALL Init_Serie(Prodcte2, Serie2%Nterm)

    Prodcte2%Coef = D * Serie2%Coef

    Return
  End Function Prodcte2

! ********************************
! *
  Type (Fourier_Serie) Function Prod(Serie1, Serie2)
! *
! ********************************
! * Calculates the product of two
! * Fourier series.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie1, Serie2

    Integer :: I, J, Ntot

    If (Serie1%Nterm > Serie2%Nterm) Then
       Ntot = Serie2%Nterm
    Else
       Ntot = Serie1%Nterm
    End If
    

    CALL Init_Serie(Prod, Ntot)
    
    Do J = -Ntot, Ntot
       Prod%Coef(J) = (0.0_DP, 0.0_DP)
       Do I = Max(J - Ntot, -Ntot), Min(J + Ntot, Ntot)
          Prod%Coef(J) = Prod%Coef(J) + &
               & Serie1%Coef(I) * Serie2%Coef(J-I)
       End Do
    End Do
    Prod%Nterm = Ntot

    Return
  End Function Prod
  
! ********************************
! *
  Type (Fourier_Serie) Function ExpS(Serie, Nex)
! *
! ********************************
! * "Exponentiates" a fourier series
! * in a fast way.
! ********************************
    

    Type (Fourier_Serie), Intent (in) :: Serie
    Integer, Intent (in) :: Nex

    Integer, Parameter :: Maxcfr = 30
    Integer :: Bin(Maxcfr), Nmax, I, Ierr
    Type (Fourier_Serie), Allocatable :: Acum(:)

    If (Nex == 0) Then
       CALL Unit(ExpS, Serie%Nterm)
       Return
    End If
    
    ! Firs we need the binary decomposition of 
    ! Nex
    Do I = 1, Maxcfr
       Bin(I) = Int(Mod(Nex, 2**I)/(2**(I-1)))
       If (Bin(I) == 1) Nmax = I
       Write(*,*)Bin(I)
    End Do

    Allocate(Acum(Nmax), STAT = Ierr)
    If (Ierr > 0) Then
       Write(0,*)'ExpS: Out of memory.'
       Stop
    End If
    
    Acum(1) = Serie
    Do I = 2, Nmax
       Acum(I) = Acum(I-1) * Acum(I-1)
    End Do
    
    ExpS = Acum(Nmax)
    Do I = Nmax - 1, 1, -1
       If (Bin(I) == 1) Then
          ExpS = ExpS * Acum(I)
       End If
    End Do
    
    Deallocate(Acum)
   
    Return
  End Function ExpS

! ********************************
! *
  Type (Fourier_Serie) Function NewExp_1D(Serie, Nex) Result (NewExp)
! *
! ********************************
! * "Exponentiates" a fourier series
! * in a fast way.
! ********************************
    

    Type (Fourier_Serie), Intent (in) :: Serie
    Integer, Intent (in) :: Nex

    Integer, Parameter :: Maxcfr = 30
    Integer :: Bin(Maxcfr), Nmax, I
    Type (Fourier_Serie) :: Acum

    If (Nex == 0) Then
       CALL Unit(NewExp, Serie%Nterm)
       Return
    End If
    
    ! Firs we need the binary decomposition of 
    ! Nex
    Do I = 1, Maxcfr
       Bin(I) = Int(Mod(Nex, 2**I)/(2**(I-1)))
       If (Bin(I) == 1) Nmax = I
    End Do

    CALL Unit(Acum, Serie%Nterm)
    Acum = Serie
    CALL Unit(NewExp, Serie%Nterm)
    Do I = 1, Nmax - 1
       If (Bin(I) == 1) Then
          NewExp = NewExp * Acum
       End If
       Acum = Acum * Acum
    End Do

    NewExp = NewExp * Acum

    Deallocate(Acum%Coef)
   
    Return
  End Function NewExp_1D

! ********************************
! *
  Type (Fourier_Serie_2D) Function NewExp_2D(Serie, Nex) Result (NewExp)
! *
! ********************************
! * "Exponentiates" a fourier series
! * in a fast way.
! ********************************
    

    Type (Fourier_Serie_2D), Intent (in) :: Serie
    Integer, Intent (in) :: Nex

    Integer, Parameter :: Maxcfr = 30
    Integer :: Bin(Maxcfr), Nmax, I
    Type (Fourier_Serie_2D) :: Acum

    If (Nex == 0) Then
       CALL Unit(NewExp, Serie%Nterm)
       Return
    End If
    
    ! Firs we need the binary decomposition of 
    ! Nex
    Do I = 1, Maxcfr
       Bin(I) = Int(Mod(Nex, 2**I)/(2**(I-1)))
       If (Bin(I) == 1) Nmax = I
    End Do

    CALL Unit(Acum, Serie%Nterm)
    Acum = Serie
    CALL Unit(NewExp, Serie%Nterm)
    Do I = 1, Nmax - 1
       If (Bin(I) == 1) Then
          NewExp = NewExp * Acum
       End If
       Acum = Acum * Acum
    End Do

    NewExp = NewExp * Acum

    Deallocate(Acum%Coef)
   
    Return
  End Function NewExp_2D
  
! ********************************
! *
  Type (Fourier_Serie_2D) Function ExpS_2D(Serie, Nex)
! *
! ********************************
! * "Exponentiates" a fourier series
! * in a fast way
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie
    Integer, Intent (in) :: Nex

    Integer, Parameter :: Maxcfr = 32
    Integer :: Bin(Maxcfr), Nmax, I, Ierr
    Type (Fourier_Serie_2D), Allocatable :: Acum(:)

!    CALL Init_Serie(ExpS, Serie%Nterm)

    If (Nex == 0) Then
       CALL Unit(ExpS_2D, Serie%Nterm)
       Return
    End If

    ! Firs we need the binary decomposition of 
    ! Nex
    Do I = 1, Maxcfr
       Bin(I) = Int(Mod(Nex, 2**I)/(2**(I-1)))
       If (Bin(I) == 1) Nmax = I
!       Write(*,*)Bin(I)
    End Do

    Allocate(Acum(Nmax), STAT = Ierr)
    If (Ierr > 0) Then
       Write(0,*)'ExpS_2D: Out of memory.'
       Stop
    End If
    
    Acum(1) = Serie
    Do I = 2, Nmax
       Acum(I) = Acum(I-1) * Acum(I-1)
    End Do
    
    ExpS_2D = Acum(Nmax)
    Do I = Nmax - 1, 1, -1
       If (Bin(I) == 1) Then
          ExpS_2D = ExpS_2D * Acum(I)
       End If
    End Do
    
    Deallocate(Acum)
   
    Return
  End Function ExpS_2D
  
! ********************************
! *
  Type (Fourier_Serie) Function Add(Serie1, Serie2)
! *
! ********************************
! * Add two fourier series.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie1, Serie2

!    Complex (kind=DPC), TARGET :: Temp(-Serie1%Nterm:Serie1%Nterm)

!    Add%Nterm = Serie1%Nterm
!    Temp = Serie1%Coef + Serie2%Coef
!    Add%Coef => Temp
    
    If (Serie1%Nterm > Serie2%Nterm) Then
       Add = Serie1
       Add%Coef(-Serie2%Nterm:Serie2%Nterm) = &
            & Add%Coef(-Serie2%Nterm:Serie2%Nterm) + &
            & Serie2%Coef(:) 
    Else
       Add = Serie2
       Add%Coef(-Serie1%Nterm:Serie1%Nterm) = &
            & Add%Coef(-Serie1%Nterm:Serie1%Nterm) + &
            & Serie1%Coef(:) 
    End If
    
    
    Return
  End Function Add      

! ********************************
! *
  Type (Fourier_Serie) Function Neg(Serie)
! *
! ********************************
! * Add two fourier series.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie

    Integer :: Ierr


    ALLOCATE(Neg%Coef(-Serie%Nterm:Serie%Nterm), STAT=Ierr)
    
    If (Ierr > 0) Then
       Write(0,*)'Equal: Out of memory.'
       Stop
    End If
    
    Neg%Coef = - Serie%Coef
    Neg%Nterm = Serie%Nterm
    
    Return
  End Function Neg

! ********************************
! *
  Type (Fourier_Serie) Function Pos(Serie)
! *
! ********************************
! * Add two fourier series.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie

    Integer :: Ierr

    
    ALLOCATE(Pos%Coef(-Serie%Nterm:Serie%Nterm), STAT=Ierr)
    
    If (Ierr > 0) Then
       Write(0,*)'Equal: Out of memory.'
       Stop
    End If
    
    Pos%Coef = Serie%Coef
    Pos%Nterm = Serie%Nterm
    
    Return
  End Function Pos

! ********************************
! *
  Type (Fourier_Serie_2D) Function Neg_2D(Serie) Result (Neg)
! *
! ********************************
! * Calclates the opposite of a fourier
! * serie.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie

    Integer :: Ierr

    
    ALLOCATE(Neg%Coef(-Serie%Nterm:Serie%Nterm, &
         & -Serie%Nterm:Serie%Nterm), STAT=Ierr)
    
    If (Ierr > 0) Then
       Write(0,*)'Equal: Out of memory.'
       Stop
    End If
    
    Neg%Coef = - Serie%Coef
    Neg%Nterm = Serie%Nterm
    
    Return
  End Function Neg_2D

! ********************************
! *
  Type (Fourier_Serie_2D) Function Pos_2D(Serie) Result (Pos)
! *
! ********************************
! * Calclates the opposite of a fourier
! * serie.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie

    Integer :: Ierr

    
    ALLOCATE(Pos%Coef(-Serie%Nterm:Serie%Nterm, &
         & -Serie%Nterm:Serie%Nterm), STAT=Ierr)
    
    If (Ierr > 0) Then
       Write(0,*)'Equal: Out of memory.'
       Stop
    End If
    
    Pos%Coef = Serie%Coef
    Pos%Nterm = Serie%Nterm
    
    Return
  End Function Pos_2D

! ********************************
! *
  Type (Fourier_Serie_2D) Function Add_2D(Serie1, Serie2)
! *
! ********************************
! * Add two fourier series.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie1, Serie2

!    Complex (kind=DPC), TARGET :: Temp(-Serie1%Nterm:Serie1%Nterm,-Serie1%Nterm:Serie1%Nterm)

!    Add_2D%Nterm = Serie1%Nterm
!    Temp = Serie1%Coef + Serie2%Coef
!    Add_2D%Coef => Temp

    
    If (Serie1%Nterm > Serie2%Nterm) Then
       Add_2D = Serie1
       Add_2D%Coef(-Serie2%Nterm:Serie2%Nterm,-Serie2%Nterm:Serie2%Nterm) = &
            & Add_2D%Coef(-Serie2%Nterm:Serie2%Nterm,-Serie2%Nterm:Serie2%Nterm)&
            & + Serie2%Coef(:,:) 
    Else
       Add_2D = Serie2
       Add_2D%Coef(-Serie1%Nterm:Serie1%Nterm,-Serie1%Nterm:Serie1%Nterm) = &
            & Add_2D%Coef(-Serie1%Nterm:Serie1%Nterm,-Serie1%Nterm:Serie1%Nterm) + &
            & Serie1%Coef(:,:) 
    End If

    
    
    Return
  End Function Add_2D

! ********************************
! *
  Type (Fourier_Serie) Function Sub(Serie1, Serie2)
! *
! ********************************
! * Substract two fourier series.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie1, Serie2
    
    If (Serie1%Nterm > Serie2%Nterm) Then
       Sub = Serie1
       Sub%Coef(-Serie2%Nterm:Serie2%Nterm) = &
            & Sub%Coef(-Serie2%Nterm:Serie2%Nterm) - &
            & Serie2%Coef(:) 
    Else
       Sub = Serie2
       Sub%Coef(-Serie1%Nterm:Serie1%Nterm) = &
            & Sub%Coef(-Serie1%Nterm:Serie1%Nterm) - &
            & Serie1%Coef(:)
       
       Sub%Coef(:) = -Sub%Coef(:)
    End If
    
    
    Return
  End Function Sub

! ********************************
! *
  Type (Fourier_Serie_2D) Function Sub_2D(Serie1, Serie2)
! *
! ********************************
! * Substract two fourier series.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie1, Serie2
    
    If (Serie1%Nterm > Serie2%Nterm) Then
       Sub_2D = Serie1
       Sub_2D%Coef(-Serie2%Nterm:Serie2%Nterm,-Serie2%Nterm:Serie2%Nterm) = &
            & Sub_2D%Coef(-Serie2%Nterm:Serie2%Nterm,-Serie2&
            & %Nterm:Serie2%Nterm) - Serie2%Coef(:,:) 

       Sub_2D%Nterm = Serie1%Nterm
    Else
       Sub_2D = Serie2
       Sub_2D%Coef(-Serie1%Nterm:Serie1%Nterm,-Serie2%Nterm:Serie2%Nterm) = &
            & Sub_2D%Coef(-Serie1%Nterm:Serie1%Nterm,-Serie2&
            & %Nterm:Serie2%Nterm) - Serie1%Coef(:,:) 
       
       Sub_2D%Coef(:,:) = -Sub_2D%Coef(:,:)
       Sub_2D%Nterm = Serie2%Nterm
    End If

    
    Return
  End Function Sub_2D

! ********************************
! *
  Subroutine Equal_Func_1D(Serie, Func, Nterm)
! *
! ********************************
! * Assign Serie2 to Serie1
! ********************************

    Integer, Intent (in) :: Nterm
    Type (Fourier_Serie), Intent (out) :: Serie

    Integer :: I

    Interface
       Function Func(N1)
         USE NumTypes

         Integer, Intent (in)  :: N1
         Complex (kind=DPC) :: Func
       End Function Func
    End Interface

    
    CALL Init_Serie(Serie, Nterm)
    
    Do I = -Nterm, Nterm
       Serie%Coef(I) = Func(I)
    End Do
    

    Return
  End Subroutine Equal_Func_1D


! ********************************
! *
  Subroutine Equal_Func_2D(Serie, Func, Nterm)
! *
! ********************************
! * Assign Serie2 to Serie1
! ********************************

    Integer, Intent (in) :: Nterm
    Type (Fourier_Serie_2D), Intent (out) :: Serie

    Integer :: I, J

    Interface
       Function Func(N1, N2)
         USE NumTypes

         Integer, Intent (in)  :: N1, N2
         Complex (kind=DPC) :: Func
       End Function Func
    End Interface

    
    CALL Init_Serie(Serie, Nterm)
    
    Do J = -Nterm, Nterm
       Do I = -Nterm, Nterm
          Serie%Coef(J, I) = Func(J, I)
       End Do
    End Do

    Return
  End Subroutine Equal_Func_2D


! ********************************
! *
  Subroutine Equal(Serie1, Serie2)
! *
! ********************************
! * Assign Serie2 to Serie1
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie2
    Type (Fourier_Serie), Intent (out) :: Serie1

    Integer :: Ierr
    

    ALLOCATE(Serie1%Coef(-Serie2%Nterm:Serie2%Nterm), STAT=Ierr)

    If (Ierr > 0) Then
       Write(0,*)'Equal: Out of memory.'
       Stop
    End If
    
    Serie1%Coef = Serie2%Coef
    Serie1%Nterm = Serie2%Nterm
    
    Return
  End Subroutine Equal

! ********************************
! *
  Subroutine Equal_2D(Serie1, Serie2)
! *
! ********************************
! * Assign Serie2 to Serie1
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie2
    Type (Fourier_Serie_2D), Intent (out) :: Serie1

    Integer :: Ierr
    

    ALLOCATE(&
         & Serie1%Coef&
         & (-Serie2%Nterm:Serie2%Nterm,-Serie2%Nterm:Serie2%Nterm),&
         & STAT=Ierr)

    If (Ierr > 0) Then
       Write(0,*)'Equal: Out of memory.'
       Stop
    End If
    
    Serie1%Coef = Serie2%Coef
    Serie1%Nterm = Serie2%Nterm
    
    Return
  End Subroutine Equal_2D


! ********************************
! *
  Subroutine Init_Serie_1D(Serie, Nespacio)
! *
! ********************************
! * Allocate memory space for  
! * the coefficients.
! ********************************
    
    Integer, Intent (in) :: Nespacio
    Type (Fourier_Serie), Intent (out) :: Serie
  
    Integer :: Ierr

    
    ALLOCATE(Serie%Coef(-Nespacio:Nespacio), STAT=Ierr)

    If (Ierr > 0) Then
       Write(0,*)'Init_Serie: Out of memory.'
       Stop
    End If
    
    Serie%Nterm = Nespacio
    Serie%Coef = (0.0_DP, 0.0_DP)

    Return
  End Subroutine Init_Serie_1D

! ********************************
! *
  Subroutine Init_Serie_2D(Serie, Nespacio)
! *
! ********************************
! * Allocate memory space for  
! * the coefficients.
! ********************************
    
    Integer, Intent (in) :: Nespacio
    Type (Fourier_Serie_2D), Intent (out) :: Serie
  
    Integer :: Ierr

    
    ALLOCATE(Serie%Coef(-Nespacio:Nespacio, -Nespacio:Nespacio),&
         & STAT=Ierr)

    If (Ierr > 0) Then
       Write(0,*)'Init_Serie: Out of memory.'
       Stop
    End If
    
    Serie%Nterm = Nespacio
    Serie%Coef = (0.0_DP, 0.0_DP)

    Return
  End Subroutine Init_Serie_2D
  
! ********************************
! *
  Complex (kind=DPC) Function Eval_Serie_1D(Serie, X, Tx) Result (Val)
! *
! ********************************
! * Evaluates the fourier Serie in
! * the point X with period
! * Tx.
! ********************************
    
    Type (Fourier_Serie), Intent (in) :: Serie
    Real (kind=DP), Intent (in) :: X, Tx

    Integer :: I

    Val = (0.0_DP, 0.0_DP)
    Do I = -Serie%Nterm, Serie%Nterm
       Val = Val + Serie%Coef(I) * exp((0.0_DP, 1.0_DP)*DPI/Tx * X * I)
    End Do


    Return
  End Function Eval_Serie_1D

! ********************************
! *
  Complex (kind=DPC) Function Eval_Serie_2D(Serie, X, Y, Tx, Ty) Result (Val)
! *
! ********************************
! * Evaluates the fourier Serie in
! * the point (X,Y) with periods 
! * Tx and Ty.
! ********************************
    
    Type (Fourier_Serie_2D), Intent (in) :: Serie
    Real (kind=DP), Intent (in) :: X, Y, Tx, Ty

    Integer :: I, J

    Val = (0.0_DP, 0.0_DP)
    Do J = -Serie%Nterm, Serie%Nterm
       Do I = -Serie%Nterm, Serie%Nterm
          Val = Val + Serie%Coef(I,J) * &
               & exp((0.0_DP, 1.0_DP)*DPI/Tx * X * I) * &
               & exp((0.0_DP, 1.0_DP)*DPI/Ty * Y * J)
       End Do
    End Do

    Return
  End Function Eval_Serie_2D

! ********************************
! *
  Type (Fourier_Serie) Function DFT_1D(Data, Isign) Result (DFT)
! *
! ********************************
! * Calculates the Discrete Fourier
! * Transformation from the complex
! * data stored in data.
! ********************************

    Complex (kind=DPC), Intent (in) :: Data(:)
    Integer, Intent (in), Optional :: Isign
    
    Complex (Kind=DPC) :: Factor, Prefactor
    Integer :: Side, I, J, Ntot, Nterm

    If (.not. Present(Isign)) Then
       Side = -1
    Else
       Side = Isign
    End If

    Ntot = Size(Data)
    Nterm = Int(Ntot/2)


    PreFactor = exp( Side * DPI * (0.0_DP, 1.0_DP) / Real(Ntot,kind=DP) )
    CALL Init_Serie(DFT, Nterm)
    
    Do I = -Nterm, Nterm
       Factor = PreFactor ** I
       
       DFT%Coef(I) = (0.0_DP, 0.0_DP)
       Do J = Ntot, 1, -1
          DFT%Coef(I) = ( DFT%Coef(I) + &
               & Data(J) ) * Factor
       End Do
    End Do
    
    If (Side == -1) Then
       DFT%Coef = DFT%Coef / Real(Ntot, kind=DP)
    End If

    Return
  End Function DFT_1D

! ********************************
! *
  Type (Fourier_Serie_2D) Function DFT_2D(Data, Isign) Result (DFT)
! *
! ********************************
! * Calculates the Discrete Fourier
! * Transformation from the complex
! * data stored in data.
! ********************************

    Complex (kind=DPC), Intent (in) :: Data(:,:)
    Integer, Intent (in), Optional :: Isign
    
    Complex (Kind=DPC) :: Factor1, Factor2, Acum(Size(Data, 1)), &
         & PreFactor
    Integer :: Side, K1, K2, N1, N2, Ntot, Nterm

    If (Present(Isign)) Then
       Side = Isign
    Else
       Side = -1
    End If

    Ntot = Size(Data, 1)
    Nterm = Int(Ntot/2)
    
    PreFactor = exp( Side * DPI * (0.0_DP, 1.0_DP) / Ntot )

    CALL Init_Serie(DFT, Nterm)

    Do K1 = -Nterm, Nterm
       Factor1 = PreFactor ** K1
       Do K2 = -Nterm, Nterm
          Factor2 = PreFactor ** K2

          Do N1 = 1, Ntot
             Acum(N1) = (0.0_DP, 0.0_DP)
             Do N2 = Ntot, 1, -1
                Acum(N1)  = ( Acum(N1) + Data(N1,N2) ) * Factor2
             End Do
          End Do

          DFT%Coef(K1,K2) = (0.0_DP, 0.0_DP)
          Do N1 = Ntot, 1, -1
             DFT%Coef(K1,K2) = ( DFT%Coef(K1,K2) + Acum(N1)) * Factor1
          End Do
       End Do
    End Do

    If (Side == -1) Then
       DFT%Coef = DFT%Coef / Real(Ntot, kind=DP)**2
    End If

    Return
  End Function DFT_2D
  
! ********************************
! *
  Subroutine Save_1D(Serie, File)
! *
! ********************************
! * Save the fourier modes of 
! * "Serie" in the File "File"
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie
    Character (len=*), Intent (in) :: File 

    Character (len=21) :: Frmt
    Integer :: N1

    Open (Unit = 23, File = File, Action = "WRITE")

    Write(frmt,'(1A1,1I12.12,1A8)')'(',2*Serie%Nterm+1,'ES33.25)'
    
    Write(23,'(I12)')Serie%Nterm
    Write(23, frmt)(Serie%Coef(N1), N1 = -Serie%Nterm, Serie%Nterm)


    Close (23)

    Return
  End Subroutine Save_1D

! ********************************
! *
  Subroutine Save_2D(Serie, File)
! *
! ********************************
! * Save the fourier modes of 
! * "Serie" in the File "File"
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie
    Character (len=*), Intent (in) :: File 

    Character (len=21) :: Frmt
    Integer :: N1, N2

    Open (Unit = 23, File = File, Action = "WRITE")

    Write(frmt,'(1A1,1I12.12,1A8)')'(',2*Serie%Nterm+1,'ES33.25)'
    
    Write(23,'(I12)')Serie%Nterm
    Do N2 = -Serie%Nterm, Serie%Nterm
       Write(23, frmt)(Serie%Coef(N1, N2), N1 = -Serie%Nterm, Serie%Nterm)
    End Do

    Close (23)

    Return
  End Subroutine Save_2D

! ********************************
! *
  Subroutine Read_1D(Serie, File)
! *
! ********************************
! * Reads the fourier modes from 
! * a the file File
! ********************************

    Character (len=*), Intent (in) :: File 
    Type (Fourier_Serie), Intent (out) :: Serie

    Character (len=21) :: Frmt
    Integer :: N1, Nterm

    Open (Unit = 23, File = File, Action = "READ")

    Read(23,'(I12)')Nterm
    CALL Init_Serie(Serie, Nterm)

    Write(frmt,'(1A1,1I12.12,1A8)')'(',2*Serie%Nterm+1,'ES33.25)'

    Read(23, frmt)(Serie%Coef(N1), N1 = -Serie%Nterm, Serie%Nterm)

    Close (23)

    Return
  End Subroutine Read_1D

! ********************************
! *
  Subroutine Read_2D(Serie, File)
! *
! ********************************
! * Reads the fourier modes from 
! * a the file File
! ********************************

    Character (len=*), Intent (in) :: File 
    Type (Fourier_Serie_2D), Intent (out) :: Serie

    Character (len=21) :: Frmt
    Integer :: N1, N2, Nterm
    
    Open (Unit = 23, File = File, Action = "READ")

    Read(23,'(I12)')Nterm
    CALL Init_Serie(Serie, Nterm)

    Write(frmt,'(1A1,1I12.12,1A8)')'(',2*Serie%Nterm+1,'ES33.25)'

    Do N2 = -Serie%Nterm, Serie%Nterm
       Read(23, frmt)(Serie%Coef(N1, N2), N1 = -Serie%Nterm, Serie%Nterm)
    End Do

    Close (23)

    Return
  End Subroutine Read_2D

! ********************************
! *
  Type (Fourier_Serie) Function ConjgFS_1D(Serie) Result (Cserie)
! *
! ********************************
! * Calculates the Conjugate of a Fourier 
! * Series.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Serie

    Integer :: K1

    CALL Init_Serie(Cserie, Serie%Nterm)
    
    Do K1 = -Serie%Nterm, Serie%Nterm
       Cserie%Coef(K1) = Conjg(Serie%Coef(-K1))
    End Do

    Return
  End Function ConjgFS_1D

! ********************************
! *
  Type (Fourier_Serie_2D) Function ConjgFS_2D(Serie) Result (Cserie)
! *
! ********************************
! * Calculates the Conjugate of a Fourier 
! * Series.
! ********************************

    Type (Fourier_Serie_2D), Intent (in) :: Serie

    Integer :: K1, K2

    CALL Init_Serie(Cserie, Serie%Nterm)
    
    Do K1 = -Serie%Nterm, Serie%Nterm
       Do K2 = -Serie%Nterm, Serie%Nterm
          Cserie%Coef(K1,K2) = Conjg(Serie%Coef(-K1,-K2))
       End Do
    End Do

    Return
  End Function ConjgFS_2D


! ********************************
! *
  Function FastFT_1D(Data, Isign) Result (FFT)
! *
! ********************************
! * Calculates the Discrete Fourier
! * Transformation from the complex
! * data stored in data.
! ********************************

    Type (Fourier_Serie), Intent (in) :: Data
    Integer, Intent (in), Optional :: Isign
    Type (Fourier_Serie) :: FFT

    Real (kind=DP), Allocatable :: Cdt(:), W(:)
    Integer :: N, Is, I
    Integer, Allocatable :: Ipt(:), Ip(:)

    N = 2*Data%Nterm
    Is = -1
    If (Present(Isign)) Is = Isign

    ALLOCATE(Cdt(0:2*N-1), Ipt(0:N-1))
    ForAll (I=0:Data%Nterm) Ipt(I) = I
    ForAll (I=1:Data%Nterm-1) Ipt(Data%Nterm+I) = -Data%Nterm+I

    
    ALLOCATE(Ip(0:3+Int(Sqrt(Real(N)))), W(0:Int(N/2-1)))
    Do I = 0, N-1
       Cdt(2*I)   = Real (Data%Coef(Ipt(I)))
       Cdt(2*I+1) = Aimag(Data%Coef(Ipt(I)))
    End Do

    Ip(0) = 0
    CALL CDFT(2*N, Is, Cdt, Ip, w)

    CALL Init_Serie(FFT, Data%Nterm)
    Deallocate(Ipt)
    Allocate(Ipt(-Data%Nterm:Data%Nterm))

    ForAll (I=0:Data%Nterm) Ipt(I) = I
    Ipt(-Data%Nterm) = Data%Nterm
    ForAll (I=1:Data%Nterm-1) Ipt(-Data%Nterm+I) = Data%Nterm+I

    Do I = -Data%Nterm, Data%Nterm
       FFT%Coef(I) = Cmplx( Cdt(2*Ipt(I)), Cdt(2*Ipt(I)+1) )
    End Do

    If (Is == -1) Then
       FFT%Coef = FFT%Coef / Real(N, kind=DP)
    End If

    Return
  End Function FASTFT_1D


! ********************************
! *
  Function FFT_1D(Data, Isign) Result (FFT)
! *
! ********************************
! * Calculates the Discrete Fourier
! * Transformation from the complex
! * data stored in data.
! ********************************

    Complex (kind=DPC), Intent (in) :: Data(:)
    Integer, Intent (in), Optional :: Isign
    Type (Fourier_Serie) :: FFT

    Real (kind=DP), Allocatable :: Cdt(:), W(:)
    Integer :: N, Is, I
    Integer, Allocatable :: Ipt(:), Ip(:)

    N = 2*Size(Data)
    Is = -1
    If (Present(Isign)) Is = Isign

    ALLOCATE(Cdt(0:N-1))
    ALLOCATE(Ip(0:3+Int(Sqrt(Real(N)))), W(0:Int(N/2-1)))
    Do I = 1, Size(Data)-1
       Cdt(2*I)   = Real(Data(I))
       Cdt(2*I+1) = Aimag(Data(I))
    End Do
    Cdt(0) = Real(Data(N/2))
    Cdt(1) = Aimag(Data(N/2))

    Ip(0) = 0
    CALL CDFT(N, Is, Cdt, Ip, w)

    CALL Init_Serie(FFT, N/4)
    Allocate(Ipt(-N/4:N/4))

    ForAll (I=0:N/4) Ipt(I) = I
    Ipt(-N/4) = N/4
    ForAll (I=1:N/4-1) Ipt(-N/4+I) = N/4+I

    Do I = -N/4, N/4
       FFT%Coef(I) = Cmplx( Cdt(2*Ipt(I)), Cdt(2*Ipt(I)+1) )
    End Do

    If (Is == -1) Then
       FFT%Coef = 2.0_DP*FFT%Coef / Real(N, kind=DP)
    End If

    Return
  End Function FFT_1D


! ******************************************
! *
! * FFT ROUTINES TAKEN FROM FFT Package:
! *
! * http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
! *
! * Copyright Takuya OOURA, 1996-2001
! *
! * LICENSE: 
! * You may use, copy, modify and distribute this code for any purpose
! * (include commercial use) and without fee. Please refer to this
! * package when you modify this code. 
! ******************************************
    
  subroutine cdft(n, isgn, a, ip, w)

    integer n, isgn, ip(0 : *), nw
    Real (kind=DP) :: a(0 : n - 1), w(0 : *)
    nw = ip(0)
    if (n .gt. 4 * nw) then
       nw = n / 4
       call makewt(nw, ip, w)
    end if
    if (isgn .ge. 0) then
       call cftfsub(n, a, ip, nw, w)
    else
       call cftbsub(n, a, ip, nw, w)
    end if
  end subroutine cdft
  !


  subroutine rdft(n, isgn, a, ip, w)

    integer n, isgn, ip(0 : *), nw, nc
    Real (kind=DP) :: a(0 : n - 1), w(0 : *), xi
    nw = ip(0)
    if (n .gt. 4 * nw) then
       nw = n / 4
       call makewt(nw, ip, w)
    end if
    nc = ip(1)
    if (n .gt. 4 * nc) then
       nc = n / 4
       call makect(nc, ip, w(nw))
    end if
    if (isgn .ge. 0) then
       if (n .gt. 4) then
          call cftfsub(n, a, ip, nw, w)
          call rftfsub(n, a, nc, w(nw))
       else if (n .eq. 4) then
          call cftfsub(n, a, ip, nw, w)
       end if
       xi = a(0) - a(1)
       a(0) = a(0) + a(1)
       a(1) = xi
    else
       a(1) = 0.5d0 * (a(0) - a(1))
       a(0) = a(0) - a(1)
       if (n .gt. 4) then
          call rftbsub(n, a, nc, w(nw))
          call cftbsub(n, a, ip, nw, w)
       else if (n .eq. 4) then
          call cftbsub(n, a, ip, nw, w)
       end if
    end if
  end subroutine rdft
  !


  subroutine ddct(n, isgn, a, ip, w)

    integer n, isgn, ip(0 : *), j, nw, nc
    Real (kind=DP) :: a(0 : n - 1), w(0 : *), xr
    nw = ip(0)
    if (n .gt. 4 * nw) then
       nw = n / 4
       call makewt(nw, ip, w)
    end if
    nc = ip(1)
    if (n .gt. nc) then
       nc = n
       call makect(nc, ip, w(nw))
    end if
    if (isgn .lt. 0) then
       xr = a(n - 1)
       do j = n - 2, 2, -2
          a(j + 1) = a(j) - a(j - 1)
          a(j) = a(j) + a(j - 1)
       end do
       a(1) = a(0) - xr
       a(0) = a(0) + xr
       if (n .gt. 4) then
          call rftbsub(n, a, nc, w(nw))
          call cftbsub(n, a, ip, nw, w)
       else if (n .eq. 4) then
          call cftbsub(n, a, ip, nw, w)
       end if
    end if
    call dctsub(n, a, nc, w(nw))
    if (isgn .ge. 0) then
       if (n .gt. 4) then
          call cftfsub(n, a, ip, nw, w)
          call rftfsub(n, a, nc, w(nw))
       else if (n .eq. 4) then
          call cftfsub(n, a, ip, nw, w)
       end if
       xr = a(0) - a(1)
       a(0) = a(0) + a(1)
       do j = 2, n - 2, 2
          a(j - 1) = a(j) - a(j + 1)
          a(j) = a(j) + a(j + 1)
       end do
       a(n - 1) = xr
    end if
  end subroutine ddct
  !


  subroutine ddst(n, isgn, a, ip, w)
    integer n, isgn, ip(0 : *), j, nw, nc
    Real (kind=DP) :: a(0 : n - 1), w(0 : *), xr
    nw = ip(0)
    if (n .gt. 4 * nw) then
       nw = n / 4
       call makewt(nw, ip, w)
    end if
    nc = ip(1)
    if (n .gt. nc) then
       nc = n
       call makect(nc, ip, w(nw))
    end if
    if (isgn .lt. 0) then
       xr = a(n - 1)
       do j = n - 2, 2, -2
          a(j + 1) = -a(j) - a(j - 1)
          a(j) = a(j) - a(j - 1)
       end do
       a(1) = a(0) + xr
       a(0) = a(0) - xr
       if (n .gt. 4) then
          call rftbsub(n, a, nc, w(nw))
          call cftbsub(n, a, ip, nw, w)
       else if (n .eq. 4) then
          call cftbsub(n, a, ip, nw, w)
       end if
    end if
    call dstsub(n, a, nc, w(nw))
    if (isgn .ge. 0) then
       if (n .gt. 4) then
          call cftfsub(n, a, ip, nw, w)
          call rftfsub(n, a, nc, w(nw))
       else if (n .eq. 4) then
          call cftfsub(n, a, ip, nw, w)
       end if
       xr = a(0) - a(1)
       a(0) = a(0) + a(1)
       do j = 2, n - 2, 2
          a(j - 1) = -a(j) - a(j + 1)
          a(j) = a(j) - a(j + 1)
       end do
       a(n - 1) = -xr
    end if
  end subroutine ddst
  !


  subroutine dfct(n, a, t, ip, w)

    integer n, ip(0 : *), j, k, l, m, mh, nw, nc
    Real (kind=DP) :: a(0 : n), t(0 : n / 2), w(0 : *), xr, xi, yr, yi
    nw = ip(0)
    if (n .gt. 8 * nw) then
       nw = n / 8
       call makewt(nw, ip, w)
    end if
    nc = ip(1)
    if (n .gt. 2 * nc) then
       nc = n / 2
       call makect(nc, ip, w(nw))
    end if
    m = n / 2
    yi = a(m)
    xi = a(0) + a(n)
    a(0) = a(0) - a(n)
    t(0) = xi - yi
    t(m) = xi + yi
    if (n .gt. 2) then
       mh = m / 2
       do j = 1, mh - 1
          k = m - j
          xr = a(j) - a(n - j)
          xi = a(j) + a(n - j)
          yr = a(k) - a(n - k)
          yi = a(k) + a(n - k)
          a(j) = xr
          a(k) = yr
          t(j) = xi - yi
          t(k) = xi + yi
       end do
       t(mh) = a(mh) + a(n - mh)
       a(mh) = a(mh) - a(n - mh)
       call dctsub(m, a, nc, w(nw))
       if (m .gt. 4) then
          call cftfsub(m, a, ip, nw, w)
          call rftfsub(m, a, nc, w(nw))
       else if (m .eq. 4) then
          call cftfsub(m, a, ip, nw, w)
       end if
       a(n - 1) = a(0) - a(1)
       a(1) = a(0) + a(1)
       do j = m - 2, 2, -2
          a(2 * j + 1) = a(j) + a(j + 1)
          a(2 * j - 1) = a(j) - a(j + 1)
       end do
       l = 2
       m = mh
       do while (m .ge. 2)
          call dctsub(m, t, nc, w(nw))
          if (m .gt. 4) then
             call cftfsub(m, t, ip, nw, w)
             call rftfsub(m, t, nc, w(nw))
          else if (m .eq. 4) then
             call cftfsub(m, t, ip, nw, w)
          end if
          a(n - l) = t(0) - t(1)
          a(l) = t(0) + t(1)
          k = 0
          do j = 2, m - 2, 2
             k = k + 4 * l
             a(k - l) = t(j) - t(j + 1)
             a(k + l) = t(j) + t(j + 1)
          end do
          l = 2 * l
          mh = m / 2
          do j = 0, mh - 1
             k = m - j
             t(j) = t(m + k) - t(m + j)
             t(k) = t(m + k) + t(m + j)
          end do
          t(mh) = t(m + mh)
          m = mh
       end do
       a(l) = t(0)
       a(n) = t(2) - t(1)
       a(0) = t(2) + t(1)
    else
       a(1) = a(0)
       a(2) = t(0)
       a(0) = t(1)
    end if
  end subroutine dfct
  !


  subroutine dfst(n, a, t, ip, w)

    integer n, ip(0 : *), j, k, l, m, mh, nw, nc
    Real (kind=DP) :: a(0 : n - 1), t(0 : n / 2 - 1), w(0 : *), xr, xi, yr, yi
    nw = ip(0)
    if (n .gt. 8 * nw) then
       nw = n / 8
       call makewt(nw, ip, w)
    end if
    nc = ip(1)
    if (n .gt. 2 * nc) then
       nc = n / 2
       call makect(nc, ip, w(nw))
    end if
    if (n .gt. 2) then
       m = n / 2
       mh = m / 2
       do j = 1, mh - 1
          k = m - j
          xr = a(j) + a(n - j)
          xi = a(j) - a(n - j)
          yr = a(k) + a(n - k)
          yi = a(k) - a(n - k)
          a(j) = xr
          a(k) = yr
          t(j) = xi + yi
          t(k) = xi - yi
       end do
       t(0) = a(mh) - a(n - mh)
       a(mh) = a(mh) + a(n - mh)
       a(0) = a(m)
       call dstsub(m, a, nc, w(nw))
       if (m .gt. 4) then
          call cftfsub(m, a, ip, nw, w)
          call rftfsub(m, a, nc, w(nw))
       else if (m .eq. 4) then
          call cftfsub(m, a, ip, nw, w)
       end if
       a(n - 1) = a(1) - a(0)
       a(1) = a(0) + a(1)
       do j = m - 2, 2, -2
          a(2 * j + 1) = a(j) - a(j + 1)
          a(2 * j - 1) = -a(j) - a(j + 1)
       end do
       l = 2
       m = mh
       do while (m .ge. 2)
          call dstsub(m, t, nc, w(nw))
          if (m .gt. 4) then
             call cftfsub(m, t, ip, nw, w)
             call rftfsub(m, t, nc, w(nw))
          else if (m .eq. 4) then
             call cftfsub(m, t, ip, nw, w)
          end if
          a(n - l) = t(1) - t(0)
          a(l) = t(0) + t(1)
          k = 0
          do j = 2, m - 2, 2
             k = k + 4 * l
             a(k - l) = -t(j) - t(j + 1)
             a(k + l) = t(j) - t(j + 1)
          end do
          l = 2 * l
          mh = m / 2
          do j = 1, mh - 1
             k = m - j
             t(j) = t(m + k) + t(m + j)
             t(k) = t(m + k) - t(m + j)
          end do
          t(0) = t(m + mh)
          m = mh
       end do
       a(l) = t(0)
    end if
    a(0) = 0
  end subroutine dfst


  ! -------- initializing routines --------
  !

  subroutine makewt(nw, ip, w)

    integer nw, ip(0 : *), j, nwh, nw0, nw1
    Real (kind=DP) :: w(0 : nw - 1), delta, wn4r, wk1r, wk1i, wk3r, wk3i
    ip(0) = nw
    ip(1) = 1
    if (nw .gt. 2) then
       nwh = nw / 2
       delta = atan(1.0d0) / nwh
       wn4r = cos(delta * nwh)
       w(0) = 1
       w(1) = wn4r
       if (nwh .eq. 4) then
          w(2) = cos(delta * 2)
          w(3) = sin(delta * 2)
       else if (nwh .gt. 4) then
          call makeipt(nw, ip)
          w(2) = 0.5d0 / cos(delta * 2)
          w(3) = 0.5d0 / cos(delta * 6)
          do j = 4, nwh - 4, 4
             w(j) = cos(delta * j)
             w(j + 1) = sin(delta * j)
             w(j + 2) = cos(3 * delta * j)
             w(j + 3) = -sin(3 * delta * j)
          end do
       end if
       nw0 = 0
       do while (nwh .gt. 2)
          nw1 = nw0 + nwh
          nwh = nwh / 2
          w(nw1) = 1
          w(nw1 + 1) = wn4r
          if (nwh .eq. 4) then
             wk1r = w(nw0 + 4)
             wk1i = w(nw0 + 5)
             w(nw1 + 2) = wk1r
             w(nw1 + 3) = wk1i
          else if (nwh .gt. 4) then
             wk1r = w(nw0 + 4)
             wk3r = w(nw0 + 6)
             w(nw1 + 2) = 0.5d0 / wk1r
             w(nw1 + 3) = 0.5d0 / wk3r
             do j = 4, nwh - 4, 4
                wk1r = w(nw0 + 2 * j)
                wk1i = w(nw0 + 2 * j + 1)
                wk3r = w(nw0 + 2 * j + 2)
                wk3i = w(nw0 + 2 * j + 3)
                w(nw1 + j) = wk1r
                w(nw1 + j + 1) = wk1i
                w(nw1 + j + 2) = wk3r
                w(nw1 + j + 3) = wk3i
             end do
          end if
          nw0 = nw1
       end do
    end if
  end subroutine makewt
  !


  subroutine makeipt(nw, ip)

    integer nw, ip(0 : *), j, l, m, m2, p, q
    ip(2) = 0
    ip(3) = 16
    m = 2
    l = nw
    do while (l .gt. 32)
       m2 = 2 * m
       q = 8 * m2
       do j = m, m2 - 1
          p = 4 * ip(j)
          ip(m + j) = p
          ip(m2 + j) = p + q
       end do
       m = m2
       l = l / 4
    end do
  end subroutine makeipt
  !


  subroutine makect(nc, ip, c)

    integer nc, ip(0 : *), j, nch
    Real (kind=DP) :: c(0 : nc - 1), delta
    ip(1) = nc
    if (nc .gt. 1) then
       nch = nc / 2
       delta = atan(1.0d0) / nch
       c(0) = cos(delta * nch)
       c(nch) = 0.5d0 * c(0)
       do j = 1, nch - 1
          c(j) = 0.5d0 * cos(delta * j)
          c(nc - j) = 0.5d0 * sin(delta * j)
       end do
    end if
  end subroutine makect
  !
  ! -------- child routines --------
  !
  function cfttreedummy(n, j, k, a, nw, w) Result (cfttree)
    
    integer n, j, k, nw, i, isplt, m, cfttree
    Real (kind=DP) :: a(0 : j - 1), w(0 : nw - 1)
    if (mod(k, 4) .ne. 0) then
       isplt = mod(k, 2)
       if (isplt .ne. 0) then
          call cftmdl1(n, a(j - n), w(nw - n / 2))
       else
          call cftmdl2(n, a(j - n), w(nw - n))
       end if
    else
       m = n
       i = k
       do while (mod(i, 4) .eq. 0)
          m = m * 4
          i = i / 4
       end do
       isplt = mod(i, 2)
       if (isplt .ne. 0) then
          do while (m .gt. 128)
             call cftmdl1(m, a(j - m), w(nw - m / 2))
             m = m / 4
          end do
       else
          do while (m .gt. 128)
             call cftmdl2(m, a(j - m), w(nw - m))
             m = m / 4
          end do
       end if
    end if
    cfttree = isplt

  end function cfttreedummy
  !
  subroutine cftfsub(n, a, ip, nw, w)
    integer n, ip(0 : *), nw
    Real (kind=DP) :: a(0 : n - 1), w(0 : nw - 1)
    if (n .gt. 8) then
       if (n .gt. 32) then
          call cftf1st(n, a, w(nw - n / 4))
          if (n .gt. 512) then
             call cftrec4(n, a, nw, w)
          else if (n .gt. 128) then
             call cftleaf(n, 1, a, nw, w)
          else
             call cftfx41(n, a, nw, w)
          end if
          call bitrv2(n, ip, a)
       else if (n .eq. 32) then
          call cftf161(a, w(nw - 8))
          call bitrv216(a)
       else
          call cftf081(a, w)
          call bitrv208(a)
       end if
    else if (n .eq. 8) then
       call cftf040(a)
    else if (n .eq. 4) then
       call cftx020(a)
    end if
  end subroutine cftfsub
  !


  subroutine cftbsub(n, a, ip, nw, w)

    integer n, ip(0 : *), nw
    Real (kind=DP) :: a(0 : n - 1), w(0 : nw - 1)
    if (n .gt. 8) then
       if (n .gt. 32) then
          call cftb1st(n, a, w(nw - n / 4))
          if (n .gt. 512) then
             call cftrec4(n, a, nw, w)
          else if (n .gt. 128) then
             call cftleaf(n, 1, a, nw, w)
          else
             call cftfx41(n, a, nw, w)
          end if
          call bitrv2conj(n, ip, a)
       else if (n .eq. 32) then
          call cftf161(a, w(nw - 8))
          call bitrv216neg(a)
       else
          call cftf081(a, w)
          call bitrv208neg(a)
       end if
    else if (n .eq. 8) then
       call cftb040(a)
    else if (n .eq. 4) then
       call cftx020(a)
    end if
  end subroutine cftbsub
  !


  subroutine bitrv2(n, ip, a)

    integer n, ip(0 : *), j, j1, k, k1, l, m, nh, nm
    Real (kind=DP) :: a(0 : n - 1), xr, xi, yr, yi
    m = 1
    l = n / 4
    do while (l .gt. 8)
       m = m * 2
       l = l / 4
    end do
    nh = n / 2
    nm = 4 * m
    if (l .eq. 8) then
       do k = 0, m - 1
          do j = 0, k - 1
             j1 = 4 * j + 2 * ip(m + k)
             k1 = 4 * k + 2 * ip(m + j)
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 - nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nh
             k1 = k1 + 2
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 + nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + 2
             k1 = k1 + nh
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 - nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nh
             k1 = k1 - 2
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 + nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
          end do
          k1 = 4 * k + 2 * ip(m + k)
          j1 = k1 + 2
          k1 = k1 + nh
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 + nm
          k1 = k1 + 2 * nm
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 + nm
          k1 = k1 - nm
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 - 2
          k1 = k1 - nh
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 + nh + 2
          k1 = k1 + nh + 2
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 - nh + nm
          k1 = k1 + 2 * nm - 2
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
       end do
    else
       do k = 0, m - 1
          do j = 0, k - 1
             j1 = 4 * j + ip(m + k)
             k1 = 4 * k + ip(m + j)
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nh
             k1 = k1 + 2
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + 2
             k1 = k1 + nh
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nh
             k1 = k1 - 2
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - nm
             xr = a(j1)
             xi = a(j1 + 1)
             yr = a(k1)
             yi = a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
          end do
          k1 = 4 * k + ip(m + k)
          j1 = k1 + 2
          k1 = k1 + nh
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 + nm
          k1 = k1 + nm
          xr = a(j1)
          xi = a(j1 + 1)
          yr = a(k1)
          yi = a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
       end do
    end if
  end subroutine bitrv2
  !
  subroutine bitrv2conj(n, ip, a)
    integer n, ip(0 : *), j, j1, k, k1, l, m, nh, nm
    Real (kind=DP) :: a(0 : n - 1), xr, xi, yr, yi
    m = 1
    l = n / 4
    do while (l .gt. 8)
       m = m * 2
       l = l / 4
    end do
    nh = n / 2
    nm = 4 * m
    if (l .eq. 8) then
       do k = 0, m - 1
          do j = 0, k - 1
             j1 = 4 * j + 2 * ip(m + k)
             k1 = 4 * k + 2 * ip(m + j)
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 - nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nh
             k1 = k1 + 2
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 + nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + 2
             k1 = k1 + nh
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 - nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nh
             k1 = k1 - 2
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 + nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - 2 * nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
          end do
          k1 = 4 * k + 2 * ip(m + k)
          j1 = k1 + 2
          k1 = k1 + nh
          a(j1 - 1) = -a(j1 - 1)
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          a(k1 + 3) = -a(k1 + 3)
          j1 = j1 + nm
          k1 = k1 + 2 * nm
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 + nm
          k1 = k1 - nm
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 - 2
          k1 = k1 - nh
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 + nh + 2
          k1 = k1 + nh + 2
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          j1 = j1 - nh + nm
          k1 = k1 + 2 * nm - 2
          a(j1 - 1) = -a(j1 - 1)
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          a(k1 + 3) = -a(k1 + 3)
       end do
    else
       do k = 0, m - 1
          do j = 0, k - 1
             j1 = 4 * j + ip(m + k)
             k1 = 4 * k + ip(m + j)
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nh
             k1 = k1 + 2
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + 2
             k1 = k1 + nh
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + nm
             k1 = k1 + nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nh
             k1 = k1 - 2
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 - nm
             k1 = k1 - nm
             xr = a(j1)
             xi = -a(j1 + 1)
             yr = a(k1)
             yi = -a(k1 + 1)
             a(j1) = yr
             a(j1 + 1) = yi
             a(k1) = xr
             a(k1 + 1) = xi
          end do
          k1 = 4 * k + ip(m + k)
          j1 = k1 + 2
          k1 = k1 + nh
          a(j1 - 1) = -a(j1 - 1)
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          a(k1 + 3) = -a(k1 + 3)
          j1 = j1 + nm
          k1 = k1 + nm
          a(j1 - 1) = -a(j1 - 1)
          xr = a(j1)
          xi = -a(j1 + 1)
          yr = a(k1)
          yi = -a(k1 + 1)
          a(j1) = yr
          a(j1 + 1) = yi
          a(k1) = xr
          a(k1 + 1) = xi
          a(k1 + 3) = -a(k1 + 3)
       end do
    end if
  end subroutine bitrv2conj
  !
  subroutine bitrv216(a)
    Real (kind=DP) :: a(0 : 31), x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i
    Real (kind=DP) :: x5r, x5i, x7r, x7i, x8r, x8i, x10r, x10i
    Real (kind=DP) :: x11r, x11i, x12r, x12i, x13r, x13i, x14r, x14i
    x1r = a(2)
    x1i = a(3)
    x2r = a(4)
    x2i = a(5)
    x3r = a(6)
    x3i = a(7)
    x4r = a(8)
    x4i = a(9)
    x5r = a(10)
    x5i = a(11)
    x7r = a(14)
    x7i = a(15)
    x8r = a(16)
    x8i = a(17)
    x10r = a(20)
    x10i = a(21)
    x11r = a(22)
    x11i = a(23)
    x12r = a(24)
    x12i = a(25)
    x13r = a(26)
    x13i = a(27)
    x14r = a(28)
    x14i = a(29)
    a(2) = x8r
    a(3) = x8i
    a(4) = x4r
    a(5) = x4i
    a(6) = x12r
    a(7) = x12i
    a(8) = x2r
    a(9) = x2i
    a(10) = x10r
    a(11) = x10i
    a(14) = x14r
    a(15) = x14i
    a(16) = x1r
    a(17) = x1i
    a(20) = x5r
    a(21) = x5i
    a(22) = x13r
    a(23) = x13i
    a(24) = x3r
    a(25) = x3i
    a(26) = x11r
    a(27) = x11i
    a(28) = x7r
    a(29) = x7i
  end subroutine bitrv216
  !
  subroutine bitrv216neg(a)
    Real (kind=DP) :: a(0 : 31), x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i
    Real (kind=DP) :: x5r, x5i, x6r, x6i, x7r, x7i, x8r, x8i
    Real (kind=DP) :: x9r, x9i, x10r, x10i, x11r, x11i, x12r, x12i
    Real (kind=DP) :: x13r, x13i, x14r, x14i, x15r, x15i
    x1r = a(2)
    x1i = a(3)
    x2r = a(4)
    x2i = a(5)
    x3r = a(6)
    x3i = a(7)
    x4r = a(8)
    x4i = a(9)
    x5r = a(10)
    x5i = a(11)
    x6r = a(12)
    x6i = a(13)
    x7r = a(14)
    x7i = a(15)
    x8r = a(16)
    x8i = a(17)
    x9r = a(18)
    x9i = a(19)
    x10r = a(20)
    x10i = a(21)
    x11r = a(22)
    x11i = a(23)
    x12r = a(24)
    x12i = a(25)
    x13r = a(26)
    x13i = a(27)
    x14r = a(28)
    x14i = a(29)
    x15r = a(30)
    x15i = a(31)
    a(2) = x15r
    a(3) = x15i
    a(4) = x7r
    a(5) = x7i
    a(6) = x11r
    a(7) = x11i
    a(8) = x3r
    a(9) = x3i
    a(10) = x13r
    a(11) = x13i
    a(12) = x5r
    a(13) = x5i
    a(14) = x9r
    a(15) = x9i
    a(16) = x1r
    a(17) = x1i
    a(18) = x14r
    a(19) = x14i
    a(20) = x6r
    a(21) = x6i
    a(22) = x10r
    a(23) = x10i
    a(24) = x2r
    a(25) = x2i
    a(26) = x12r
    a(27) = x12i
    a(28) = x4r
    a(29) = x4i
    a(30) = x8r
    a(31) = x8i
  end subroutine bitrv216neg
  !
  subroutine bitrv208(a)
    Real (kind=DP) :: a(0 : 15), x1r, x1i, x3r, x3i, x4r, x4i, x6r, x6i
    x1r = a(2)
    x1i = a(3)
    x3r = a(6)
    x3i = a(7)
    x4r = a(8)
    x4i = a(9)
    x6r = a(12)
    x6i = a(13)
    a(2) = x4r
    a(3) = x4i
    a(6) = x6r
    a(7) = x6i
    a(8) = x1r
    a(9) = x1i
    a(12) = x3r
    a(13) = x3i
  end subroutine bitrv208
  !
  subroutine bitrv208neg(a)
    Real (kind=DP) :: a(0 : 15), x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i
    Real (kind=DP) :: x5r, x5i, x6r, x6i, x7r, x7i
    x1r = a(2)
    x1i = a(3)
    x2r = a(4)
    x2i = a(5)
    x3r = a(6)
    x3i = a(7)
    x4r = a(8)
    x4i = a(9)
    x5r = a(10)
    x5i = a(11)
    x6r = a(12)
    x6i = a(13)
    x7r = a(14)
    x7i = a(15)
    a(2) = x7r
    a(3) = x7i
    a(4) = x3r
    a(5) = x3i
    a(6) = x5r
    a(7) = x5i
    a(8) = x1r
    a(9) = x1i
    a(10) = x6r
    a(11) = x6i
    a(12) = x2r
    a(13) = x2i
    a(14) = x4r
    a(15) = x4i
  end subroutine bitrv208neg
  !
  subroutine cftf1st(n, a, w)
    integer n, j, j0, j1, j2, j3, k, m, mh
    Real (kind=DP) :: a(0 : n - 1), w(0 : *)
    Real (kind=DP) :: wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i
    Real (kind=DP) :: wd1r, wd1i, wd3r, wd3i
    Real (kind=DP) :: x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    Real (kind=DP) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    mh = n / 8
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(0) + a(j2)
    x0i = a(1) + a(j2 + 1)
    x1r = a(0) - a(j2)
    x1i = a(1) - a(j2 + 1)
    x2r = a(j1) + a(j3)
    x2i = a(j1 + 1) + a(j3 + 1)
    x3r = a(j1) - a(j3)
    x3i = a(j1 + 1) - a(j3 + 1)
    a(0) = x0r + x2r
    a(1) = x0i + x2i
    a(j1) = x0r - x2r
    a(j1 + 1) = x0i - x2i
    a(j2) = x1r - x3i
    a(j2 + 1) = x1i + x3r
    a(j3) = x1r + x3i
    a(j3 + 1) = x1i - x3r
    wn4r = w(1)
    csc1 = w(2)
    csc3 = w(3)
    wd1r = 1
    wd1i = 0
    wd3r = 1
    wd3i = 0
    k = 0
    do j = 2, mh - 6, 4
       k = k + 4
       wk1r = csc1 * (wd1r + w(k))
       wk1i = csc1 * (wd1i + w(k + 1))
       wk3r = csc3 * (wd3r + w(k + 2))
       wk3i = csc3 * (wd3i + w(k + 3))
       wd1r = w(k)
       wd1i = w(k + 1)
       wd3r = w(k + 2)
       wd3i = w(k + 3)
       j1 = j + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j) + a(j2)
       x0i = a(j + 1) + a(j2 + 1)
       x1r = a(j) - a(j2)
       x1i = a(j + 1) - a(j2 + 1)
       y0r = a(j + 2) + a(j2 + 2)
       y0i = a(j + 3) + a(j2 + 3)
       y1r = a(j + 2) - a(j2 + 2)
       y1i = a(j + 3) - a(j2 + 3)
       x2r = a(j1) + a(j3)
       x2i = a(j1 + 1) + a(j3 + 1)
       x3r = a(j1) - a(j3)
       x3i = a(j1 + 1) - a(j3 + 1)
       y2r = a(j1 + 2) + a(j3 + 2)
       y2i = a(j1 + 3) + a(j3 + 3)
       y3r = a(j1 + 2) - a(j3 + 2)
       y3i = a(j1 + 3) - a(j3 + 3)
       a(j) = x0r + x2r
       a(j + 1) = x0i + x2i
       a(j + 2) = y0r + y2r
       a(j + 3) = y0i + y2i
       a(j1) = x0r - x2r
       a(j1 + 1) = x0i - x2i
       a(j1 + 2) = y0r - y2r
       a(j1 + 3) = y0i - y2i
       x0r = x1r - x3i
       x0i = x1i + x3r
       a(j2) = wk1r * x0r - wk1i * x0i
       a(j2 + 1) = wk1r * x0i + wk1i * x0r
       x0r = y1r - y3i
       x0i = y1i + y3r
       a(j2 + 2) = wd1r * x0r - wd1i * x0i
       a(j2 + 3) = wd1r * x0i + wd1i * x0r
       x0r = x1r + x3i
       x0i = x1i - x3r
       a(j3) = wk3r * x0r + wk3i * x0i
       a(j3 + 1) = wk3r * x0i - wk3i * x0r
       x0r = y1r + y3i
       x0i = y1i - y3r
       a(j3 + 2) = wd3r * x0r + wd3i * x0i
       a(j3 + 3) = wd3r * x0i - wd3i * x0r
       j0 = m - j
       j1 = j0 + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j0) + a(j2)
       x0i = a(j0 + 1) + a(j2 + 1)
       x1r = a(j0) - a(j2)
       x1i = a(j0 + 1) - a(j2 + 1)
       y0r = a(j0 - 2) + a(j2 - 2)
       y0i = a(j0 - 1) + a(j2 - 1)
       y1r = a(j0 - 2) - a(j2 - 2)
       y1i = a(j0 - 1) - a(j2 - 1)
       x2r = a(j1) + a(j3)
       x2i = a(j1 + 1) + a(j3 + 1)
       x3r = a(j1) - a(j3)
       x3i = a(j1 + 1) - a(j3 + 1)
       y2r = a(j1 - 2) + a(j3 - 2)
       y2i = a(j1 - 1) + a(j3 - 1)
       y3r = a(j1 - 2) - a(j3 - 2)
       y3i = a(j1 - 1) - a(j3 - 1)
       a(j0) = x0r + x2r
       a(j0 + 1) = x0i + x2i
       a(j0 - 2) = y0r + y2r
       a(j0 - 1) = y0i + y2i
       a(j1) = x0r - x2r
       a(j1 + 1) = x0i - x2i
       a(j1 - 2) = y0r - y2r
       a(j1 - 1) = y0i - y2i
       x0r = x1r - x3i
       x0i = x1i + x3r
       a(j2) = wk1i * x0r - wk1r * x0i
       a(j2 + 1) = wk1i * x0i + wk1r * x0r
       x0r = y1r - y3i
       x0i = y1i + y3r
       a(j2 - 2) = wd1i * x0r - wd1r * x0i
       a(j2 - 1) = wd1i * x0i + wd1r * x0r
       x0r = x1r + x3i
       x0i = x1i - x3r
       a(j3) = wk3i * x0r + wk3r * x0i
       a(j3 + 1) = wk3i * x0i - wk3r * x0r
       x0r = y1r + y3i
       x0i = y1i - y3r
       a(j3 - 2) = wd3i * x0r + wd3r * x0i
       a(j3 - 1) = wd3i * x0i - wd3r * x0r
    end do
    wk1r = csc1 * (wd1r + wn4r)
    wk1i = csc1 * (wd1i + wn4r)
    wk3r = csc3 * (wd3r - wn4r)
    wk3i = csc3 * (wd3i - wn4r)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(j0 - 2) + a(j2 - 2)
    x0i = a(j0 - 1) + a(j2 - 1)
    x1r = a(j0 - 2) - a(j2 - 2)
    x1i = a(j0 - 1) - a(j2 - 1)
    x2r = a(j1 - 2) + a(j3 - 2)
    x2i = a(j1 - 1) + a(j3 - 1)
    x3r = a(j1 - 2) - a(j3 - 2)
    x3i = a(j1 - 1) - a(j3 - 1)
    a(j0 - 2) = x0r + x2r
    a(j0 - 1) = x0i + x2i
    a(j1 - 2) = x0r - x2r
    a(j1 - 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(j2 - 2) = wk1r * x0r - wk1i * x0i
    a(j2 - 1) = wk1r * x0i + wk1i * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(j3 - 2) = wk3r * x0r + wk3i * x0i
    a(j3 - 1) = wk3r * x0i - wk3i * x0r
    x0r = a(j0) + a(j2)
    x0i = a(j0 + 1) + a(j2 + 1)
    x1r = a(j0) - a(j2)
    x1i = a(j0 + 1) - a(j2 + 1)
    x2r = a(j1) + a(j3)
    x2i = a(j1 + 1) + a(j3 + 1)
    x3r = a(j1) - a(j3)
    x3i = a(j1 + 1) - a(j3 + 1)
    a(j0) = x0r + x2r
    a(j0 + 1) = x0i + x2i
    a(j1) = x0r - x2r
    a(j1 + 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(j2) = wn4r * (x0r - x0i)
    a(j2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(j3) = -wn4r * (x0r + x0i)
    a(j3 + 1) = -wn4r * (x0i - x0r)
    x0r = a(j0 + 2) + a(j2 + 2)
    x0i = a(j0 + 3) + a(j2 + 3)
    x1r = a(j0 + 2) - a(j2 + 2)
    x1i = a(j0 + 3) - a(j2 + 3)
    x2r = a(j1 + 2) + a(j3 + 2)
    x2i = a(j1 + 3) + a(j3 + 3)
    x3r = a(j1 + 2) - a(j3 + 2)
    x3i = a(j1 + 3) - a(j3 + 3)
    a(j0 + 2) = x0r + x2r
    a(j0 + 3) = x0i + x2i
    a(j1 + 2) = x0r - x2r
    a(j1 + 3) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(j2 + 2) = wk1i * x0r - wk1r * x0i
    a(j2 + 3) = wk1i * x0i + wk1r * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(j3 + 2) = wk3i * x0r + wk3r * x0i
    a(j3 + 3) = wk3i * x0i - wk3r * x0r
  end subroutine cftf1st
  !
  subroutine cftb1st(n, a, w)
    integer n, j, j0, j1, j2, j3, k, m, mh
    Real (kind=DP) :: a(0 : n - 1), w(0 : *)
    Real (kind=DP) :: wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i
    Real (kind=DP) :: wd1r, wd1i, wd3r, wd3i
    Real (kind=DP) :: x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    Real (kind=DP) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    mh = n / 8
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(0) + a(j2)
    x0i = -a(1) - a(j2 + 1)
    x1r = a(0) - a(j2)
    x1i = -a(1) + a(j2 + 1)
    x2r = a(j1) + a(j3)
    x2i = a(j1 + 1) + a(j3 + 1)
    x3r = a(j1) - a(j3)
    x3i = a(j1 + 1) - a(j3 + 1)
    a(0) = x0r + x2r
    a(1) = x0i - x2i
    a(j1) = x0r - x2r
    a(j1 + 1) = x0i + x2i
    a(j2) = x1r + x3i
    a(j2 + 1) = x1i + x3r
    a(j3) = x1r - x3i
    a(j3 + 1) = x1i - x3r
    wn4r = w(1)
    csc1 = w(2)
    csc3 = w(3)
    wd1r = 1
    wd1i = 0
    wd3r = 1
    wd3i = 0
    k = 0
    do j = 2, mh - 6, 4
       k = k + 4
       wk1r = csc1 * (wd1r + w(k))
       wk1i = csc1 * (wd1i + w(k + 1))
       wk3r = csc3 * (wd3r + w(k + 2))
       wk3i = csc3 * (wd3i + w(k + 3))
       wd1r = w(k)
       wd1i = w(k + 1)
       wd3r = w(k + 2)
       wd3i = w(k + 3)
       j1 = j + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j) + a(j2)
       x0i = -a(j + 1) - a(j2 + 1)
       x1r = a(j) - a(j2)
       x1i = -a(j + 1) + a(j2 + 1)
       y0r = a(j + 2) + a(j2 + 2)
       y0i = -a(j + 3) - a(j2 + 3)
       y1r = a(j + 2) - a(j2 + 2)
       y1i = -a(j + 3) + a(j2 + 3)
       x2r = a(j1) + a(j3)
       x2i = a(j1 + 1) + a(j3 + 1)
       x3r = a(j1) - a(j3)
       x3i = a(j1 + 1) - a(j3 + 1)
       y2r = a(j1 + 2) + a(j3 + 2)
       y2i = a(j1 + 3) + a(j3 + 3)
       y3r = a(j1 + 2) - a(j3 + 2)
       y3i = a(j1 + 3) - a(j3 + 3)
       a(j) = x0r + x2r
       a(j + 1) = x0i - x2i
       a(j + 2) = y0r + y2r
       a(j + 3) = y0i - y2i
       a(j1) = x0r - x2r
       a(j1 + 1) = x0i + x2i
       a(j1 + 2) = y0r - y2r
       a(j1 + 3) = y0i + y2i
       x0r = x1r + x3i
       x0i = x1i + x3r
       a(j2) = wk1r * x0r - wk1i * x0i
       a(j2 + 1) = wk1r * x0i + wk1i * x0r
       x0r = y1r + y3i
       x0i = y1i + y3r
       a(j2 + 2) = wd1r * x0r - wd1i * x0i
       a(j2 + 3) = wd1r * x0i + wd1i * x0r
       x0r = x1r - x3i
       x0i = x1i - x3r
       a(j3) = wk3r * x0r + wk3i * x0i
       a(j3 + 1) = wk3r * x0i - wk3i * x0r
       x0r = y1r - y3i
       x0i = y1i - y3r
       a(j3 + 2) = wd3r * x0r + wd3i * x0i
       a(j3 + 3) = wd3r * x0i - wd3i * x0r
       j0 = m - j
       j1 = j0 + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j0) + a(j2)
       x0i = -a(j0 + 1) - a(j2 + 1)
       x1r = a(j0) - a(j2)
       x1i = -a(j0 + 1) + a(j2 + 1)
       y0r = a(j0 - 2) + a(j2 - 2)
       y0i = -a(j0 - 1) - a(j2 - 1)
       y1r = a(j0 - 2) - a(j2 - 2)
       y1i = -a(j0 - 1) + a(j2 - 1)
       x2r = a(j1) + a(j3)
       x2i = a(j1 + 1) + a(j3 + 1)
       x3r = a(j1) - a(j3)
       x3i = a(j1 + 1) - a(j3 + 1)
       y2r = a(j1 - 2) + a(j3 - 2)
       y2i = a(j1 - 1) + a(j3 - 1)
       y3r = a(j1 - 2) - a(j3 - 2)
       y3i = a(j1 - 1) - a(j3 - 1)
       a(j0) = x0r + x2r
       a(j0 + 1) = x0i - x2i
       a(j0 - 2) = y0r + y2r
       a(j0 - 1) = y0i - y2i
       a(j1) = x0r - x2r
       a(j1 + 1) = x0i + x2i
       a(j1 - 2) = y0r - y2r
       a(j1 - 1) = y0i + y2i
       x0r = x1r + x3i
       x0i = x1i + x3r
       a(j2) = wk1i * x0r - wk1r * x0i
       a(j2 + 1) = wk1i * x0i + wk1r * x0r
       x0r = y1r + y3i
       x0i = y1i + y3r
       a(j2 - 2) = wd1i * x0r - wd1r * x0i
       a(j2 - 1) = wd1i * x0i + wd1r * x0r
       x0r = x1r - x3i
       x0i = x1i - x3r
       a(j3) = wk3i * x0r + wk3r * x0i
       a(j3 + 1) = wk3i * x0i - wk3r * x0r
       x0r = y1r - y3i
       x0i = y1i - y3r
       a(j3 - 2) = wd3i * x0r + wd3r * x0i
       a(j3 - 1) = wd3i * x0i - wd3r * x0r
    end do
    wk1r = csc1 * (wd1r + wn4r)
    wk1i = csc1 * (wd1i + wn4r)
    wk3r = csc3 * (wd3r - wn4r)
    wk3i = csc3 * (wd3i - wn4r)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(j0 - 2) + a(j2 - 2)
    x0i = -a(j0 - 1) - a(j2 - 1)
    x1r = a(j0 - 2) - a(j2 - 2)
    x1i = -a(j0 - 1) + a(j2 - 1)
    x2r = a(j1 - 2) + a(j3 - 2)
    x2i = a(j1 - 1) + a(j3 - 1)
    x3r = a(j1 - 2) - a(j3 - 2)
    x3i = a(j1 - 1) - a(j3 - 1)
    a(j0 - 2) = x0r + x2r
    a(j0 - 1) = x0i - x2i
    a(j1 - 2) = x0r - x2r
    a(j1 - 1) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(j2 - 2) = wk1r * x0r - wk1i * x0i
    a(j2 - 1) = wk1r * x0i + wk1i * x0r
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(j3 - 2) = wk3r * x0r + wk3i * x0i
    a(j3 - 1) = wk3r * x0i - wk3i * x0r
    x0r = a(j0) + a(j2)
    x0i = -a(j0 + 1) - a(j2 + 1)
    x1r = a(j0) - a(j2)
    x1i = -a(j0 + 1) + a(j2 + 1)
    x2r = a(j1) + a(j3)
    x2i = a(j1 + 1) + a(j3 + 1)
    x3r = a(j1) - a(j3)
    x3i = a(j1 + 1) - a(j3 + 1)
    a(j0) = x0r + x2r
    a(j0 + 1) = x0i - x2i
    a(j1) = x0r - x2r
    a(j1 + 1) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(j2) = wn4r * (x0r - x0i)
    a(j2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(j3) = -wn4r * (x0r + x0i)
    a(j3 + 1) = -wn4r * (x0i - x0r)
    x0r = a(j0 + 2) + a(j2 + 2)
    x0i = -a(j0 + 3) - a(j2 + 3)
    x1r = a(j0 + 2) - a(j2 + 2)
    x1i = -a(j0 + 3) + a(j2 + 3)
    x2r = a(j1 + 2) + a(j3 + 2)
    x2i = a(j1 + 3) + a(j3 + 3)
    x3r = a(j1 + 2) - a(j3 + 2)
    x3i = a(j1 + 3) - a(j3 + 3)
    a(j0 + 2) = x0r + x2r
    a(j0 + 3) = x0i - x2i
    a(j1 + 2) = x0r - x2r
    a(j1 + 3) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(j2 + 2) = wk1i * x0r - wk1r * x0i
    a(j2 + 3) = wk1i * x0i + wk1r * x0r
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(j3 + 2) = wk3i * x0r + wk3r * x0i
    a(j3 + 3) = wk3i * x0i - wk3r * x0r
  end subroutine cftb1st
  !
  subroutine cftrec4(n, a, nw, w)
    integer n, nw, isplt, j, k, m
    Real (kind=DP) :: a(0 : n - 1), w(0 : nw - 1)
    m = n
    do while (m .gt. 512)
       m = m / 4
       call cftmdl1(m, a(n - m), w(nw - m / 2))
    end do
    call cftleaf(m, 1, a(n - m), nw, w)
    k = 0
    do j = n - m, m, -m
       k = k + 1
       isplt = cfttree(m, j, k, a, nw, w)
       call cftleaf(m, isplt, a(j - m), nw, w)
    end do
  end subroutine cftrec4
  !
  subroutine cftleaf(n, isplt, a, nw, w)
    integer n, isplt, nw
    Real (kind=DP) :: a(0 : n - 1), w(0 : nw - 1)
    if (n .eq. 512) then
       call cftmdl1(128, a, w(nw - 64))
       call cftf161(a, w(nw - 8))
       call cftf162(a(32), w(nw - 32))
       call cftf161(a(64), w(nw - 8))
       call cftf161(a(96), w(nw - 8))
       call cftmdl2(128, a(128), w(nw - 128))
       call cftf161(a(128), w(nw - 8))
       call cftf162(a(160), w(nw - 32))
       call cftf161(a(192), w(nw - 8))
       call cftf162(a(224), w(nw - 32))
       call cftmdl1(128, a(256), w(nw - 64))
       call cftf161(a(256), w(nw - 8))
       call cftf162(a(288), w(nw - 32))
       call cftf161(a(320), w(nw - 8))
       call cftf161(a(352), w(nw - 8))
       if (isplt .ne. 0) then
          call cftmdl1(128, a(384), w(nw - 64))
          call cftf161(a(480), w(nw - 8))
       else
          call cftmdl2(128, a(384), w(nw - 128))
          call cftf162(a(480), w(nw - 32))
       end if
       call cftf161(a(384), w(nw - 8))
       call cftf162(a(416), w(nw - 32))
       call cftf161(a(448), w(nw - 8))
    else
       call cftmdl1(64, a, w(nw - 32))
       call cftf081(a, w(nw - 8))
       call cftf082(a(16), w(nw - 8))
       call cftf081(a(32), w(nw - 8))
       call cftf081(a(48), w(nw - 8))
       call cftmdl2(64, a(64), w(nw - 64))
       call cftf081(a(64), w(nw - 8))
       call cftf082(a(80), w(nw - 8))
       call cftf081(a(96), w(nw - 8))
       call cftf082(a(112), w(nw - 8))
       call cftmdl1(64, a(128), w(nw - 32))
       call cftf081(a(128), w(nw - 8))
       call cftf082(a(144), w(nw - 8))
       call cftf081(a(160), w(nw - 8))
       call cftf081(a(176), w(nw - 8))
       if (isplt .ne. 0) then
          call cftmdl1(64, a(192), w(nw - 32))
          call cftf081(a(240), w(nw - 8))
       else
          call cftmdl2(64, a(192), w(nw - 64))
          call cftf082(a(240), w(nw - 8))
       end if
       call cftf081(a(192), w(nw - 8))
       call cftf082(a(208), w(nw - 8))
       call cftf081(a(224), w(nw - 8))
    end if
  end subroutine cftleaf
  !
  subroutine cftmdl1(n, a, w)
    integer n, j, j0, j1, j2, j3, k, m, mh
    Real (kind=DP) :: a(0 : n - 1), w(0 : *)
    Real (kind=DP) :: wn4r, wk1r, wk1i, wk3r, wk3i
    Real (kind=DP) :: x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    mh = n / 8
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(0) + a(j2)
    x0i = a(1) + a(j2 + 1)
    x1r = a(0) - a(j2)
    x1i = a(1) - a(j2 + 1)
    x2r = a(j1) + a(j3)
    x2i = a(j1 + 1) + a(j3 + 1)
    x3r = a(j1) - a(j3)
    x3i = a(j1 + 1) - a(j3 + 1)
    a(0) = x0r + x2r
    a(1) = x0i + x2i
    a(j1) = x0r - x2r
    a(j1 + 1) = x0i - x2i
    a(j2) = x1r - x3i
    a(j2 + 1) = x1i + x3r
    a(j3) = x1r + x3i
    a(j3 + 1) = x1i - x3r
    wn4r = w(1)
    k = 0
    do j = 2, mh - 2, 2
       k = k + 4
       wk1r = w(k)
       wk1i = w(k + 1)
       wk3r = w(k + 2)
       wk3i = w(k + 3)
       j1 = j + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j) + a(j2)
       x0i = a(j + 1) + a(j2 + 1)
       x1r = a(j) - a(j2)
       x1i = a(j + 1) - a(j2 + 1)
       x2r = a(j1) + a(j3)
       x2i = a(j1 + 1) + a(j3 + 1)
       x3r = a(j1) - a(j3)
       x3i = a(j1 + 1) - a(j3 + 1)
       a(j) = x0r + x2r
       a(j + 1) = x0i + x2i
       a(j1) = x0r - x2r
       a(j1 + 1) = x0i - x2i
       x0r = x1r - x3i
       x0i = x1i + x3r
       a(j2) = wk1r * x0r - wk1i * x0i
       a(j2 + 1) = wk1r * x0i + wk1i * x0r
       x0r = x1r + x3i
       x0i = x1i - x3r
       a(j3) = wk3r * x0r + wk3i * x0i
       a(j3 + 1) = wk3r * x0i - wk3i * x0r
       j0 = m - j
       j1 = j0 + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j0) + a(j2)
       x0i = a(j0 + 1) + a(j2 + 1)
       x1r = a(j0) - a(j2)
       x1i = a(j0 + 1) - a(j2 + 1)
       x2r = a(j1) + a(j3)
       x2i = a(j1 + 1) + a(j3 + 1)
       x3r = a(j1) - a(j3)
       x3i = a(j1 + 1) - a(j3 + 1)
       a(j0) = x0r + x2r
       a(j0 + 1) = x0i + x2i
       a(j1) = x0r - x2r
       a(j1 + 1) = x0i - x2i
       x0r = x1r - x3i
       x0i = x1i + x3r
       a(j2) = wk1i * x0r - wk1r * x0i
       a(j2 + 1) = wk1i * x0i + wk1r * x0r
       x0r = x1r + x3i
       x0i = x1i - x3r
       a(j3) = wk3i * x0r + wk3r * x0i
       a(j3 + 1) = wk3i * x0i - wk3r * x0r
    end do
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(j0) + a(j2)
    x0i = a(j0 + 1) + a(j2 + 1)
    x1r = a(j0) - a(j2)
    x1i = a(j0 + 1) - a(j2 + 1)
    x2r = a(j1) + a(j3)
    x2i = a(j1 + 1) + a(j3 + 1)
    x3r = a(j1) - a(j3)
    x3i = a(j1 + 1) - a(j3 + 1)
    a(j0) = x0r + x2r
    a(j0 + 1) = x0i + x2i
    a(j1) = x0r - x2r
    a(j1 + 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(j2) = wn4r * (x0r - x0i)
    a(j2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(j3) = -wn4r * (x0r + x0i)
    a(j3 + 1) = -wn4r * (x0i - x0r)
  end subroutine cftmdl1
  !
  subroutine cftmdl2(n, a, w)
    integer n, j, j0, j1, j2, j3, k, kr, m, mh
    Real (kind=DP) :: a(0 : n - 1), w(0 : *)
    Real (kind=DP) :: wn4r, wk1r, wk1i, wk3r, wk3i, wd1r, wd1i, wd3r, wd3i
    Real (kind=DP) :: x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    Real (kind=DP) :: y0r, y0i, y2r, y2i
    mh = n / 8
    m = 2 * mh
    wn4r = w(1)
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(0) - a(j2 + 1)
    x0i = a(1) + a(j2)
    x1r = a(0) + a(j2 + 1)
    x1i = a(1) - a(j2)
    x2r = a(j1) - a(j3 + 1)
    x2i = a(j1 + 1) + a(j3)
    x3r = a(j1) + a(j3 + 1)
    x3i = a(j1 + 1) - a(j3)
    y0r = wn4r * (x2r - x2i)
    y0i = wn4r * (x2i + x2r)
    a(0) = x0r + y0r
    a(1) = x0i + y0i
    a(j1) = x0r - y0r
    a(j1 + 1) = x0i - y0i
    y0r = wn4r * (x3r - x3i)
    y0i = wn4r * (x3i + x3r)
    a(j2) = x1r - y0i
    a(j2 + 1) = x1i + y0r
    a(j3) = x1r + y0i
    a(j3 + 1) = x1i - y0r
    k = 0
    kr = 2 * m
    do j = 2, mh - 2, 2
       k = k + 4
       wk1r = w(k)
       wk1i = w(k + 1)
       wk3r = w(k + 2)
       wk3i = w(k + 3)
       kr = kr - 4
       wd1i = w(kr)
       wd1r = w(kr + 1)
       wd3i = w(kr + 2)
       wd3r = w(kr + 3)
       j1 = j + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j) - a(j2 + 1)
       x0i = a(j + 1) + a(j2)
       x1r = a(j) + a(j2 + 1)
       x1i = a(j + 1) - a(j2)
       x2r = a(j1) - a(j3 + 1)
       x2i = a(j1 + 1) + a(j3)
       x3r = a(j1) + a(j3 + 1)
       x3i = a(j1 + 1) - a(j3)
       y0r = wk1r * x0r - wk1i * x0i
       y0i = wk1r * x0i + wk1i * x0r
       y2r = wd1r * x2r - wd1i * x2i
       y2i = wd1r * x2i + wd1i * x2r
       a(j) = y0r + y2r
       a(j + 1) = y0i + y2i
       a(j1) = y0r - y2r
       a(j1 + 1) = y0i - y2i
       y0r = wk3r * x1r + wk3i * x1i
       y0i = wk3r * x1i - wk3i * x1r
       y2r = wd3r * x3r + wd3i * x3i
       y2i = wd3r * x3i - wd3i * x3r
       a(j2) = y0r + y2r
       a(j2 + 1) = y0i + y2i
       a(j3) = y0r - y2r
       a(j3 + 1) = y0i - y2i
       j0 = m - j
       j1 = j0 + m
       j2 = j1 + m
       j3 = j2 + m
       x0r = a(j0) - a(j2 + 1)
       x0i = a(j0 + 1) + a(j2)
       x1r = a(j0) + a(j2 + 1)
       x1i = a(j0 + 1) - a(j2)
       x2r = a(j1) - a(j3 + 1)
       x2i = a(j1 + 1) + a(j3)
       x3r = a(j1) + a(j3 + 1)
       x3i = a(j1 + 1) - a(j3)
       y0r = wd1i * x0r - wd1r * x0i
       y0i = wd1i * x0i + wd1r * x0r
       y2r = wk1i * x2r - wk1r * x2i
       y2i = wk1i * x2i + wk1r * x2r
       a(j0) = y0r + y2r
       a(j0 + 1) = y0i + y2i
       a(j1) = y0r - y2r
       a(j1 + 1) = y0i - y2i
       y0r = wd3i * x1r + wd3r * x1i
       y0i = wd3i * x1i - wd3r * x1r
       y2r = wk3i * x3r + wk3r * x3i
       y2i = wk3i * x3i - wk3r * x3r
       a(j2) = y0r + y2r
       a(j2 + 1) = y0i + y2i
       a(j3) = y0r - y2r
       a(j3 + 1) = y0i - y2i
    end do
    wk1r = w(m)
    wk1i = w(m + 1)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    x0r = a(j0) - a(j2 + 1)
    x0i = a(j0 + 1) + a(j2)
    x1r = a(j0) + a(j2 + 1)
    x1i = a(j0 + 1) - a(j2)
    x2r = a(j1) - a(j3 + 1)
    x2i = a(j1 + 1) + a(j3)
    x3r = a(j1) + a(j3 + 1)
    x3i = a(j1 + 1) - a(j3)
    y0r = wk1r * x0r - wk1i * x0i
    y0i = wk1r * x0i + wk1i * x0r
    y2r = wk1i * x2r - wk1r * x2i
    y2i = wk1i * x2i + wk1r * x2r
    a(j0) = y0r + y2r
    a(j0 + 1) = y0i + y2i
    a(j1) = y0r - y2r
    a(j1 + 1) = y0i - y2i
    y0r = wk1i * x1r - wk1r * x1i
    y0i = wk1i * x1i + wk1r * x1r
    y2r = wk1r * x3r - wk1i * x3i
    y2i = wk1r * x3i + wk1i * x3r
    a(j2) = y0r - y2r
    a(j2 + 1) = y0i - y2i
    a(j3) = y0r + y2r
    a(j3 + 1) = y0i + y2i
  end subroutine cftmdl2
  !
  subroutine cftfx41(n, a, nw, w)
    integer n, nw
    Real (kind=DP) :: a(0 : n - 1), w(0 : nw - 1)
    if (n .eq. 128) then
       call cftf161(a, w(nw - 8))
       call cftf162(a(32), w(nw - 32))
       call cftf161(a(64), w(nw - 8))
       call cftf161(a(96), w(nw - 8))
    else
       call cftf081(a, w(nw - 8))
       call cftf082(a(16), w(nw - 8))
       call cftf081(a(32), w(nw - 8))
       call cftf081(a(48), w(nw - 8))
    end if
  end subroutine cftfx41
  !
  subroutine cftf161(a, w)
    Real (kind=DP) :: a(0 : 31), w(0 : *), wn4r, wk1r, wk1i
    Real (kind=DP) :: x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    Real (kind=DP) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    Real (kind=DP) :: y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
    Real (kind=DP) :: y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i
    Real (kind=DP) :: y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i
    wn4r = w(1)
    wk1r = w(2)
    wk1i = w(3)
    x0r = a(0) + a(16)
    x0i = a(1) + a(17)
    x1r = a(0) - a(16)
    x1i = a(1) - a(17)
    x2r = a(8) + a(24)
    x2i = a(9) + a(25)
    x3r = a(8) - a(24)
    x3i = a(9) - a(25)
    y0r = x0r + x2r
    y0i = x0i + x2i
    y4r = x0r - x2r
    y4i = x0i - x2i
    y8r = x1r - x3i
    y8i = x1i + x3r
    y12r = x1r + x3i
    y12i = x1i - x3r
    x0r = a(2) + a(18)
    x0i = a(3) + a(19)
    x1r = a(2) - a(18)
    x1i = a(3) - a(19)
    x2r = a(10) + a(26)
    x2i = a(11) + a(27)
    x3r = a(10) - a(26)
    x3i = a(11) - a(27)
    y1r = x0r + x2r
    y1i = x0i + x2i
    y5r = x0r - x2r
    y5i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y9r = wk1r * x0r - wk1i * x0i
    y9i = wk1r * x0i + wk1i * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    y13r = wk1i * x0r - wk1r * x0i
    y13i = wk1i * x0i + wk1r * x0r
    x0r = a(4) + a(20)
    x0i = a(5) + a(21)
    x1r = a(4) - a(20)
    x1i = a(5) - a(21)
    x2r = a(12) + a(28)
    x2i = a(13) + a(29)
    x3r = a(12) - a(28)
    x3i = a(13) - a(29)
    y2r = x0r + x2r
    y2i = x0i + x2i
    y6r = x0r - x2r
    y6i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y10r = wn4r * (x0r - x0i)
    y10i = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    y14r = wn4r * (x0r + x0i)
    y14i = wn4r * (x0i - x0r)
    x0r = a(6) + a(22)
    x0i = a(7) + a(23)
    x1r = a(6) - a(22)
    x1i = a(7) - a(23)
    x2r = a(14) + a(30)
    x2i = a(15) + a(31)
    x3r = a(14) - a(30)
    x3i = a(15) - a(31)
    y3r = x0r + x2r
    y3i = x0i + x2i
    y7r = x0r - x2r
    y7i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y11r = wk1i * x0r - wk1r * x0i
    y11i = wk1i * x0i + wk1r * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    y15r = wk1r * x0r - wk1i * x0i
    y15i = wk1r * x0i + wk1i * x0r
    x0r = y12r - y14r
    x0i = y12i - y14i
    x1r = y12r + y14r
    x1i = y12i + y14i
    x2r = y13r - y15r
    x2i = y13i - y15i
    x3r = y13r + y15r
    x3i = y13i + y15i
    a(24) = x0r + x2r
    a(25) = x0i + x2i
    a(26) = x0r - x2r
    a(27) = x0i - x2i
    a(28) = x1r - x3i
    a(29) = x1i + x3r
    a(30) = x1r + x3i
    a(31) = x1i - x3r
    x0r = y8r + y10r
    x0i = y8i + y10i
    x1r = y8r - y10r
    x1i = y8i - y10i
    x2r = y9r + y11r
    x2i = y9i + y11i
    x3r = y9r - y11r
    x3i = y9i - y11i
    a(16) = x0r + x2r
    a(17) = x0i + x2i
    a(18) = x0r - x2r
    a(19) = x0i - x2i
    a(20) = x1r - x3i
    a(21) = x1i + x3r
    a(22) = x1r + x3i
    a(23) = x1i - x3r
    x0r = y5r - y7i
    x0i = y5i + y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    x0r = y5r + y7i
    x0i = y5i - y7r
    x3r = wn4r * (x0r - x0i)
    x3i = wn4r * (x0i + x0r)
    x0r = y4r - y6i
    x0i = y4i + y6r
    x1r = y4r + y6i
    x1i = y4i - y6r
    a(8) = x0r + x2r
    a(9) = x0i + x2i
    a(10) = x0r - x2r
    a(11) = x0i - x2i
    a(12) = x1r - x3i
    a(13) = x1i + x3r
    a(14) = x1r + x3i
    a(15) = x1i - x3r
    x0r = y0r + y2r
    x0i = y0i + y2i
    x1r = y0r - y2r
    x1i = y0i - y2i
    x2r = y1r + y3r
    x2i = y1i + y3i
    x3r = y1r - y3r
    x3i = y1i - y3i
    a(0) = x0r + x2r
    a(1) = x0i + x2i
    a(2) = x0r - x2r
    a(3) = x0i - x2i
    a(4) = x1r - x3i
    a(5) = x1i + x3r
    a(6) = x1r + x3i
    a(7) = x1i - x3r
  end subroutine cftf161
  !
  subroutine cftf162(a, w)
    Real (kind=DP) :: a(0 : 31), w(0 : *)
    Real (kind=DP) :: wn4r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
    Real (kind=DP) :: x0r, x0i, x1r, x1i, x2r, x2i
    Real (kind=DP) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    Real (kind=DP) :: y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
    Real (kind=DP) :: y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i
    Real (kind=DP) :: y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i
    wn4r = w(1)
    wk1r = w(4)
    wk1i = w(5)
    wk3r = w(6)
    wk3i = -w(7)
    wk2r = w(8)
    wk2i = w(9)
    x1r = a(0) - a(17)
    x1i = a(1) + a(16)
    x0r = a(8) - a(25)
    x0i = a(9) + a(24)
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    y0r = x1r + x2r
    y0i = x1i + x2i
    y4r = x1r - x2r
    y4i = x1i - x2i
    x1r = a(0) + a(17)
    x1i = a(1) - a(16)
    x0r = a(8) + a(25)
    x0i = a(9) - a(24)
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    y8r = x1r - x2i
    y8i = x1i + x2r
    y12r = x1r + x2i
    y12i = x1i - x2r
    x0r = a(2) - a(19)
    x0i = a(3) + a(18)
    x1r = wk1r * x0r - wk1i * x0i
    x1i = wk1r * x0i + wk1i * x0r
    x0r = a(10) - a(27)
    x0i = a(11) + a(26)
    x2r = wk3i * x0r - wk3r * x0i
    x2i = wk3i * x0i + wk3r * x0r
    y1r = x1r + x2r
    y1i = x1i + x2i
    y5r = x1r - x2r
    y5i = x1i - x2i
    x0r = a(2) + a(19)
    x0i = a(3) - a(18)
    x1r = wk3r * x0r - wk3i * x0i
    x1i = wk3r * x0i + wk3i * x0r
    x0r = a(10) + a(27)
    x0i = a(11) - a(26)
    x2r = wk1r * x0r + wk1i * x0i
    x2i = wk1r * x0i - wk1i * x0r
    y9r = x1r - x2r
    y9i = x1i - x2i
    y13r = x1r + x2r
    y13i = x1i + x2i
    x0r = a(4) - a(21)
    x0i = a(5) + a(20)
    x1r = wk2r * x0r - wk2i * x0i
    x1i = wk2r * x0i + wk2i * x0r
    x0r = a(12) - a(29)
    x0i = a(13) + a(28)
    x2r = wk2i * x0r - wk2r * x0i
    x2i = wk2i * x0i + wk2r * x0r
    y2r = x1r + x2r
    y2i = x1i + x2i
    y6r = x1r - x2r
    y6i = x1i - x2i
    x0r = a(4) + a(21)
    x0i = a(5) - a(20)
    x1r = wk2i * x0r - wk2r * x0i
    x1i = wk2i * x0i + wk2r * x0r
    x0r = a(12) + a(29)
    x0i = a(13) - a(28)
    x2r = wk2r * x0r - wk2i * x0i
    x2i = wk2r * x0i + wk2i * x0r
    y10r = x1r - x2r
    y10i = x1i - x2i
    y14r = x1r + x2r
    y14i = x1i + x2i
    x0r = a(6) - a(23)
    x0i = a(7) + a(22)
    x1r = wk3r * x0r - wk3i * x0i
    x1i = wk3r * x0i + wk3i * x0r
    x0r = a(14) - a(31)
    x0i = a(15) + a(30)
    x2r = wk1i * x0r - wk1r * x0i
    x2i = wk1i * x0i + wk1r * x0r
    y3r = x1r + x2r
    y3i = x1i + x2i
    y7r = x1r - x2r
    y7i = x1i - x2i
    x0r = a(6) + a(23)
    x0i = a(7) - a(22)
    x1r = wk1i * x0r + wk1r * x0i
    x1i = wk1i * x0i - wk1r * x0r
    x0r = a(14) + a(31)
    x0i = a(15) - a(30)
    x2r = wk3i * x0r - wk3r * x0i
    x2i = wk3i * x0i + wk3r * x0r
    y11r = x1r + x2r
    y11i = x1i + x2i
    y15r = x1r - x2r
    y15i = x1i - x2i
    x1r = y0r + y2r
    x1i = y0i + y2i
    x2r = y1r + y3r
    x2i = y1i + y3i
    a(0) = x1r + x2r
    a(1) = x1i + x2i
    a(2) = x1r - x2r
    a(3) = x1i - x2i
    x1r = y0r - y2r
    x1i = y0i - y2i
    x2r = y1r - y3r
    x2i = y1i - y3i
    a(4) = x1r - x2i
    a(5) = x1i + x2r
    a(6) = x1r + x2i
    a(7) = x1i - x2r
    x1r = y4r - y6i
    x1i = y4i + y6r
    x0r = y5r - y7i
    x0i = y5i + y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(8) = x1r + x2r
    a(9) = x1i + x2i
    a(10) = x1r - x2r
    a(11) = x1i - x2i
    x1r = y4r + y6i
    x1i = y4i - y6r
    x0r = y5r + y7i
    x0i = y5i - y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(12) = x1r - x2i
    a(13) = x1i + x2r
    a(14) = x1r + x2i
    a(15) = x1i - x2r
    x1r = y8r + y10r
    x1i = y8i + y10i
    x2r = y9r - y11r
    x2i = y9i - y11i
    a(16) = x1r + x2r
    a(17) = x1i + x2i
    a(18) = x1r - x2r
    a(19) = x1i - x2i
    x1r = y8r - y10r
    x1i = y8i - y10i
    x2r = y9r + y11r
    x2i = y9i + y11i
    a(20) = x1r - x2i
    a(21) = x1i + x2r
    a(22) = x1r + x2i
    a(23) = x1i - x2r
    x1r = y12r - y14i
    x1i = y12i + y14r
    x0r = y13r + y15i
    x0i = y13i - y15r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(24) = x1r + x2r
    a(25) = x1i + x2i
    a(26) = x1r - x2r
    a(27) = x1i - x2i
    x1r = y12r + y14i
    x1i = y12i - y14r
    x0r = y13r - y15i
    x0i = y13i + y15r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(28) = x1r - x2i
    a(29) = x1i + x2r
    a(30) = x1r + x2i
    a(31) = x1i - x2r
  end subroutine cftf162
  !
  subroutine cftf081(a, w)
    Real (kind=DP) :: a(0 : 15), w(0 : *)
    Real (kind=DP) :: wn4r, x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    Real (kind=DP) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    Real (kind=DP) :: y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
    wn4r = w(1)
    x0r = a(0) + a(8)
    x0i = a(1) + a(9)
    x1r = a(0) - a(8)
    x1i = a(1) - a(9)
    x2r = a(4) + a(12)
    x2i = a(5) + a(13)
    x3r = a(4) - a(12)
    x3i = a(5) - a(13)
    y0r = x0r + x2r
    y0i = x0i + x2i
    y2r = x0r - x2r
    y2i = x0i - x2i
    y1r = x1r - x3i
    y1i = x1i + x3r
    y3r = x1r + x3i
    y3i = x1i - x3r
    x0r = a(2) + a(10)
    x0i = a(3) + a(11)
    x1r = a(2) - a(10)
    x1i = a(3) - a(11)
    x2r = a(6) + a(14)
    x2i = a(7) + a(15)
    x3r = a(6) - a(14)
    x3i = a(7) - a(15)
    y4r = x0r + x2r
    y4i = x0i + x2i
    y6r = x0r - x2r
    y6i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    x2r = x1r + x3i
    x2i = x1i - x3r
    y5r = wn4r * (x0r - x0i)
    y5i = wn4r * (x0r + x0i)
    y7r = wn4r * (x2r - x2i)
    y7i = wn4r * (x2r + x2i)
    a(8) = y1r + y5r
    a(9) = y1i + y5i
    a(10) = y1r - y5r
    a(11) = y1i - y5i
    a(12) = y3r - y7i
    a(13) = y3i + y7r
    a(14) = y3r + y7i
    a(15) = y3i - y7r
    a(0) = y0r + y4r
    a(1) = y0i + y4i
    a(2) = y0r - y4r
    a(3) = y0i - y4i
    a(4) = y2r - y6i
    a(5) = y2i + y6r
    a(6) = y2r + y6i
    a(7) = y2i - y6r
  end subroutine cftf081
  !
  subroutine cftf082(a, w)
    Real (kind=DP) :: a(0 : 15), w(0 : *)
    Real (kind=DP) :: wn4r, wk1r, wk1i, x0r, x0i, x1r, x1i
    Real (kind=DP) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    Real (kind=DP) :: y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
    wn4r = w(1)
    wk1r = w(2)
    wk1i = w(3)
    y0r = a(0) - a(9)
    y0i = a(1) + a(8)
    y1r = a(0) + a(9)
    y1i = a(1) - a(8)
    x0r = a(4) - a(13)
    x0i = a(5) + a(12)
    y2r = wn4r * (x0r - x0i)
    y2i = wn4r * (x0i + x0r)
    x0r = a(4) + a(13)
    x0i = a(5) - a(12)
    y3r = wn4r * (x0r - x0i)
    y3i = wn4r * (x0i + x0r)
    x0r = a(2) - a(11)
    x0i = a(3) + a(10)
    y4r = wk1r * x0r - wk1i * x0i
    y4i = wk1r * x0i + wk1i * x0r
    x0r = a(2) + a(11)
    x0i = a(3) - a(10)
    y5r = wk1i * x0r - wk1r * x0i
    y5i = wk1i * x0i + wk1r * x0r
    x0r = a(6) - a(15)
    x0i = a(7) + a(14)
    y6r = wk1i * x0r - wk1r * x0i
    y6i = wk1i * x0i + wk1r * x0r
    x0r = a(6) + a(15)
    x0i = a(7) - a(14)
    y7r = wk1r * x0r - wk1i * x0i
    y7i = wk1r * x0i + wk1i * x0r
    x0r = y0r + y2r
    x0i = y0i + y2i
    x1r = y4r + y6r
    x1i = y4i + y6i
    a(0) = x0r + x1r
    a(1) = x0i + x1i
    a(2) = x0r - x1r
    a(3) = x0i - x1i
    x0r = y0r - y2r
    x0i = y0i - y2i
    x1r = y4r - y6r
    x1i = y4i - y6i
    a(4) = x0r - x1i
    a(5) = x0i + x1r
    a(6) = x0r + x1i
    a(7) = x0i - x1r
    x0r = y1r - y3i
    x0i = y1i + y3r
    x1r = y5r - y7r
    x1i = y5i - y7i
    a(8) = x0r + x1r
    a(9) = x0i + x1i
    a(10) = x0r - x1r
    a(11) = x0i - x1i
    x0r = y1r + y3i
    x0i = y1i - y3r
    x1r = y5r + y7r
    x1i = y5i + y7i
    a(12) = x0r - x1i
    a(13) = x0i + x1r
    a(14) = x0r + x1i
    a(15) = x0i - x1r
  end subroutine cftf082
  !
  subroutine cftf040(a)
    Real (kind=DP) :: a(0 : 7), x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    x0r = a(0) + a(4)
    x0i = a(1) + a(5)
    x1r = a(0) - a(4)
    x1i = a(1) - a(5)
    x2r = a(2) + a(6)
    x2i = a(3) + a(7)
    x3r = a(2) - a(6)
    x3i = a(3) - a(7)
    a(0) = x0r + x2r
    a(1) = x0i + x2i
    a(2) = x1r - x3i
    a(3) = x1i + x3r
    a(4) = x0r - x2r
    a(5) = x0i - x2i
    a(6) = x1r + x3i
    a(7) = x1i - x3r
  end subroutine cftf040
  !
  subroutine cftb040(a)
    Real (kind=DP) :: a(0 : 7), x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    x0r = a(0) + a(4)
    x0i = a(1) + a(5)
    x1r = a(0) - a(4)
    x1i = a(1) - a(5)
    x2r = a(2) + a(6)
    x2i = a(3) + a(7)
    x3r = a(2) - a(6)
    x3i = a(3) - a(7)
    a(0) = x0r + x2r
    a(1) = x0i + x2i
    a(2) = x1r + x3i
    a(3) = x1i - x3r
    a(4) = x0r - x2r
    a(5) = x0i - x2i
    a(6) = x1r - x3i
    a(7) = x1i + x3r
  end subroutine cftb040
  !
  subroutine cftx020(a)
    Real (kind=DP) :: a(0 : 3), x0r, x0i
    x0r = a(0) - a(2)
    x0i = a(1) - a(3)
    a(0) = a(0) + a(2)
    a(1) = a(1) + a(3)
    a(2) = x0r
    a(3) = x0i
  end subroutine cftx020
  !
  subroutine rftfsub(n, a, nc, c)
    integer n, nc, j, k, kk, ks, m
    Real (kind=DP) :: a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
    m = n / 2
    ks = 2 * nc / m
    kk = 0
    do j = 2, m - 2, 2
       k = n - j
       kk = kk + ks
       wkr = 0.5d0 - c(nc - kk)
       wki = c(kk)
       xr = a(j) - a(k)
       xi = a(j + 1) + a(k + 1)
       yr = wkr * xr - wki * xi
       yi = wkr * xi + wki * xr
       a(j) = a(j) - yr
       a(j + 1) = a(j + 1) - yi
       a(k) = a(k) + yr
       a(k + 1) = a(k + 1) - yi
    end do
  end subroutine rftfsub
  !
  subroutine rftbsub(n, a, nc, c)
    integer n, nc, j, k, kk, ks, m
    Real (kind=DP) :: a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
    m = n / 2
    ks = 2 * nc / m
    kk = 0
    do j = 2, m - 2, 2
       k = n - j
       kk = kk + ks
       wkr = 0.5d0 - c(nc - kk)
       wki = c(kk)
       xr = a(j) - a(k)
       xi = a(j + 1) + a(k + 1)
       yr = wkr * xr + wki * xi
       yi = wkr * xi - wki * xr
       a(j) = a(j) - yr
       a(j + 1) = a(j + 1) - yi
       a(k) = a(k) + yr
       a(k + 1) = a(k + 1) - yi
    end do
  end subroutine rftbsub
  !
  subroutine dctsub(n, a, nc, c)
    integer n, nc, j, k, kk, ks, m
    Real (kind=DP) :: a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
    m = n / 2
    ks = nc / n
    kk = 0
    do j = 1, m - 1
       k = n - j
       kk = kk + ks
       wkr = c(kk) - c(nc - kk)
       wki = c(kk) + c(nc - kk)
       xr = wki * a(j) - wkr * a(k)
       a(j) = wkr * a(j) + wki * a(k)
       a(k) = xr
    end do
    a(m) = c(0) * a(m)
  end subroutine dctsub
  !
  subroutine dstsub(n, a, nc, c)
    integer n, nc, j, k, kk, ks, m
    Real (kind=DP) :: a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
    m = n / 2
    ks = nc / n
    kk = 0
    do j = 1, m - 1
       k = n - j
       kk = kk + ks
       wkr = c(kk) - c(nc - kk)
       wki = c(kk) + c(nc - kk)
       xr = wki * a(k) - wkr * a(j)
       a(k) = wkr * a(k) + wki * a(j)
       a(j) = xr
    end do
    a(m) = c(0) * a(m)
  end subroutine dstsub
  !


End MODULE Fourier

