!
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
  Recursive Function FastFT_1D(Data, Isign) Result (FFT)
! *
! ********************************
! * Calculates the Discrete Fourier
! * Transformation from the complex
! * data stored in data.
! ********************************

    Complex (kind=DPC), Intent (in) :: Data(:)
    Integer, Intent (in), Optional :: Isign
    Complex (kind=DPC) :: FFT(Size(Data))

    Complex (Kind=DPC) :: Factor, Prefactor
    Complex (kind=DPC), Allocatable :: DfE(:), DfO(:), TW(:)
    Integer :: Side, I, J, Ntot, Nterm, Nhalf


    If (.not. Present(Isign)) Then
       Side = -1
    Else
       Side = Isign
    End If

    Ntot = Size(Data)
    PreFactor = exp(Side*DPI*UnitImag_DPC / Cmplx(Ntot,kind=DP) )
    Write(*,*)Ntot, Prefactor
    Write(*,*)Data
    Write(*,*)

    If (Ntot == 2) Then
       FFT(1) = Data(2) - Data(1)
       FFT(2) = Data(2) + Data(1)
    Else If (Mod(Ntot,2) == 0) Then
       Nhalf = Int(Ntot/2)
       ALLOCATE(DfE(Nhalf), DfO(Nhalf), Tw(Nhalf))
       DfE = FastFT(Data(2:Ntot:2))
       DfO = FastFT(Data(1:Ntot:2))
       
       Tw(1) = Prefactor
       Do I = 2, Nhalf
          Tw(I) = Tw(I-1)*Prefactor
       End Do
       Do I = 1, Nhalf
          FFT(I)       = DfE(I) + TW(I)*DfO(I)
          FFT(I+Nhalf) = DfE(I) - TW(I)*DfO(I)
       End Do
       Deallocate(DfE, DfO, Tw)
    Else
       Do I = 1, Ntot
          Factor = PreFactor ** I
          
          FFT(I) = (0.0_DP, 0.0_DP)
          Do J = Ntot, 1, -1
             FFT(I) = ( FFT(I) + Data(J) ) * Factor
          End Do
       End Do
    End If

    If (Side == -1) Then
!       FFT(:) = FFT(:) / Cmplx(Ntot, kind=DP)
    End If

    Write(*,*)FFT
    Write(*,*)
    Write(*,*)

    Return
  End Function FASTFT_1D


End MODULE Fourier

