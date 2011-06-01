!
! MODULE to solve linear systems of equations
!
! Copyright (C) 2006  Alberto Ramos <alberto@martin.ft.uam.es>
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

! $ v. 1.0; Released: 15/09/2006; $

! ***********************************************************
! *
MODULE Linear
! *
! ***********************************************************

  USE NumTypes
  USE Error

  IMPLICIT NONE

  Interface Pivoting
     Module Procedure Pivoting_SP, Pivoting_DP, Pivoting_SPC, &
          & Pivoting_DPC
  End Interface

  Interface lu
     Module Procedure lu_SP, lu_DP, lu_SPC, lu_DPC
  End Interface

  Interface Cholesky
     Module Procedure Cholesky_DP, Cholesky_SP, Cholesky_DPC, Cholesky_SPC
  End Interface

  Interface lusolve
     Module Procedure lusolve_SP, lusolve_DP, lusolve_SPC, &
          & lusolve_DPC
  End Interface

  Interface Det
     Module Procedure Det_SP, Det_DP, Det_SPC, Det_DPC
  End Interface

  
  Private Pivoting_SP, lu_SP, lusolve_SP, Det_SP, &
       & Pivoting_DP, lu_DP, lusolve_DP, Det_DP, &
       & Pivoting_SPC, lu_SPC, lusolve_SPC, Det_SPC, &
       & Pivoting_DPC, lu_DPC, lusolve_DPC, Det_DPC, &
       & Cholesky_DP, Cholesky_SP, Cholesky_DPC, Cholesky_SPC

CONTAINS

!  *********************************************
!  *                                           *
  Subroutine Pivoting_SP(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Pivote M to arrange the elemnts in the
!  * diag. big.
!  *********************************************

    Real (kind=SP), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Real (kind=SP) :: Rval(Size(M,2))
    Integer :: Ipos(1), Kval, I

    ! Set  initial Idet, Ipiv
    Idet = 1
    forall (I=1:Size(Ipiv)) Ipiv(I) = I

    Do I = 1, Size(M,1) - 1
       Ipos = MaxLoc(Abs(M(I:,I))) + (I-1)

       Rval(:) = M(I,:)
       M(I,:) = M(Ipos(1),:)
       M(Ipos(1),:) = Rval(:)

       Kval = Ipiv(I)
       Ipiv(I) = Ipiv(Ipos(1))
       Ipiv(Ipos(1)) = Kval

       If (Ipiv(I) .ne. I) Idet = -Idet
    End Do

    Return
  End Subroutine Pivoting_SP

!  *********************************************
!  *                                           *
  Subroutine lu_SP(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Makes LU simple precision decomposition of 
!  * matrix M, with pivoting Ipiv. It returns  
!  * in Idet if the number of premutations is 
!  * even or odd.
!  *********************************************

    Real (kind=SP), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Integer :: I, J, Idim

    Idim = Size(M,1)
    ! First make pivoting
    CALL Pivoting(M, Ipiv, Idet)

    
    ! NOW: LU Decomposition
    ! The first step is done apart
    Do J = 2, Idim
       M(J, 1) = M(J, 1) / M(1, 1)
    End Do

    Do I = 2, Idim
       M(I, I) = M(I, I) - &
            & Dot_Product(M(I,1:I-1),M(1:I-1, I)) 

       Do J = i+1, Idim
          M(I, J) = M(I, J) - Dot_Product(M(I, 1:I-1),&
               & M(1:I-1, J))
          M(J, I) = M(J, I) - Dot_Product(M(J, 1:I-1)&
               &,M(1:I-1, I))

          If (Abs(M(I,I)) < Epsilon(1.0_SP)) &
               &CALL Abort('lu_SP','Singular matrix.')
          M(J, I) = M(J, I) / M(I, I)
       End Do
    End Do

    Return
  End Subroutine lu_SP

!  *********************************************
!  *                                           *
  Subroutine lusolve_SP(M, b)
!  *                                           *
!  *********************************************
!  * Solve a linear set of equations using LU 
!  * decomposition. Both M and b are overwritten
!  *********************************************

    Real (kind=SP), Intent (inout) :: M(:,:), b(:)

    Real (kind=SP) :: bcp(Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, Idim, I
    

    Idim = Size(M,1)
    CALL LU(M, Ipiv, Id)

    bcp = b
    ! Now we Permutate b
    Do I = 1, Idim
       b(I) = bcp(Ipiv(I))
    End Do

    ! First solve Lx = b
    Do I = 2, Idim
       b(I) = b(I) - Dot_Product(M(I,1:I-1),b(1:I-1))
    End Do
    
    ! Now solve Ux = b
    b(Idim) = b(Idim) / M(Idim, Idim)
    Do I = Idim - 1, 1, -1
       b(I) = b(I) - Dot_Product(M(I, I+1:Idim),b(I+1:Idim))
       b(I) = b(I) / M(I, I)
    End Do

    Return
  End Subroutine lusolve_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function Det_SP(M)
!  *                                           *
!  *********************************************
!  * Compute the determinant of the matrix M.
!  *********************************************

    Real (kind=SP), Intent (in) :: M(:,:)
    
    Real (kind=SP) :: Mcp(Size(M,1),Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, I, J, Idim

    Idim = Size(M,1)
    Mcp = M
    ! First make pivoting
    CALL Pivoting(Mcp, Ipiv, Id)

    
    ! NOW: LU Decomposition
    ! Separamos el paso I = 1
    Do J = 2, Idim
       Mcp(J, 1) = Mcp(J, 1) / Mcp(1, 1)
    End Do

    Det_SP = Mcp(1,1)
    Do I = 2, Idim
       Mcp(I, I) = Mcp(I, I) - &
            & Dot_Product(Mcp(I,1:I-1), Mcp(1:I-1, I)) 

       Do J = i+1, Idim
          Mcp(I, J) = Mcp(I, J) - Dot_Product(Mcp(I, 1:I-1),&
               & Mcp(1:I-1, J))
          Mcp(J, I) = Mcp(J, I) - Dot_Product(Mcp(J, 1:I-1)&
               &,Mcp(1:I-1, I))

          If (Abs(Mcp(I,I)) < Epsilon(1.0_SP)) Then
             Det_SP = 0.0_SP
             Return
          End If
          Mcp(J, I) = Mcp(J, I) / Mcp(I, I)
       End Do
       Det_SP = Det_SP * Mcp(I,I)
    End Do

    Det_SP = Real(Id,kind=SP)*Det_SP

    Return
  End Function Det_SP

!  *********************************************
!  *                                           *
  Subroutine Pivoting_DP(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Pivote M to arrange the elemnts in the
!  * diag. big.
!  *********************************************

    Real (kind=DP), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Real (kind=DP) :: Rval(Size(M,2))
    Integer :: Ipos(1), Kval, I

    ! Set  initial Idet, Ipiv
    Idet = 1
    forall (I=1:Size(Ipiv)) Ipiv(I) = I

    Do I = 1, Size(M,1) - 1
       Ipos = MaxLoc(Abs(M(I:,I))) + (I-1)

       Rval(:) = M(I,:)
       M(I,:) = M(Ipos(1),:)
       M(Ipos(1),:) = Rval(:)

       Kval = Ipiv(I)
       Ipiv(I) = Ipiv(Ipos(1))
       Ipiv(Ipos(1)) = Kval

       If (Ipiv(I) .ne. I) Idet = -Idet
    End Do

    Return
  End Subroutine Pivoting_DP

!  *********************************************
!  *                                           *
  Subroutine lu_DP(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Makes LU simple precision decomposition of 
!  * matrix M, with pivoting Ipiv. It returns  
!  * in Idet if the number of premutations is 
!  * even or odd.
!  *********************************************

    Real (kind=DP), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Integer :: I, J, Idim

    Idim = Size(M,1)
    ! First make pivoting
    CALL Pivoting(M, Ipiv, Idet)

    
    ! NOW: LU Decomposition
    ! Separamos el paso I = 1
    Do J = 2, Idim
       M(J, 1) = M(J, 1) / M(1, 1)
    End Do

    Do I = 2, Idim
       M(I, I) = M(I, I) - &
            & Dot_Product(M(I,1:I-1),M(1:I-1, I)) 

       Do J = i+1, Idim
          M(I, J) = M(I, J) - Dot_Product(M(I, 1:I-1),&
               & M(1:I-1, J))
          M(J, I) = M(J, I) - Dot_Product(M(J, 1:I-1)&
               &,M(1:I-1, I))

          If (Abs(M(I,I)) < Epsilon(1.0_DP)) &
               & CALL Abort('lu_DP','Singular matrix.')
          M(J, I) = M(J, I) / M(I, I)
       End Do
    End Do

    Return
  End Subroutine lu_DP

!  *********************************************
!  *                                           *
  Subroutine lusolve_DP(M, b)
!  *                                           *
!  *********************************************
!  * Solve a linear set of equations using LU 
!  * decomposition. Both M and b are overwritten
!  *********************************************

    Real (kind=DP), Intent (inout) :: M(:,:), b(:)

    Real (kind=DP) :: bcp(Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, Idim, I
    

    Idim = Size(M,1)
    CALL LU(M, Ipiv, Id)

    bcp = b
    ! Now we Permutate b
    Do I = 1, Idim
       b(I) = bcp(Ipiv(I))
    End Do

    ! First solve Lx = b
    Do I = 2, Idim
       b(I) = b(I) - Dot_Product(M(I,1:I-1),b(1:I-1))
    End Do
    
    ! Now solve Ux = b
    b(Idim) = b(Idim) / M(Idim, Idim)
    Do I = Idim - 1, 1, -1
       b(I) = b(I) - Dot_Product(M(I, I+1:Idim),b(I+1:Idim))
       b(I) = b(I) / M(I, I)
    End Do

    Return
  End Subroutine lusolve_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Det_DP(M)
!  *                                           *
!  *********************************************
!  * Compute the determinant of the matrix M.
!  *********************************************

    Real (kind=DP), Intent (in) :: M(:,:)
    
    Real (kind=DP) :: Mcp(Size(M,1),Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, I, J, Idim

    Idim = Size(M,1)
    Mcp = M
    ! First make pivoting
    CALL Pivoting(Mcp, Ipiv, Id)

    
    ! NOW: LU Decomposition
    ! Separamos el paso I = 1
    Do J = 2, Idim
       Mcp(J, 1) = Mcp(J, 1) / Mcp(1, 1)
    End Do

    Det_DP = Mcp(1,1)
    Do I = 2, Idim
       Mcp(I, I) = Mcp(I, I) - &
            & Dot_Product(Mcp(I,1:I-1), Mcp(1:I-1, I)) 

       Do J = i+1, Idim
          Mcp(I, J) = Mcp(I, J) - Dot_Product(Mcp(I, 1:I-1),&
               & Mcp(1:I-1, J))
          Mcp(J, I) = Mcp(J, I) - Dot_Product(Mcp(J, 1:I-1)&
               &,Mcp(1:I-1, I))

          If (Abs(Mcp(I,I)) < Epsilon(1.0_DP)) Then
             Det_DP = 0.0_DP
             Return
          End If
          Mcp(J, I) = Mcp(J, I) / Mcp(I, I)
       End Do
       Det_DP = Det_DP * Mcp(I,I)
    End Do

    Det_DP = Real(Id,kind=DP)*Det_DP

    Return
  End Function Det_DP

!  *********************************************
!  *                                           *
  Subroutine Pivoting_SPC(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Pivote M to arrange the elemnts in the
!  * diag. big.
!  *********************************************

    Complex (kind=SPC), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Complex (kind=SPC) :: Rval(Size(M,2))
    Integer :: Ipos(1), Kval, I

    ! Set  initial Idet, Ipiv
    Idet = 1
    forall (I=1:Size(Ipiv)) Ipiv(I) = I

    Do I = 1, Size(M,1) - 1
       Ipos = MaxLoc(Abs(M(I:,I))) + (I-1)

       Rval(:) = M(I,:)
       M(I,:) = M(Ipos(1),:)
       M(Ipos(1),:) = Rval(:)

       Kval = Ipiv(I)
       Ipiv(I) = Ipiv(Ipos(1))
       Ipiv(Ipos(1)) = Kval

       If (Ipiv(I) .ne. I) Idet = -Idet
    End Do

    Return
  End Subroutine Pivoting_SPC

!  *********************************************
!  *                                           *
  Subroutine lu_SPC(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Makes LU simple precision decomposition of 
!  * matrix M, with pivoting Ipiv. It returns  
!  * in Idet if the number of premutations is 
!  * even or odd.
!  *********************************************

    Complex (kind=SPC), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Integer :: I, J, Idim

    Idim = Size(M,1)
    ! First make pivoting
    CALL Pivoting(M, Ipiv, Idet)

    
    ! NOW: LU Decomposition
    ! The first step is done apart
    Do J = 2, Idim
       M(J, 1) = M(J, 1) / M(1, 1)
    End Do

    Do I = 2, Idim
       M(I, I) = M(I, I) - &
            & Dot_Product(Conjg(M(I,1:I-1)),M(1:I-1, I)) 

       Do J = i+1, Idim
          M(I, J) = M(I, J) - Dot_Product(Conjg(M(I, 1:I-1)),&
               & M(1:I-1, J))
          M(J, I) = M(J, I) - Dot_Product(Conjg(M(J, 1:I-1)),&
               & M(1:I-1, I))

          If (Abs(M(I,I)) < Epsilon(1.0_SP)) &
               &CALL Abort('lu_SP','Singular matrix.')
          M(J, I) = M(J, I) / M(I, I)
       End Do
    End Do

    Return
  End Subroutine lu_SPC

!  *********************************************
!  *                                           *
  Subroutine lusolve_SPC(M, b)
!  *                                           *
!  *********************************************
!  * Solve a linear set of equations using LU 
!  * decomposition. Both M and b are overwritten
!  *********************************************

    Complex (kind=SPC), Intent (inout) :: M(:,:), b(:)

    Complex (kind=SPC) :: bcp(Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, Idim, I
    

    Idim = Size(M,1)
    CALL LU(M, Ipiv, Id)

    bcp = b
    ! Now we Permutate b
    Do I = 1, Idim
       b(I) = bcp(Ipiv(I))
    End Do

    ! First solve Lx = b
    Do I = 2, Idim
       b(I) = b(I) - Dot_Product(Conjg(M(I,1:I-1)),b(1:I-1))
    End Do
    
    ! Now solve Ux = b
    b(Idim) = b(Idim) / M(Idim, Idim)
    Do I = Idim - 1, 1, -1
       b(I) = b(I) - Dot_Product(Conjg(M(I, I+1:Idim)),b(I+1:Idim))
       b(I) = b(I) / M(I, I)
    End Do

    Return
  End Subroutine lusolve_SPC

!  *********************************************
!  *                                           *
  Complex (kind=SPC) Function Det_SPC(M)
!  *                                           *
!  *********************************************
!  * Compute the determinant of the matrix M.
!  *********************************************

    Complex (kind=SPC), Intent (in) :: M(:,:)
    
    Complex (kind=SPC) :: Mcp(Size(M,1),Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, I, J, Idim

    Idim = Size(M,1)
    Mcp = M
    ! First make pivoting
    CALL Pivoting(Mcp, Ipiv, Id)

    
    ! NOW: LU Decomposition
    ! Separamos el paso I = 1
    Do J = 2, Idim
       Mcp(J, 1) = Mcp(J, 1) / Mcp(1, 1)
    End Do

    Det_SPC = Mcp(1,1)
    Do I = 2, Idim
       Mcp(I, I) = Mcp(I, I) - &
            & Dot_Product(Conjg(Mcp(I,1:I-1)), Mcp(1:I-1, I)) 

       Do J = i+1, Idim
          Mcp(I, J) = Mcp(I, J) - Dot_Product(Conjg(Mcp(I, 1:I-1)),&
               & Mcp(1:I-1, J))
          Mcp(J, I) = Mcp(J, I) - Dot_Product(Conjg(Mcp(J, 1:I-1)),&
               & Mcp(1:I-1, I))

          If (Abs(Mcp(I,I)) < Epsilon(1.0_SP)) Then
             Det_SPC = 0.0_SP
             Return
          End If
          Mcp(J, I) = Mcp(J, I) / Mcp(I, I)
       End Do
       Det_SPC = Det_SPC * Mcp(I,I)
    End Do

    Det_SPC = Cmplx(Real(Id,kind=SP),0.0_SP)*Det_SPC

    Return
  End Function Det_SPC

!  *********************************************
!  *                                           *
  Subroutine Pivoting_DPC(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Pivote M to arrange the elemnts in the
!  * diag. big.
!  *********************************************

    Complex (kind=DPC), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Complex (kind=DPC) :: Rval(Size(M,2))
    Integer :: Ipos(1), Kval, I

    ! Set  initial Idet, Ipiv
    Idet = 1
    forall (I=1:Size(Ipiv)) Ipiv(I) = I

    Do I = 1, Size(M,1) - 1
       Ipos = MaxLoc(Abs(M(I:,I))) + (I-1)

       Rval(:) = M(I,:)
       M(I,:) = M(Ipos(1),:)
       M(Ipos(1),:) = Rval(:)

       Kval = Ipiv(I)
       Ipiv(I) = Ipiv(Ipos(1))
       Ipiv(Ipos(1)) = Kval

       If (Ipiv(I) .ne. I) Idet = -Idet
    End Do

    Return
  End Subroutine Pivoting_DPC

!  *********************************************
!  *                                           *
  Subroutine lu_DPC(M,Ipiv,Idet)
!  *                                           *
!  *********************************************
!  * Makes LU simple precision decomposition of 
!  * matrix M, with pivoting Ipiv. It returns  
!  * in Idet if the number of premutations is 
!  * even or odd.
!  *********************************************

    Complex (kind=DPC), Intent (inout) :: M(:,:)
    Integer, Intent (out) :: Ipiv(:), Idet

    Integer :: I, J, Idim

    Idim = Size(M,1)
    ! First make pivoting
    CALL Pivoting(M, Ipiv, Idet)

    
    ! NOW: LU Decomposition
    ! The first step is done apart
    Do J = 2, Idim
       M(J, 1) = M(J, 1) / M(1, 1)
    End Do

    Do I = 2, Idim
       M(I, I) = M(I, I) - &
            & Dot_Product(Conjg(M(I,1:I-1)),M(1:I-1, I)) 

       Do J = i+1, Idim
          M(I, J) = M(I, J) - Dot_Product(Conjg(M(I, 1:I-1)),&
               & M(1:I-1, J))
          M(J, I) = M(J, I) - Dot_Product(Conjg(M(J, 1:I-1)),&
               & M(1:I-1, I))

          If (Abs(M(I,I)) < Epsilon(1.0_SP)) &
               &CALL Abort('lu_SP','Singular matrix.')
          M(J, I) = M(J, I) / M(I, I)
       End Do
    End Do

    Return
  End Subroutine lu_DPC

!  *********************************************
!  *                                           *
  Subroutine lusolve_DPC(M, b)
!  *                                           *
!  *********************************************
!  * Solve a linear set of equations using LU 
!  * decomposition. Both M and b are overwritten
!  *********************************************

    Complex (kind=DPC), Intent (inout) :: M(:,:), b(:)

    Complex (kind=DPC) :: bcp(Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, Idim, I
    

    Idim = Size(M,1)
    CALL LU(M, Ipiv, Id)

    bcp = b
    ! Now we Permutate b
    Do I = 1, Idim
       b(I) = bcp(Ipiv(I))
    End Do

    ! First solve Lx = b
    Do I = 2, Idim
       b(I) = b(I) - Dot_Product(Conjg(M(I,1:I-1)),b(1:I-1))
    End Do
    
    ! Now solve Ux = b
    b(Idim) = b(Idim) / M(Idim, Idim)
    Do I = Idim - 1, 1, -1
       b(I) = b(I) - Dot_Product(Conjg(M(I, I+1:Idim)),b(I+1:Idim))
       b(I) = b(I) / M(I, I)
    End Do

    Return
  End Subroutine lusolve_DPC

!  *********************************************
!  *                                           *
  Complex (kind=DPC) Function Det_DPC(M)
!  *                                           *
!  *********************************************
!  * Compute the determinant of the matrix M.
!  *********************************************

    Complex (kind=DPC), Intent (in) :: M(:,:)
    
    Complex (kind=DPC) :: Mcp(Size(M,1),Size(M,1))
    Integer :: Ipiv(Size(M,1)), Id, I, J, Idim

    Idim = Size(M,1)
    Mcp = M
    ! First make pivoting
    CALL Pivoting(Mcp, Ipiv, Id)

    
    ! NOW: LU Decomposition
    ! Separamos el paso I = 1
    Do J = 2, Idim
       Mcp(J, 1) = Mcp(J, 1) / Mcp(1, 1)
    End Do

    Det_DPC = Mcp(1,1)
    Do I = 2, Idim
       Mcp(I, I) = Mcp(I, I) - &
            & Dot_Product(Conjg(Mcp(I,1:I-1)), Mcp(1:I-1, I)) 

       Do J = i+1, Idim
          Mcp(I, J) = Mcp(I, J) - Dot_Product(Conjg(Mcp(I, 1:I-1)),&
               & Mcp(1:I-1, J))
          Mcp(J, I) = Mcp(J, I) - Dot_Product(Conjg(Mcp(J, 1:I-1)),&
               & Mcp(1:I-1, I))

          If (Abs(Mcp(I,I)) < Epsilon(1.0_SP)) Then
             Det_DPC = 0.0_SP
             Return
          End If
          Mcp(J, I) = Mcp(J, I) / Mcp(I, I)
       End Do
       Det_DPC = Det_DPC * Mcp(I,I)
    End Do

    Det_DPC = Cmplx(Real(Id,kind=SP),0.0_SP)*Det_DPC

    Return
  End Function Det_DPC

!  *********************************************
!  *                                           *
  Function Cholesky_DP(M) Result (L)
!  *                                           *
!  *********************************************
!  * Makes Cholesky decomposition of Matrix M
!  *********************************************

    Real (kind=DP), Intent (inout) :: M(:,:)

    Real (kind=DP) :: L(Size(M,1),Size(M,2)), R
    Integer :: I, J

    L = 0.0_DP
    Do J = 1, Size(M,2)
       R = M(J,J) - Sum(L(J,1:J-1)**2)
       If (R > 0.0_DP) Then
          L(J,J) = Sqrt(R)
       Else
          CALL abort("Cholesky", "Matrix not positive definite")
       End If
       Do I = J+1, Size(M,1)
        L(I,J) = 1.0_DP/L(J,J) * ( M(I,J) - Sum(L(I,1:J-1)*L(J,1:J-1)) )  
       End Do
    End Do

    Return
  End Function Cholesky_DP

!  *********************************************
!  *                                           *
  Function Cholesky_SP(M) Result (L)
!  *                                           *
!  *********************************************
!  * Makes Cholesky decomposition of Matrix M
!  *********************************************

    Real (kind=SP), Intent (inout) :: M(:,:)

    Real (kind=SP) :: L(Size(M,1),Size(M,2)), R
    Integer :: I, J

    L = 0.0_SP
    Do J = 1, Size(M,2)
       R = M(J,J) - Sum(L(J,1:J-1)**2)
       If (R > 0.0_SP) Then
          L(J,J) = Sqrt(R)
       Else
          CALL abort("Cholesky", "Matrix not positive definite")
       End If
       Do I = J+1, Size(M,1)
        L(I,J) = 1.0_SP/L(J,J) * ( M(I,J) - Sum(L(I,1:J-1)*L(J,1:J-1)) )  
       End Do
    End Do

    Return
  End Function Cholesky_SP

!  *********************************************
!  *                                           *
  Function Cholesky_DPC(M) Result (L)
!  *                                           *
!  *********************************************
!  * Makes Cholesky decomposition of Matrix M
!  *********************************************

    Complex (kind=DPC), Intent (inout) :: M(:,:)

    Complex (kind=DPC) :: L(Size(M,1),Size(M,2))
    Real (kind=DP) :: R
    Integer :: I, J

    L = 0.0_DPC
    Do J = 1, Size(M,2)
       R = M(J,J) - Sum(Abs(L(J,1:J-1))**2)
       If (R > 0.0_DPC) Then
          L(J,J) = Sqrt(R)
       Else
          CALL abort("Cholesky", "Matrix not positive definite")
       End If
       Do I = J+1, Size(M,1)
        L(I,J) = 1.0_DPC/L(J,J) * ( M(I,J) - Sum(L(I,1:J-1)*Conjg(L(J,1:J-1))) )  
       End Do
    End Do

    Return
  End Function Cholesky_DPC


!  *********************************************
!  *                                           *
  Function Cholesky_SPC(M) Result (L)
!  *                                           *
!  *********************************************
!  * Makes Cholesky decomposition of Matrix M
!  *********************************************

    Complex (kind=SPC), Intent (inout) :: M(:,:)

    Complex (kind=SPC) :: L(Size(M,1),Size(M,2))
    Real (kind=SP) :: R
    Integer :: I, J

    L = 0.0_SPC
    Do J = 1, Size(M,2)
       R = M(J,J) - Sum(Abs(L(J,1:J-1))**2)
       If (R > 0.0_SPC) Then
          L(J,J) = Sqrt(R)
       Else
          CALL abort("Cholesky", "Matrix not positive definite")
       End If
       Do I = J+1, Size(M,1)
        L(I,J) = 1.0_SPC/L(J,J) * ( M(I,J) - Sum(L(I,1:J-1)*Conjg(L(J,1:J-1))) )  
       End Do
    End Do

    Return
  End Function Cholesky_SPC



End MODULE Linear
