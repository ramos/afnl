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


  Interface Pivoting
     Module Procedure Pivoting_SP, Pivoting_DP
  End Interface

  Interface lu
     Module Procedure lu_SP, lu_DP
  End Interface

  Interface lusolve
     Module Procedure lusolve_SP, lusolve_DP
  End Interface

  Interface Det
     Module Procedure Det_SP, Det_DP
  End Interface

  
  Private Pivoting_SP, lu_SP, lusolve_SP, Det_SP, &
       & Pivoting_DP, lu_DP, lusolve_DP, Det_DP

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
    Integer :: Ipos(1), Kval

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
    Integer :: Ipiv(Size(M,1)), Id, Idim
    

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
    Integer :: Ipiv(Size(M,1)), Id

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
       Det_SP = Det_SP * M(I,I)
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
    Integer :: Ipos(1), Kval

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
    Integer :: Ipiv(Size(M,1)), Id, Idim
    

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
    Integer :: Ipiv(Size(M,1)), Id

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
       Det_DP = Det_DP * M(I,I)
    End Do

    Det_DP = Real(Id,kind=DP)*Det_DP

    Return
  End Function Det_DP

End MODULE Linear
