!
! MODULE to minimize functions
!
! Copyright (C) 2004  Alberto Ramos <alberto@martin.ft.uam.es>
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

! $ v. 1.0; Released: 10/03/2004; $


MODULE Optimization

  USE NumTypes
  USE Error

  IMPLICIT NONE

  Real (kind=DP), Parameter :: DEFTOL_DP = 1.0E-3_DP
  Real (kind=SP), Parameter :: DEFTOL_SP = 1.0E-3_SP

  Interface Step
     Module Procedure Step_DP, MultiStep_DP, Step_SP, MultiStep_SP
  End Interface

  Interface MaxPosition
     Module Procedure MaxPosition_2D_DP!, MaxPosition_2D_SP
  End Interface

  Interface InterpolMinMax
     Module Procedure InterpolMinMax_2D_DP
  End Interface

  Private Step_DP, MultiStep_DP, Step_SP, MultiStep_SP, &
       & MaxPosition_2D_DP!, MaxPosition_2D_SP, InterpolMinMax_2D_DP

CONTAINS

! **********************************************
! *                                            *
  Real (kind=DP) Function Step_DP(Xo, FStep, Tolerance) 
! *                                            *
! **********************************************
! * Minimiza una funcion, FStep() de una       *
! * variable por el metodo de Stepping.        *
! **********************************************

    Real (kind=DP) :: h, Rr, Tol
    Real (kind=DP), Intent (in) :: Xo
    Real (kind=DP), Intent (in), Optional :: Tolerance
    
    Interface
       Function Fstep(Xo)
         USE NumTypes
         
         Real (kind=DP), Intent (in) :: Xo
         Real (kind=DP) :: Fstep
       End Function Fstep
    End Interface
    
    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_DP
    Else
       Tol = Tolerance
    End If

    h = 0.3_DP
    Step_DP = Xo
    Do While (h.gt.Tol)
       Rr = FStep(Step_DP)
       If (Rr.gt.FStep(Step_DP+h)) Then
          Step_DP = Step_DP + h
       Else If (Rr.gt.FStep(Step_DP-h)) Then
          Step_DP = Step_DP - h
       Else 
          h = h/2.D0
       End If
    End Do
    
    Return
  End Function Step_DP

! **********************************************
! *                                            *
  Function MultiStep_DP(Xo, FStep, Tolerance) Result (Mstep)
! *                                            *
! **********************************************
! * Minimiza una funcion, FStep() de varias    *
! * variables por el metodo de Stepping.       *
! **********************************************

    Real (kind=DP), Intent (in) :: Xo(:)
    Real (kind=DP), Intent (in), Optional :: Tolerance

    Real (kind=DP) :: h, Tol, Rr
    Real (kind=DP) :: Mstep(Size(Xo)), hvec(Size(Xo))
    Logical :: Ok
    Integer :: I
    
    Interface
       Function Fstep(Xo)
         USE NumTypes

         Real (kind=DP), Intent (in) :: Xo(:)
         Real (kind=DP) :: Fstep
       End Function Fstep
    End Interface
    
    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_DP
    Else
       Tol = Tolerance
    End If

    h = 0.5_DP
    Mstep = Xo

    Do While (h.gt.Tol)
       Do I = 1, Size(Xo)
          hvec = 0.0_DP
          hvec(I) = h
          Ok = .True.
          Do While (Ok)
             Rr = FStep(Mstep)
             If (Rr.gt.FStep(Mstep+hvec)) Then
                Mstep(I) = Mstep(I) + h
             Else If (Rr.gt.FStep(Mstep-hvec)) Then
                Mstep(I) = Mstep(I) - h
             Else
                Ok = .False.
             End If
          End Do
       End Do
       h = h / 2.0_DP
    End Do
       

    Return
  End Function MultiStep_DP


! **********************************************
! *                                            *
  Real (kind=SP) Function Step_SP(Xo, FStep, Tolerance) 
! *                                            *
! **********************************************
! * Minimiza una funcion, FStep() de una       *
! * variable por el metodo de Stepping.        *
! **********************************************

    Real (kind=SP) :: h, Rr, Tol
    Real (kind=SP), Intent (in) :: Xo
    Real (kind=SP), Intent (in), Optional :: Tolerance
    
    Interface
       Function Fstep(Xo)
         USE NumTypes
         
         Real (kind=SP), Intent (in) :: Xo
         Real (kind=SP) :: Fstep
       End Function Fstep
    End Interface
    
    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_SP
    Else
       Tol = Tolerance
    End If

    h = 0.3_SP
    Step_SP = Xo
    Do While (h.gt.Tol)
       Rr = FStep(Step_SP)
       If (Rr.gt.FStep(Step_SP+h)) Then
          Step_SP = Step_SP + h
       Else If (Rr.gt.FStep(Step_SP-h)) Then
          Step_SP = Step_SP - h
       Else 
          h = h/2.D0
       End If
    End Do
    
    Return
  End Function Step_SP

! **********************************************
! *                                            *
  Function MultiStep_SP(Xo, FStep, Tolerance) Result (Mstep)
! *                                            *
! **********************************************
! * Minimiza una funcion, FStep() de una       *
! * variable por el metodo de Stepping.        *
! **********************************************

    Real (kind=SP), Intent (in) :: Xo(:)
    Real (kind=SP), Intent (in), Optional :: Tolerance

    Real (kind=SP) :: h, Tol, Rr
    Real (kind=SP) :: Mstep(Size(Xo)), hvec(Size(Xo))
    Logical :: Ok
    Integer :: I
    
    Interface
       Function Fstep(Xo)
         USE NumTypes

         Real (kind=SP), Intent (in) :: Xo(:)
         Real (kind=SP) :: Fstep
       End Function Fstep
    End Interface
    
    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_SP
    Else
       Tol = Tolerance
    End If

    h = 0.5_SP
    Mstep = Xo

    Do While (h.gt.Tol)
       Do I = 1, Size(Xo)
          hvec = 0.0_SP
          hvec(I) = h
          Ok = .True.
          Do While (Ok)
             Rr = FStep(Mstep)
             If (Rr.gt.FStep(Mstep+hvec)) Then
                Mstep(I) = Mstep(I) + h
             Else If (Rr.gt.FStep(Mstep-hvec)) Then
                Mstep(I) = Mstep(I) - h
             Else
                Ok = .False.
             End If
          End Do
       End Do
       h = h / 2.0_SP
    End Do
       

    Return
  End Function MultiStep_SP

! **********************************************
! *                                            *
  Integer Function MaxPosition_2D_DP(Z, IposX, IposY) 
! *                                            *
! **********************************************
! * Locates the number of maximums and its 
! * positions of the data Z(:,:).
! **********************************************

    Real (kind=DP), Intent (in) :: Z(:,:)
    Integer, Intent (out) :: IposX(:), IposY(:)

    Real (kind=DP) :: Val, Zcp(0:Size(Z,1)+1, 0:Size(Z,2)+1)
    Integer :: I, J, Ipos(2)

    MaxPosition_2D_DP = 0
    
    Zcp(1:Size(Z,1),1:Size(Z,2)) = Z
    Zcp(0,1:Size(Z,2)) = Z(Size(Z,1),:)
    Zcp(Size(Z,1)+1,1:Size(Z,2)) = Z(1,:)
    Zcp(1:Size(Z,1),0) = Z(:,Size(Z,2))
    Zcp(1:Size(Z,1),Size(Z,2)+1) = Z(:,1)
    
    Zcp(0,0) = Z(Size(Z,1),Size(Z,2))
    Zcp(0,Size(Z,2)+1) = Z(Size(Z,1),1)
    Zcp(Size(Z,1)+1,0) = Z(1,Size(Z,2))
    Zcp(Size(Z,1)+1,Size(Z,2)+1) = Z(1,1)
    
    Do I = 1, Size(Z,1)
       Do J = 1, Size(Z,2)
          Ipos = MaxLoc(Z(I-1:I+1,J-1:J+1))
          If ( (Ipos(1) == 2).and.(Ipos(2) == 2) ) Then
             MaxPosition_2D_DP = MaxPosition_2D_DP + 1
             IposX(MaxPosition_2D_DP) = I
             IposY(MaxPosition_2D_DP) = J
          End If
       End Do
    End Do

    Return
  End Function MaxPosition_2D_DP

! **********************************************
! *                                            *
  Function InterpolMinMax_2D_DP(X, Y, Z) Result (PosMax) 
! *                                            *
! **********************************************
! * Locates the number of maximums and its 
! * positions of the data Z(:,:).
! **********************************************

    USE Linear
    
    Real (kind=DP), Intent (in) :: X(:), Y(:), Z(:,:)
    Real (kind=DP) :: PosMax(2)

    Integer :: I, J
    Real (kind=DP) :: Ma(5,5), b(5)
    
    
    Ma(1,1) = X(1)**2
    Ma(1,2) = X(1)
    Ma(1,3) = Y(2)**2
    Ma(1,4) = Y(2)
    Ma(1,5) = 1.0_DP
    b(1) = Z(1,2)

    Ma(2,1) = X(2)**2
    Ma(2,2) = X(2)
    Ma(2,3) = Y(3)**2
    Ma(2,4) = Y(3)
    Ma(2,5) = 1.0_DP
    b(2) = Z(2,3)

    Ma(3,1) = X(2)**2
    Ma(3,2) = X(2)
    Ma(3,3) = Y(2)**2
    Ma(3,4) = Y(2)
    Ma(3,5) = 1.0_DP
    b(3) = Z(2,2)

    Ma(4,1) = X(2)**2
    Ma(4,2) = X(2)
    Ma(4,3) = Y(1)**2
    Ma(4,4) = Y(1)
    Ma(4,5) = 1.0_DP
    b(4) = Z(2,1)

    Ma(5,1) = X(3)**2
    Ma(5,2) = X(3)
    Ma(5,3) = Y(2)**2
    Ma(5,4) = Y(2)
    Ma(5,5) = 1.0_DP
    b(5) = Z(3,2)

!    Do I = 1, 5
!       Write(stderr,'(100ES10.2)')(Ma(I,J), J = 1, 5)
!    End Do

    CALL LUSolve(Ma, b)
    PosMax(1) = -b(2)/(2.0_DP*b(1))
    PosMax(2) = -b(4)/(2.0_DP*b(3))
    
    Return
  End Function InterpolMinMax_2D_DP


End MODULE Optimization
