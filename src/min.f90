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
  Real (kind=DP), Parameter :: EPSILON_DP = 1.0E-10_DP
  Real (kind=SP), Parameter :: EPSILON_SP = 1.0E-10_SP
    


  Interface Step
     Module Procedure Step_DP, MultiStep_DP, Step_SP, MultiStep_SP
  End Interface

  Interface MaxPosition
     Module Procedure MaxPosition_2D_DP!, MaxPosition_2D_SP
  End Interface

  Interface InterpolMinMax
     Module Procedure InterpolMinMax_2D_DP
  End Interface

  Interface Bracket
     Module Procedure Bracket_DP, Bracket_SP
  End Interface

  Interface LineSrch
     Module Procedure LineSrch_DP, LineSrch_SP, LineSrchMulti_DP, &
          & LineSrchMulti_SP
  End Interface

  Interface ConjGrad
     Module Procedure ConjGrad_DP, ConjGrad_SP
  End Interface


  Private Step_DP, MultiStep_DP, Step_SP, MultiStep_SP, &
       & MaxPosition_2D_DP, Bracket_DP, Bracket_SP, LineSrch_DP, &
       & LineSrch_SP, ConjGrad_DP, ConjGrad_SP, LineSrchMulti_DP, &
       & LineSrchMulti_SP
       !, MaxPosition_2D_SP, InterpolMinMax_2D_DP

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

    Real (kind=DP) :: Zcp(0:Size(Z,1)+1, 0:Size(Z,2)+1)
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

! **********************************************
! *                                            *
  Subroutine LineSrch_DP(X2, Func, Tolerance) 
! *                                            *
! **********************************************
! * Minimize funcion Func with precision Tol   *
! * (if present), using golden rule line       *
! * minimisation. The minimum should be        *
! * bracketed by [X1,X2,X3]                    * 
! **********************************************

    Real (kind=DP), Intent (inout) :: X2
    Real (kind=DP), Intent (in), Optional :: Tolerance
    
    Real (kind=DP) :: Tol, X4, F1, F2, F3, F4, X1, X3
    
    Interface 
       Function Func(X)
         USE NumTypes
         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (Present(Tolerance)) Then
       Tol = Tolerance
    Else
       Tol = DEFTOL_DP
    End If

    CALL Bracket(X1,X2,X3,Func)
    X4 = X1 - X2 + X3
    F1 = Func(X1)
    F2 = Func(X2)
    F3 = Func(X3)
    F4 = Func(X4)

    Do While (X3 - X1 > Tol*(Abs(X3) + Abs(X1) + EPSILON_DP))
       If (F4 < F2) Then
          If (X2 < X4) Then
             X1 = X2
             F1 = F2
          Else
             X3 = X2
             F3 = F2
          End If
          X2 = X4
          F2 = F4
       Else
          If (X2 < X4) Then
             X3 = X4
             F3 = F4
          Else
             X1 = X4
             F1 = F4
          End If
       End If
       X4 = X1 - X2 + X3
       F4 = Func(X4)
    End Do
    
    If (Func( (X3+X2)/2.0_DP) < F4) Then
       X2 = (X3+X2)/2.0_DP
    Else
       X2 = X4
    End If

    Return
  End Subroutine LineSrch_DP


 ! **********************************************
! *                                            *
  Subroutine Bracket_DP(X1, X2, X3, Func) 
! *                                            *
! **********************************************
! * Brcket a minimum of Func with the golden   *
! * proportion. X2 is used as a guess of the   *
! * region of the minimum can be introduced    *
! * in X2.                                     *
! **********************************************

    Real (kind=DP), Parameter :: GOLD = 1.6180339887498948482_DP

    Real (kind=DP), Intent (inout) :: X1, X2, X3

    Integer :: I, kont, Maxstep = 20
    Real (kind=DP) :: F2


    Interface 
       Function Func(X)
         USE NumTypes
         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    
    kont = 0
    Do While (.True.)
       F2 = Func(X2)
       Do I = 0, Maxstep
          X1 = X2 - 2.0_DP**I
          X3 = X2 + 2.0_DP**I * GOLD
          If ( (Func(X3) > F2).and.(Func(X1) > F2) ) Return
       End Do
       ! If dont suceed in braketing, we have to move
       X2 = X2 + (-2.0_DP)**kont
       kont = Kont + 1
    End Do
    
    Return
  End Subroutine Bracket_DP

! **********************************************
! *                                            *
  Subroutine ConjGrad_DP(X, Func, FuncD, Tolerance) 
! *                                            *
! **********************************************
! * Minimises a function of several variables  *
! * using the conjugate gradient.              *
! **********************************************

    Real (kind=DP), Intent (inout) :: X(:)
    Real (kind=DP), Intent (in), Optional :: Tolerance

    Real (kind=DP) :: B(Size(X), Size(X)), beta(Size(X)), &
         & Fd(Size(X)), Tol, P(Size(X)), D(Size(X)), &
         & G(Size(X), Size(X)), X2, Fold, Fnew
    Integer :: I


    Interface
       Function Func(X)
         USE NumTypes
         Real (kind=DP), Intent (in)  :: X(:)
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    Interface
       Subroutine FuncD(X, Fd)
         USE NumTypes
         Real (kind=DP), Intent (in)  :: X(:)
         Real (kind=DP), Intent (out) :: Fd(Size(X))
       End Subroutine FuncD
    End Interface

    If (Present(Tolerance)) Then
       Tol = Tolerance
    Else
       Tol = DEFTOL_DP
    End If


    Fold = 235623763.0_DP
    Do While (.True.)
       CALL FuncD(X, Fd)
       G(1,:) = - Fd(:)
       If (Sum(G(1,:)**2) < Tol) Return
       B(1,:) = G(1,:)
       
       X2 = 0.0_DP
       P(:) = X(:)
       D(:) = B(1,:)
       CALL LineSrch(X2,Fdir, Tol)
       X(:) = P(:) + X2 * D(:)
       
       Do I = 2, Size(X)
          CALL FuncD(X, Fd)
          G(I,:) = - Fd(:)
          If (Sum(G(I,:)**2) < Tol) Return
          beta(I) = Sum(G(I,:)*(G(I,:)-G(I-1,:)) ) / Sum(G(I-1,:)*G(I-1,:))
          beta(I) = Max(Beta(I), 0.0_DP)
          B(I,:) = G(I,:) + beta(I)*B(I-1,:)
          
          X2 = 0.0_DP
          P(:) = X(:)
          D(:) = B(I,:)
          CALL LineSrch(X2,Fdir, Tol)
          X(:) = P(:) + X2 * D(:)
       End Do
       Fnew = Func(X)
       If (2.0_DP*Abs(Fnew-Fold) < Tol*(Abs(Fnew)+Abs(Fold)+Epsilon_DP))&
          & Return
       Fold = Fnew
    End Do


    Return

    CONTAINS

      ! *************************************
      Real (kind=DP) Function Fdir(h)
      ! *************************************
      ! * Computes the value of Func in the
      ! * direction Dir(:)
      ! *************************************


        Real (kind=DP), Intent (in) :: h
        
        Fdir = Func(P(:) + h*D(:))
        
        Return
      End Function Fdir

  End Subroutine ConjGrad_DP

! **********************************************
! *                                            *
  Subroutine LineSrchMulti_DP(X, Func, Tolerance) 
! *                                            *
! **********************************************
! * Minimize funcion Func with precision Tol   *
! * (if present), using golden rule line       *
! * minimisation. The minimum should be        *
! * bracketed by [X1,X2,X3]                    * 
! **********************************************

    Real (kind=DP), Intent (inout) :: X(:)
    Real (kind=DP), Intent (in), Optional :: Tolerance
    
    Real (kind=DP) :: Tol, X2, Fold, Fnew
    
    Integer :: Idir

    Interface 
       Function Func(X)
         USE NumTypes
         Real (kind=DP), Intent (in) :: X(:)
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (Present(Tolerance)) Then
       Tol = Tolerance
    Else
       Tol = DEFTOL_DP
    End If

    Fnew = 0.0_DP
    Do While (.True.)
       Fold = Fnew
       Do Idir = 1, Size(X)
          CALL LineSrch(X2, Fdir, Tol)
          X(Idir) = X2
       End Do
       Fnew = Func(X)
       If (2.0_SP*Abs(Fnew-Fold) < Tol*(Abs(Fnew)+Abs(Fold)+Epsilon_SP))&
          & Return
    End Do


    Return
    CONTAINS

      ! *************************************
      Real (kind=DP) Function Fdir(h)
      ! *************************************
      ! * Computes the value of Func in the
      ! * direction Dir(:)
      ! *************************************


        Real (kind=DP), Intent (in) :: h
        Real (kind=DP) :: XX(Size(X))

        XX = X
        Xx(Idir) = h
        Fdir = Func(Xx)
        
        Return
      End Function Fdir


  End Subroutine LineSrchMulti_DP

! **********************************************
! *                                            *
  Subroutine LineSrch_SP(X2, Func, Tolerance) 
! *                                            *
! **********************************************
! * Minimize funcion Func with precision Tol   *
! * (if present), using golden rule line       *
! * minimisation. The minimum should be        *
! * bracketed by [X1,X2,X3]                    * 
! **********************************************

    Real (kind=SP), Intent (inout) :: X2
    Real (kind=SP), Intent (in), Optional :: Tolerance
    
    Real (kind=SP) :: Tol, X4, F1, F2, F3, F4, X1, X3
    
    Interface 
       Function Func(X)
         USE NumTypes
         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (Present(Tolerance)) Then
       Tol = Tolerance
    Else
       Tol = DEFTOL_SP
    End If

    CALL Bracket(X1,X2,X3,Func)
    X4 = X1 - X2 + X3
    F1 = Func(X1)
    F2 = Func(X2)
    F3 = Func(X3)
    F4 = Func(X4)

    Do While (X3 - X1 > Tol*(Abs(X3) + Abs(X1) + EPSILON_SP))
       If (F4 < F2) Then
          If (X2 < X4) Then
             X1 = X2
             F1 = F2
          Else
             X3 = X2
             F3 = F2
          End If
          X2 = X4
          F2 = F4
       Else
          If (X2 < X4) Then
             X3 = X4
             F3 = F4
          Else
             X1 = X4
             F1 = F4
          End If
       End If
       X4 = X1 - X2 + X3
       F4 = Func(X4)
    End Do
    
    If (Func( (X3+X2)/2.0_SP) < F4) Then
       X2 = (X3+X2)/2.0_SP
    Else
       X2 = X4
    End If

    Return
  End Subroutine LineSrch_SP


 ! **********************************************
! *                                            *
  Subroutine Bracket_SP(X1, X2, X3, Func) 
! *                                            *
! **********************************************
! * Brcket a minimum of Func with the golden   *
! * proportion. X2 is used as a guess of the   *
! * region of the minimum can be introduced    *
! * in X2.                                     *
! **********************************************

    Real (kind=SP), Parameter :: GOLD = 1.6180339887498948482_SP

    Real (kind=SP), Intent (inout) :: X1, X2, X3

    Integer :: I, kont, Maxstep = 20
    Real (kind=SP) :: F2


    Interface 
       Function Func(X)
         USE NumTypes
         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    
    kont = 0
    Do While (.True.)
       F2 = Func(X2)
       Do I = 0, Maxstep
          X1 = X2 - 2.0_SP**I
          X3 = X2 + 2.0_SP**I * GOLD
          If ( (Func(X3) > F2).and.(Func(X1) > F2) ) Return
       End Do
       ! If dont suceed in braketing, we have to move
       X2 = X2 + (-2.0_SP)**kont
       kont = Kont + 1
    End Do
    
    Return
  End Subroutine Bracket_SP

! **********************************************
! *                                            *
  Subroutine ConjGrad_SP(X, Func, FuncD, Tolerance) 
! *                                            *
! **********************************************
! * Minimises a function of several variables  *
! * using the conjugate gradient.              *
! **********************************************

    Real (kind=SP), Intent (inout) :: X(:)
    Real (kind=SP), Intent (in), Optional :: Tolerance

    Real (kind=SP) :: B(Size(X), Size(X)), beta(Size(X)), &
         & Fd(Size(X)), Tol, P(Size(X)), D(Size(X)), &
         & G(Size(X), Size(X)), X2, Fold, Fnew
    Integer :: I


    Interface
       Function Func(X)
         USE NumTypes
         Real (kind=SP), Intent (in)  :: X(:)
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    Interface
       Subroutine FuncD(X, Fd)
         USE NumTypes
         Real (kind=SP), Intent (in)  :: X(:)
         Real (kind=SP), Intent (out) :: Fd(Size(X))
       End Subroutine FuncD
    End Interface

    If (Present(Tolerance)) Then
       Tol = Tolerance
    Else
       Tol = DEFTOL_SP
    End If


    Fold = 235623763.0_SP
    Do While (.True.)
       CALL FuncD(X, Fd)
       G(1,:) = - Fd(:)
       If (Sum(G(1,:)**2) < Tol) Return
       B(1,:) = G(1,:)
       
       X2 = 0.0_SP
       P(:) = X(:)
       D(:) = B(1,:)
       CALL LineSrch(X2,Fdir, Tol)
       X(:) = P(:) + X2 * D(:)
       
       Do I = 2, Size(X)
          CALL FuncD(X, Fd)
          G(I,:) = - Fd(:)
          If (Sum(G(I,:)**2) < Tol) Return
          beta(I) = Sum(G(I,:)*(G(I,:)-G(I-1,:)) ) / Sum(G(I-1,:)*G(I-1,:))
          beta(I) = Max(Beta(I), 0.0_SP)
          B(I,:) = G(I,:) + beta(I)*B(I-1,:)
          
          X2 = 0.0_SP
          P(:) = X(:)
          D(:) = B(I,:)
          CALL LineSrch(X2,Fdir, Tol)
          X(:) = P(:) + X2 * D(:)
       End Do
       Fnew = Func(X)
       If (2.0_SP*Abs(Fnew-Fold) < Tol*(Abs(Fnew)+Abs(Fold)+Epsilon_SP))&
          & Return
       Fold = Fnew
    End Do


    Return

    CONTAINS

      ! *************************************
      Real (kind=SP) Function Fdir(h)
      ! *************************************
      ! * Computes the value of Func in the
      ! * direction Dir(:)
      ! *************************************


        Real (kind=SP), Intent (in) :: h
        
        Fdir = Func(P(:) + h*D(:))
        
        Return
      End Function Fdir

  End Subroutine ConjGrad_SP

! **********************************************
! *                                            *
  Subroutine LineSrchMulti_SP(X, Func, Tolerance) 
! *                                            *
! **********************************************
! * Minimize funcion Func with precision Tol   *
! * (if present), using golden rule line       *
! * minimisation. The minimum should be        *
! * bracketed by [X1,X2,X3]                    * 
! **********************************************

    Real (kind=SP), Intent (inout) :: X(:)
    Real (kind=SP), Intent (in), Optional :: Tolerance
    
    Real (kind=SP) :: Tol, X2, Fold, Fnew
    
    Integer :: Idir

    Interface 
       Function Func(X)
         USE NumTypes
         Real (kind=SP), Intent (in) :: X(:)
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (Present(Tolerance)) Then
       Tol = Tolerance
    Else
       Tol = DEFTOL_SP
    End If

    Fnew = 0.0_SP
    Do While (.True.)
       Fold = Fnew
       Do Idir = 1, Size(X)
          CALL LineSrch(X2, Fdir, Tol)
          X(Idir) = X2
       End Do
       Fnew = Func(X)
       If (2.0_SP*Abs(Fnew-Fold) < Tol*(Abs(Fnew)+Abs(Fold)+Epsilon_SP))&
          & Return
    End Do


    Return
    CONTAINS

      ! *************************************
      Real (kind=SP) Function Fdir(h)
      ! *************************************
      ! * Computes the value of Func in the
      ! * direction Dir(:)
      ! *************************************


        Real (kind=SP), Intent (in) :: h
        Real (kind=SP) :: XX(Size(X))

        XX = X
        Xx(Idir) = h
        Fdir = Func(Xx)
        
        Return
      End Function Fdir


  End Subroutine LineSrchMulti_SP

 
End MODULE Optimization
