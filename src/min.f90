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

  Private Step_DP, MultiStep_DP, Step_SP, MultiStep_SP

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



End MODULE Optimization
