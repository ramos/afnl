
! A MODULE for numerical Integration and to solve ODE's 
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

! $ v. 1.0; Released: 21/09/2003; $

! ***************************************************
! *
MODULE Integration
! *
! ***************************************************

 
  USE NumTypes
  USE Linear

  IMPLICIT NONE

  ! The maximum number of iterations before
  ! Printing error msg and exit.

  Integer, Parameter :: MaxIter = 26

  Real (kind=DP), Parameter :: DEFTOL_DP = 1.E-2_DP
  Real (kind=SP), Parameter :: DEFTOL_SP = 1.E-2_SP

  Real (kind=DP), Parameter :: DEFINTSTEPS_DP = 10.0_DP
  Real (kind=SP), Parameter :: DEFINTSTEPS_SP = 10.0_SP

  ! Here we provide Interfaces, so that SP
  ! functions and DP functions, could be 
  ! called in the same way.

  Interface Trapecio
     Module Procedure Trapecio_DP, Trapecio_SP
  End Interface
  
  Interface Simpson
     Module Procedure Simpson_DP, Simpson_SP
  End Interface

  Interface TrapecioAB
     Module Procedure TrapecioAB_DP, TrapecioAB_SP
  End Interface
  
  Interface SimpsonAB
     Module Procedure SimpsonAB_DP, SimpsonAB_SP
  End Interface

  Interface TrapInfUP
     Module Procedure TrapInfUP_DP, TrapInfUP_SP
  End Interface
  
  Interface SimpsonInfUP
     Module Procedure SimpsonInfUP_DP, SimpsonInfUP_SP
  End Interface

  Interface TrapInfDW
     Module Procedure TrapInfDW_DP, TrapInfDW_SP
  End Interface
  
  Interface SimpsonInfDW
     Module Procedure SimpsonInfDW_DP, SimpsonInfDW_SP
  End Interface

  Interface TrapSingUP
     Module Procedure TrapSingUP_DP, TrapSingUP_SP
  End Interface
  
  Interface SimpsonSingUP
     Module Procedure SimpsonSingUP_DP, SimpsonSingUP_SP
  End Interface

  Interface TrapSingDW
     Module Procedure TrapSingDW_DP, TrapSingDW_SP
  End Interface
  
  Interface Trap
     Module Procedure Trap_DP, Trap_SP
  End Interface
  
  Interface TrapAB
     Module Procedure TrapAB_DP, TrapAB_SP
  End Interface
  
  Interface SimpsonSingDW
     Module Procedure SimpsonSingDW_DP, SimpsonSingDW_SP
  End Interface

  Interface Euler
     Module Procedure Euler_DP, Euler_SP
  End Interface
  
  Interface Rgnkta
     Module Procedure Rgnkta_DP, Rgnkta_SP
  End Interface

  Interface IntQuadrilateral
     Module Procedure IntQuadrilateral_DP, IntQuadrilateral_SP
  End Interface

  ! Trap-like routines are "low level", that must
  ! be used only by functions of this module.
  ! There is no reason why someone must call 
  ! Simpson_DP...

  Private Trap_DP, TrapAb_DP, TrapInfUp_DP, TrapInfDw_DP, &
       & TrapSingUp_DP, TrapSingDw_DP, Trapecio_DP, Trapecio_SP, &
       & Simpson_DP, Simpson_SP, TrapecioAB_DP, TrapecioAB_SP, &
       & SimpsonAB_DP, SimpsonAB_SP, SimpsonInfUP_DP, SimpsonInfUP_SP,&
       & SimpsonInfDW_DP, SimpsonInfDW_SP, SimpsonSingUP_DP, &
       & SimpsonSingUP_SP, SimpsonSingDW_DP, SimpsonSingDW_SP, &
       & Euler_DP, Euler_SP, Rgnkta_DP, Rgnkta_SP, &
       & IntQuadrilateral_DP, IntQuadrilateral_SP

  Private Maxiter, DEFTOL_DP, DEFTOL_SP, DEFINTSTEPS_DP, DEFINTSTEPS_SP

CONTAINS


!  *********************************************
!  *                                           *
  Subroutine Trap_DP(a, b, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. DP 
!  * version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: a, b
    Integer, Intent (in) :: N
    Real (kind=DP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=DP) :: Val, h, X

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface
    

    If (N == 1) Then
       Sum = (b - a) * ( Func(a) + Func(b) ) / 2.0_DP
    Else 
       Ntrap = 2**(N-2)
       Val = 0.0_DP
       h = (b - a) / Real(Ntrap, kind=DP)

       X = a + h / 2.0_DP
       Do I = 1, Ntrap
          Val = Val + Func(X)
          X = X + h
       End Do

       Sum = ( Sum + (b-a)*Val/Real(Ntrap, DP) ) / 2.0_DP
       
    End If
    
    Return
  End Subroutine Trap_DP
    

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Trapecio_DP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the trapezoid rule. 
!  * DP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: a, b
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, Old, ainf, bsup
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)

    CALL Trap(ainf, bsup, Func, Trapecio_DP, 1)
    Old = Trapecio_DP
    Do I = 2, Maxiter
       CALL Trap(ainf, bsup, Func, Trapecio_DP, I)

       If (Abs(Old - Trapecio_DP) < Tol) Then
          If (b < a) Trapecio_DP = -Trapecio_DP 
          Return
       End If
       
       Old = Trapecio_DP
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(0,*)'Trapecio: Problem with convergence'
    Stop

    Return
  End Function Trapecio_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Simpson_DP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the Simpson rule. 
!  * DP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: a, b
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, OldN, Old2N, ainf, bsup, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)
  
    CALL Trap(ainf, bsup, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter
       CALL Trap(ainf, bsup, Func, Old2N, I)
       Simpson_DP = (4.0_DP * Old2N - OldN) / 3.0_DP

       If (Abs(Simpson_DP - OldSum) < Tol) Then
          If (b < a) Simpson_DP = - Simpson_DP
          Return
       End If
       
       OldSum = Simpson_DP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(0,*)'Simpson: Problem with convergence.'
    Stop
    
    Return
  End Function Simpson_DP

!  *********************************************
!  *                                           *
  Subroutine TrapAb_DP(a, b, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * version is a open traps rule (the function 
!  * is never evaluated in a or b). DP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=DP), intent (in) :: a, b
    Integer, Intent (in) :: N
    Real (kind=DP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=DP) :: Val, h, X

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface
    
    
    If (N == 1) Then
       Sum = (b-a) * Func((a+b)/2.0_DP)
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_DP
       h = (b - a) / ( 3.0_DP*Real(Ntrap, DP) )

       X = a + h / 2.0_DP
       Do I = 1, Ntrap
          Val = Val + Func(X)
          X = X + 2.0_DP*h
          Val = Val + Func(X)
          X = X + h
       End Do

       Sum = ( Sum + (b-a)*Val/Real(Ntrap, DP) ) / 3.0_DP
    End If
    

    Return
  End Subroutine TrapAb_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function TrapecioAb_DP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the open trapezoid 
!  * rule. DP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: a, b
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, Old, ainf, bsup
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)

    CALL TrapAb(ainf, bsup, Func, TrapecioAb_DP, 1)
    Old = TrapecioAb_DP
    Do I = 2, Maxiter
       CALL TrapAb(ainf, bsup, Func, TrapecioAb_DP, I)

       If (Abs(Old - TrapecioAb_DP) < Tol) Then
          If (b < a) TrapecioAb_DP = - TrapecioAb_DP
          Return
       End If
       
       Old = TrapecioAb_DP
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'TrapecioAb: Problem with convergence.'
    Stop

    Return
  End Function TrapecioAb_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function SimpsonAb_DP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the open Simpson rule. 
!  * DP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: a, b
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, OldN, Old2N, ainf, bsup, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)
  
    CALL TrapAb(ainf, bsup, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter
       CALL TrapAb(ainf, bsup, Func, Old2N, I)
       SimpsonAb_DP = (9.0_DP * Old2N - OldN) / 8.0_DP

       If (Abs(SimpsonAb_DP - OldSum) < Tol) Then
          If (b < a) SimpsonAb_DP = - SimpsonAb_DP
          Return
       End If
       
       OldSum = SimpsonAb_DP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(0,*)'SimpsonAb: Problem with convergence'
    Stop
    
    Return
  End Function SimpsonAb_DP

!  *********************************************
!  *                                           *
  Subroutine TrapInfUp_DP(aa, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to claculates ingrals with infinite 
!  * upper limit. DP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=DP), intent (in) :: aa
    Integer, Intent (in) :: N
    Real (kind=DP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=DP) :: Val, h, X, b

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface
    
    b = 1.0_DP / aa

    If (N == 1) Then
       Sum = 4.0_DP * Func(2.0_DP / b) / b**2
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_DP
       h = b / ( 3.0_DP*Real(Ntrap, DP) )

       X = h / 2.0_DP
       Do I = 1, Ntrap
          Val = Val + Func(1.0_DP / X) / X**2
          X = X + 2.0_DP*h
          Val = Val + Func(1.0_DP / X) / X**2
          X = X + h
       End Do

       Sum = ( Sum + (b-aa)*Val/Real(Ntrap, DP) ) / 3.0_DP
    End If
    

    Return
  End Subroutine TrapInfUp_DP
  
!  *********************************************
!  *                                           *
  Real (kind=DP) Function SimpsonInfUP_DP(aa, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of Fint between 
!  * aa and infinite. DP version.
!  *
!  * $ Last Mod: 21/09/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: aa
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, OldN, Old2N, Val, OldSum, a
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    If (aa < 1.0_DP) Then
       Val = Simpson(a, 1.0_DP, Func, Eps)
       a = 1.0_DP
    Else
       a = aa
       Val = 0.0_DP
    End If
    

    CALL TrapInfUP(a, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter
       CALL TrapInfUP(a, Func, Old2N, I)
       SimpsonInfUP_DP = (9.0_DP * Old2N - OldN) / 8.0_DP

       If (Abs(SimpsonInfUP_DP - OldSum) < Tol) Then
          SimpsonInfUP_DP = SimpsonInfUP_DP + Val
          Return
       End If
       
       OldSum = SimpsonInfUP_DP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonInfUP: Problem with convergence.'
    Stop

    
  End Function SimpsonInfUP_DP

!  *********************************************
!  *                                           *
  Subroutine TrapInfDw_DP(aa, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to claculates ingrals with infinite 
!  * lower limit. DP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: aa
    Integer, Intent (in) :: N
    Real (kind=DP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=DP) :: Val, h, X, a

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    
    a = 1.0_DP / aa
    
    If (N == 1) Then
       Sum = - 4.0_DP * Func( 2.0_DP/a ) / a**2
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_DP
       h = - a / ( 3.0_DP*Real(Ntrap, DP) )

       X = a + h / 2.0_DP
       Do I = 1, Ntrap
          Val = Val + Func(1.0_DP / X) / X**2
          X = X + 2.0_DP*h
          Val = Val + Func(1.0_DP / X) / X**2
          X = X + h
       End Do

       Sum = ( Sum - a*Val/Real(Ntrap, DP) ) / 3.0_DP
    End If
    

    Return
  End Subroutine TrapInfDw_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function SimpsonInfDw_DP(aa, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of Fint between 
!  * -infinite and aa. DP version.
!  *
!  * $ Last Mod: 21/09/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: aa
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, OldN, Old2N, Val, OldSum, a
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    If (aa > - 1.0_DP) Then
       Val = Simpson(-1.0_DP, aa, Func, Eps)
       a = -1.0_DP
    Else
       a = aa
       Val = 0.0_DP
    End If
    

    Old2N = 0.0_DP
    CALL TrapInfDw(a, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter
       CALL TrapInfDw(a, Func, Old2N, I)
       SimpsonInfDw_DP = (9.0_DP * Old2N - OldN) / 8.0_DP

       If (Abs(SimpsonInfDw_DP - OldSum) < Tol) Then
          SimpsonInfDw_DP = SimpsonInfDw_DP + Val
          Return
       End If
       
       OldSum = SimpsonInfDw_DP
       OldN = Old2N

    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonInfDw: Problem with convergence.'
    Stop

    
  End Function SimpsonInfDw_DP

!  *********************************************
!  *                                           *
  Subroutine TrapSingDw_DP(aa, bb, Func, Sum, N, gamma)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to claculate integrals with a 
!  * singular lower limit. DP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=DP), intent (in) :: aa, bb, gamma
    Integer, Intent (in) :: N
    Real (kind=DP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=DP) :: Val, h, X, b, Exp

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface
    
    b = (bb - aa)**(1.0_DP-gamma)
    Exp = 1.0_DP / (1.0_DP - gamma)
    
    If (N == 1) Then
       Sum = Exp * b * Func( (b/2.0_DP)**(gamma*Exp) + aa) * &
            & (b/2.0_DP)**(gamma*Exp)
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_DP
       h = b / ( 3.0_DP*Real(Ntrap, DP) )

       X = h / 2.0_DP
       Do I = 1, Ntrap
          Val = Val + Func(X**Exp + aa) * X**(gamma*Exp)
          X = X + 2.0_DP*h
          Val = Val + Func(X**Exp + aa) * X**(gamma*Exp)
          X = X + h
       End Do

       Sum = ( Sum + Exp * b * Val/Real(Ntrap, DP) ) / 3.0_DP
    End If

    Return
  End Subroutine TrapSingDw_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function SimpsonSingDw_DP(a, b, Func, Eps, gamma)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the Simpson rule. 
!  * Gamma is the singularity in the lower limit
!  * of Func. DP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: a, b, gamma
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, OldN, Old2N, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    If (gamma .ge. 1.0_DP) Then
       Write(*,*)'SimpsonSingDw: gamma > 1.0_DP, esto NADIE sabe hacerlo.'
       Stop
    End If
    
    CALL TrapSingDw(a, b, Func, Old2N, 1, gamma)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter

       CALL TrapSingDw(a, b, Func, Old2N, I, gamma)
       SimpsonSingDw_DP = (9.0_DP * Old2N - OldN) / 8.0_DP

       If (Abs(SimpsonSingDw_DP - OldSum) < Tol) Then
          SimpsonSingDw_DP = SimpsonSingDw_DP 
          Return
       End If
       
       OldSum = SimpsonSingDw_DP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonSingDw: No hemos convergido.'
    Stop

    
  End Function SimpsonSingDw_DP

!  *********************************************
!  *                                           *
  Subroutine TrapSingUp_DP(aa, bb, Func, Sum, N, gamma)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to claculate integrals with a 
!  * singular upper limit. DP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=DP), intent (in) :: aa, bb, gamma
    Integer, Intent (in) :: N
    Real (kind=DP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=DP) :: Val, h, X, b, Exp

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface
    
    b = (bb - aa)**(1.0_DP-gamma)
    Exp = 1.0_DP / (1.0_DP - gamma)
    
    If (N == 1) Then
       Sum = b * Func( bb - (b/2.0_DP)**(gamma*Exp)) * &
            & (b/2.0_DP)**(gamma*Exp)
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_DP
       h = b / ( 3.0_DP*Real(Ntrap, DP) )

       X = h / 2.0_DP
       Do I = 1, Ntrap
          Val = Val + Func(bb - X**Exp) * X**(gamma*Exp)
          X = X + 2.0_DP*h
          Val = Val + Func(bb - X**Exp) * X**(gamma*Exp)
          X = X + h
       End Do

       Sum = ( Sum + Exp * b * Val/Real(Ntrap, DP) ) / 3.0_DP
    End If

    Return
  End Subroutine TrapSingUp_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function SimpsonSingUp_DP(a, b, Func, Eps, gamma)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the Simpson rule. 
!  * Gamma is the singularity in the upper limit
!  * of Func. DP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=DP), intent (in) :: a, b, gamma
    Real (kind=DP), Intent (in), Optional :: Eps

    Real (kind=DP) :: Tol, OldN, Old2N, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_DP
    Else 
       Tol = Eps
    End If

    If (gamma .ge. 1.0_DP) Then
       Write(*,*)'SimpsonSingUp_DP: gamma > 1.0_DP, esto NADIE sabe hacerlo.'
       Stop
    End If
    
    CALL TrapSingUp(a, b, Func, Old2N, 1, gamma)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter

       CALL TrapSingUp(a, b, Func, Old2N, I, gamma)
       SimpsonSingUp_DP = (9.0_DP * Old2N - OldN) / 8.0_DP

       If (Abs(SimpsonSingUp_DP - OldSum) < Tol) Then
          SimpsonSingUp_DP = SimpsonSingUp_DP 
          Return
       End If
       
       OldSum = SimpsonSingUp_DP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonSingUp_DP: No hemos convergido.'
    Stop

    
  End Function SimpsonSingUp_DP

!  *********************************************
!  *                                           *
  Function Euler_DP(Iniciales,Xo,Xfin,Feuler,Tol) Result (Finales)
!  *                                           *
!  *********************************************
!  * Integrates numerically a set of ODE's 
!  * defined by the function Feuler by the Euler
!  * method. The initial conditions are stored
!  * in iniciales, and integrates between Xo and 
!  * Xfin with tolerance Tol.
!  * 
!  * To integrate a set of ODE's, it is better 
!  * to use Rgnkta. DP version.
!  * 
!  * $ Last Mod: 21/09/2003 $
!  *********************************************


    Real (kind=DP), Intent (in) :: Iniciales(:), Xo, Xfin
    Real (kind=DP), Intent (in), Optional :: Tol
    Real (kind=DP) :: Finales(Size(Iniciales))
    Real (kind=DP) :: R, h

    Interface
       Function Feuler(X, Y) Result (Func)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X, Y(:)
         Real (kind=DP) :: Func(Size(Y))
       End Function Feuler
    End Interface

    If (Present(Tol)) Then
       h = Tol
    Else
       h = (Xfin - Xo)/DEFINTSTEPS_DP
    End If

    R = Xo
    Finales = Iniciales
    Do While (R .lt. Xfin)
       Finales = Finales + h * Feuler(R, Finales)
       R = R + h
    End Do

    Return
  End Function Euler_DP

!  *********************************************
!  *                                           *
  Function Rgnkta_DP(Iniciales,Xo,Xfin,Frgnkta,Tol) Result (Finales)
!  *                                           *
!  *********************************************
!  * Integrates numerically a set of ODE's 
!  * defined by the function Frgnkta by the 
!  * Runge-Kutta method. The initial conditions 
!  * are stored in iniciales, and integrates 
!  * between Xo and Xfin with tolerance Tol.
!  * 
!  * $ Last Mod: 21/09/2003 $
!  *********************************************


    Real (kind=DP), Intent (in) :: Iniciales(:), Xo, Xfin
    Real (kind=DP), Intent (in), Optional :: Tol
    Real (kind=DP) :: Finales(Size(Iniciales))
    Real (kind=DP) :: R, h
    Real (kind=DP) :: R1(Size(Iniciales)), R2(Size(Iniciales)), &
         & R3(Size(Iniciales)), R4(Size(Iniciales)), &
         & Y1(Size(Iniciales)), Aux(Size(Iniciales))
    
    Interface
       Function FRgnkta(X, Y) Result (Func)
         USE NumTypes

         Real (kind=DP), Intent (in) :: X, Y(:)
         Real (kind=DP) :: Func(Size(Y))
       End Function FRgnkta
    End Interface
    
    If (Present(Tol)) Then
       h = Tol
    Else
       h = (Xfin - Xo)/DEFINTSTEPS_DP
    End If

    Y1 = Iniciales
    
    R = Xo
!    Do While (R.lt.Xfin - h)
    Do While (R.lt.Xfin)
       R1 = FRgnkta(R,Y1)
       Aux = Y1 + h*R1/3.0_DP
                
       R2 = FRgnkta(R + h/3.0_DP,Aux)
       Aux = Y1 + h*( -R1/3.0_DP + R2 )
       
       R3 = FRgnkta(R + 2.0_DP*h/3.0_DP,Aux)
       Aux = Y1 + h*( R1 - R2 + R3 )
       
       R4 = FRgnkta(R + h,Aux)
       Y1 = Y1 + (h/8.0_DP)*(R1 + 3.0_DP*R2 + 3.0_DP*R3 + R4)

       R = R + h
    End Do
    
    Finales = Y1
        
    Return
  End Function Rgnkta_DP


!  *********************************************
!  *                                           *
  Subroutine Trap_SP(a, b, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. SP 
!  * version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: a, b
    Integer, Intent (in) :: N
    Real (kind=SP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=SP) :: Val, h, X

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface
    

    If (N == 1) Then
       Sum = (b - a) * ( Func(a) + Func(b) ) / 2.0_SP
    Else 
       Ntrap = 2**(N-2)
       Val = 0.0_SP
       h = (b - a) / Real(Ntrap, SP)

       X = a + h / 2.0_SP
       Do I = 1, Ntrap
          Val = Val + Func(X)
          X = X + h
       End Do

       Sum = ( Sum + (b-a)*Val/Real(Ntrap, SP) ) / 2.0_SP
       
    End If
    
    Return
  End Subroutine Trap_SP
    

!  *********************************************
!  *                                           *
  Real (kind=SP) Function Trapecio_SP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the trapezoid rule. 
!  * SP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: a, b
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, Old, ainf, bsup
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)

    CALL Trap(ainf, bsup, Func, Trapecio_SP, 1)
    Old = Trapecio_SP
    Do I = 2, Maxiter
       CALL Trap(ainf, bsup, Func, Trapecio_SP, I)

       If (Abs(Old - Trapecio_SP) < Tol) Then
          If (b < a) Trapecio_SP = - Trapecio_SP
          Return
       End If
       
       Old = Trapecio_SP
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(0,*)'Trapecio: Problem with convergence'
    Stop

    Return
  End Function Trapecio_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function Simpson_SP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the Simpson rule. 
!  * SP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: a, b
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, OldN, Old2N, ainf, bsup, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)
  
    CALL Trap(ainf, bsup, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter
       CALL Trap(ainf, bsup, Func, Old2N, I)
       Simpson_SP = (4.0_SP * Old2N - OldN) / 3.0_SP

       If (Abs(Simpson_SP - OldSum) < Tol) Then
          If (b < a) Simpson_SP = - Simpson_SP
          Return
       End If
       
       OldSum = Simpson_SP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(0,*)'Simpson: Problem with convergence.'
    Stop
    
    Return
  End Function Simpson_SP

!  *********************************************
!  *                                           *
  Subroutine TrapAb_SP(a, b, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * version is a open traps rule (the function 
!  * is never evaluated in a or b). SP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=SP), intent (in) :: a, b
    Integer, Intent (in) :: N
    Real (kind=SP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=SP) :: Val, h, X

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface
    
    
    If (N == 1) Then
       Sum = (b-a) * Func((a+b)/2.0_SP)
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_SP
       h = (b - a) / ( 3.0_SP*Real(Ntrap, SP) )

       X = a + h / 2.0_SP
       Do I = 1, Ntrap
          Val = Val + Func(X)
          X = X + 2.0_SP*h
          Val = Val + Func(X)
          X = X + h
       End Do

       Sum = ( Sum + (b-a)*Val/Real(Ntrap, SP) ) / 3.0_SP
    End If
    

    Return
  End Subroutine TrapAb_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function TrapecioAb_SP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the open trapezoid 
!  * rule. SP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: a, b
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, Old, ainf, bsup
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)

    CALL TrapAb(ainf, bsup, Func, TrapecioAb_SP, 1)
    Old = TrapecioAb_SP
    Do I = 2, Maxiter
       CALL TrapAb(ainf, bsup, Func, TrapecioAb_SP, I)

       If (Abs(Old - TrapecioAb_SP) < Tol) Then
          If (b < a) TrapecioAb_SP = - TrapecioAb_SP
          Return
       End If
       
       Old = TrapecioAb_SP
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'TrapecioAb: Problem with convergence.'
    Stop

    Return
  End Function TrapecioAb_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function SimpsonAb_SP(a, b, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the open Simpson rule. 
!  * SP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: a, b
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, OldN, Old2N, ainf, bsup, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    ainf = Min(a,b)
    bsup = Max(a,b)
  
    CALL TrapAb(ainf, bsup, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter
       CALL TrapAb(ainf, bsup, Func, Old2N, I)
       SimpsonAb_SP = (9.0_SP * Old2N - OldN) / 8.0_SP

       If (Abs(SimpsonAb_SP - OldSum) < Tol) Then
          If (b < a) SimpsonAb_SP = - SimpsonAb_SP
          Return
       End If
       
       OldSum = SimpsonAb_SP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(0,*)'SimpsonAb: Problem with convergence'
    Stop
    
    Return
  End Function SimpsonAb_SP

!  *********************************************
!  *                                           *
  Subroutine TrapInfUp_SP(aa, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to claculates ingrals with infinite 
!  * upper limit. SP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=SP), intent (in) :: aa
    Integer, Intent (in) :: N
    Real (kind=SP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=SP) :: Val, h, X, a, b

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface
    
    b = 1.0_SP / aa

    If (N == 1) Then
       Sum = 4.0_SP * Func(2.0_SP / b) / b**2
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_SP
       h = b / ( 3.0_SP*Real(Ntrap, SP) )

       X = h / 2.0_SP
       Do I = 1, Ntrap
          Val = Val + Func(1.0_SP / X) / X**2
          X = X + 2.0_SP*h
          Val = Val + Func(1.0_SP / X) / X**2
          X = X + h
       End Do

       Sum = ( Sum + (b-a)*Val/Real(Ntrap, SP) ) / 3.0_SP
    End If
    

    Return
  End Subroutine TrapInfUp_SP
  
!  *********************************************
!  *                                           *
  Real (kind=SP) Function SimpsonInfUP_SP(aa, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of Fint between 
!  * aa and infinite. SP version.
!  *
!  * $ Last Mod: 21/09/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: aa
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, OldN, Old2N, Val, OldSum, a
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    If (aa < 1.0_SP) Then
       Val = Simpson(a, 1.0_SP, Func, Eps)
       a = 1.0_SP
    Else
       a = aa
       Val = 0.0_SP
    End If
    

    CALL TrapInfUP(a, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter
       CALL TrapInfUP(a, Func, Old2N, I)
       SimpsonInfUP_SP = (9.0_SP * Old2N - OldN) / 8.0_SP

       If (Abs(SimpsonInfUP_SP - OldSum) < Tol) Then
          SimpsonInfUP_SP = SimpsonInfUP_SP + Val
          Return
       End If
       
       OldSum = SimpsonInfUP_SP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonInfUP: Problem with convergence.'
    Stop

    
  End Function SimpsonInfUP_SP

!  *********************************************
!  *                                           *
  Subroutine TrapInfDw_SP(aa, Func, Sum, N)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to claculates ingrals with infinite 
!  * lower limit. SP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: aa
    Integer, Intent (in) :: N
    Real (kind=SP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=SP) :: Val, h, X, a

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    
    a = 1.0_SP / aa
    
    If (N == 1) Then
       Sum = - 4.0_SP * Func( 2.0_SP/a ) / a**2
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_SP
       h = - a / ( 3.0_SP*Real(Ntrap, SP) )

       X = a + h / 2.0_SP
       Do I = 1, Ntrap
          Val = Val + Func(1.0_SP / X) / X**2
          X = X + 2.0_SP*h
          Val = Val + Func(1.0_SP / X) / X**2
          X = X + h
       End Do

       Sum = ( Sum - a*Val/Real(Ntrap, SP) ) / 3.0_SP
    End If
    

    Return
  End Subroutine TrapInfDw_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function SimpsonInfDw_SP(aa, Func, Eps)
!  *                                           *
!  *********************************************
!  * Calculates the integral of Fint between 
!  * -infinite and aa. SP version.
!  *
!  * $ Last Mod: 21/09/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: aa
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, OldN, Old2N, Val, OldSum, a
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    If (aa > - 1.0_SP) Then
       Val = Simpson(-1.0_SP, aa, Func, Eps)
       a = -1.0_SP
    Else
       a = aa
       Val = 0.0_SP
    End If
    

    CALL TrapInfDw(a, Func, Old2N, 1)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter

       CALL TrapInfDw(a, Func, Old2N, I)
       SimpsonInfDw_SP = (9.0_SP * Old2N - OldN) / 8.0_SP

       If (Abs(SimpsonInfDw_SP - OldSum) < Tol) Then
          SimpsonInfDw_SP = SimpsonInfDw_SP + Val
          Return
       End If
       
       OldSum = SimpsonInfDw_SP
       OldN = Old2N

    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonInfDw: Problem with convergence.'
    Stop

    
  End Function SimpsonInfDw_SP

!  *********************************************
!  *                                           *
  Subroutine TrapSingDw_SP(aa, bb, Func, Sum, N, gamma)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to calculate integrals with a 
!  * singular lower limit. SP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=SP), intent (in) :: aa, bb, gamma
    Integer, Intent (in) :: N
    Real (kind=SP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=SP) :: Val, h, X, b, Exp

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface
    
    b = (bb - aa)**(1.0_SP-gamma)
    Exp = 1.0_SP / (1.0_SP - gamma)
    
    If (N == 1) Then
       Sum = Exp * b * Func( (b/2.0_SP)**(gamma*Exp) + aa) * &
            & (b/2.0_SP)**(gamma*Exp)
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_SP
       h = b / ( 3.0_SP*Real(Ntrap, SP) )

       X = h / 2.0_SP
       Do I = 1, Ntrap
          Val = Val + Func(X**Exp + aa) * X**(gamma*Exp)
          X = X + 2.0_SP*h
          Val = Val + Func(X**Exp + aa) * X**(gamma*Exp)
          X = X + h
       End Do

       Sum = ( Sum + Exp * b * Val/Real(Ntrap, SP) ) / 3.0_SP
    End If

    Return
  End Subroutine TrapSingDw_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function SimpsonSingDw_SP(a, b, Func, Eps, gamma)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the Simpson rule. 
!  * Gamma is the singularity in the lower limit
!  * of Func. SP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: a, b, gamma
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, OldN, Old2N, Val, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    If (gamma .ge. 1.0_SP) Then
       Write(*,*)'SimpsonSingDw: gamma > 1.0_SP, esto NADIE sabe hacerlo.'
       Stop
    End If
    
    CALL TrapSingDw(a, b, Func, Old2N, 1, gamma)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter

       CALL TrapSingDw(a, b, Func, Old2N, I, gamma)
       SimpsonSingDw_SP = (9.0_SP * Old2N - OldN) / 8.0_SP

       If (Abs(SimpsonSingDw_SP - OldSum) < Tol) Then
          SimpsonSingDw_SP = SimpsonSingDw_SP + Val
          Return
       End If
       
       OldSum = SimpsonSingDw_SP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonSingDw: No hemos convergido.'
    Stop

    
  End Function SimpsonSingDw_SP

!  *********************************************
!  *                                           *
  Subroutine TrapSingUp_SP(aa, bb, Func, Sum, N, gamma)
!  *                                           *
!  *********************************************
!  * If Sum is an estimation of the value of the
!  * integral using 2**N trap's, a call to this
!  * routine will refine this approximation with
!  * an estimation using 2**(N+1) trap's. This 
!  * routine performs a change of variables to
!  * be able to claculate integrals with a 
!  * singular upper limit. SP version.
!  *
!  * This is an "internal routine".
!  * 
!  * $ Last mod: 21/09/2003 $
!  *********************************************
  
    Real (kind=SP), intent (in) :: aa, bb, gamma
    Integer, Intent (in) :: N
    Real (kind=SP), Intent (inout) :: Sum

    Integer :: Ntrap, I
    Real (kind=SP) :: Val, h, X, b, Exp

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface
    
    b = (bb - aa)**(1.0_SP-gamma)
    Exp = 1.0_SP / (1.0_SP - gamma)
    
    If (N == 1) Then
       Sum = b * Func( bb - (b/2.0_SP)**(gamma*Exp)) * &
            & (b/2.0_SP)**(gamma*Exp)
    Else 
       Ntrap = 3**(N-2)
       Val = 0.0_SP
       h = b / ( 3.0_SP*Real(Ntrap, SP) )

       X = h / 2.0_SP
       Do I = 1, Ntrap
          Val = Val + Func(bb - X**Exp) * X**(gamma*Exp)
          X = X + 2.0_SP*h
          Val = Val + Func(bb - X**Exp) * X**(gamma*Exp)
          X = X + h
       End Do

       Sum = ( Sum + Exp * b * Val/Real(Ntrap, SP) ) / 3.0_SP
    End If

    Return
  End Subroutine TrapSingUp_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function SimpsonSingUp_SP(a, b, Func, Eps, gamma)
!  *                                           *
!  *********************************************
!  * Calculates the integral of the function
!  * Func with limits (a,b), and a precision 
!  * Eps (optional) with the Simpson rule. 
!  * Gamma is the singularity in the upper limit
!  * of Func. SP Version.
!  *
!  * $ Last Mod: 21/07/2003 $
!  *********************************************

    Real (kind=SP), intent (in) :: a, b, gamma
    Real (kind=SP), Intent (in), Optional :: Eps

    Real (kind=SP) :: Tol, OldN, Old2N, Val, OldSum
    Integer :: I

    Interface 
       Function Func(X)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: Func
       End Function Func
    End Interface

    If (.not. PRESENT(Eps)) Then
       Tol = DEFTOL_SP
    Else 
       Tol = Eps
    End If

    If (gamma .ge. 1.0_SP) Then
       Write(*,*)'SimpsonSingUp_SP: gamma > 1.0_SP, esto NADIE sabe hacerlo.'
       Stop
    End If
    
    CALL TrapSingUp(a, b, Func, Old2N, 1, gamma)

    OldSum = Old2N + 1.D39*Old2N
    Do I = 2, Maxiter

       CALL TrapSingUp(a, b, Func, Old2N, I, gamma)
       SimpsonSingUp_SP = (9.0_SP * Old2N - OldN) / 8.0_SP

       If (Abs(SimpsonSingUp_SP - OldSum) < Tol) Then
          SimpsonSingUp_SP = SimpsonSingUp_SP + Val
          Return
       End If
       
       OldSum = SimpsonSingUp_SP
       OldN = Old2N
    End Do
    
    ! If we reach this point, we have a problem with
    ! the convergence.

    Write(*,*)'SimpsonSingUp_SP: No hemos convergido.'
    Stop

    
  End Function SimpsonSingUp_SP

!  *********************************************
!  *                                           *
  Function Euler_SP(Iniciales,Xo,Xfin,Feuler,Tol) Result (Finales)
!  *                                           *
!  *********************************************
!  * Integrates numerically a set of ODE's 
!  * defined by the function Feuler by the Euler
!  * method. The initial conditions are stored
!  * in iniciales, and integrates between Xo and 
!  * Xfin with tolerance Tol.
!  * 
!  * To integrate a set of ODE's, it is better 
!  * to use Rgnkta. SP version.
!  * 
!  * $ Last Mod: 21/09/2003 $
!  *********************************************

    Real (kind=SP), Intent (in) :: Iniciales(:), Xo, Xfin
    Real (kind=SP), Intent (in), Optional :: Tol
    Real (kind=SP) :: Finales(Size(Iniciales))
    Real (kind=SP) :: R, h

    Interface
       Function Feuler(X, Y) Result (Func)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X, Y(:)
         Real (kind=SP) :: Func(Size(Y))
       End Function Feuler
    End Interface

    If (Present(Tol)) Then
       h = Tol
    Else
       h = (Xfin - Xo)/DEFINTSTEPS_SP
    End If

    R = Xo
    Finales = Iniciales
    Do While (R .lt. Xfin)
       Finales = Finales + Tol * Feuler(R, Finales)
       R = R + Tol
    End Do

    Return
  End Function Euler_SP

!  *********************************************
!  *                                           *
  Function Rgnkta_SP(Iniciales,Xo,Xfin,Frgnkta,Tol) Result (Finales)
!  *                                           *
!  *********************************************
!  * Integrates numerically a set of ODE's 
!  * defined by the function Frgnkta by the 
!  * Runge-Kutta method. The initial conditions 
!  * are stored in iniciales, and integrates 
!  * between Xo and Xfin with tolerance Tol.
!  * 
!  * $ Last Mod: 21/09/2003 $
!  *********************************************

    Real (kind=SP), Intent (in) :: Iniciales(:), Xo, Xfin
    Real (kind=SP), Intent (in), Optional :: Tol
    Real (kind=SP) :: Finales(Size(Iniciales))
    Real (kind=SP) :: R, h
    Real (kind=SP) :: R1(Size(Iniciales)), R2(Size(Iniciales)), &
         & R3(Size(Iniciales)), R4(Size(Iniciales)), &
         & Y1(Size(Iniciales)), Aux(Size(Iniciales))
    
    Interface
       Function FRgnkta(X, Y) Result (Func)
         USE NumTypes

         Real (kind=SP), Intent (in) :: X, Y(:)
         Real (kind=SP) :: Func(Size(Y))
       End Function FRgnkta
    End Interface
    
    If (Present(Tol)) Then
       h = Tol
    Else
       h = (Xfin - Xo)/DEFINTSTEPS_SP
    End If

    Y1 = Iniciales
    
    R = Xo
    Do While (R.lt.Xfin - h)
       R1 = FRgnkta(R,Y1)
       Aux = Y1 + h*R1/3.0_SP
                
       R2 = FRgnkta(R + h/3.0_SP,Aux)
       Aux = Y1 + h*( -R1/3.0_SP + R2 )
       
       R3 = FRgnkta(R + 2.0_SP*h/3.0_SP,Aux)
       Aux = Y1 + h*( R1 - R2 + R3 )
       
       R4 = FRgnkta(R + h,Aux)
       Y1 = Y1 + (h/8.0_SP)*(R1 + 3.0_SP*R2 + 3.0_SP*R3 + R4)

       R = R + h
    End Do
    
    Finales = Y1
        
    Return
  End Function Rgnkta_SP


!  *********************************************
!  *                                           *
  Function IntQuadrilateral_DP(P1,P2,P3,P4,Fval) Result (ValI)
!  *                                           *
!  *********************************************
!  * 
!  *********************************************

    Real (kind=DP), Intent (in) :: P1(2), P2(2), P3(2), &
         & P4(2), Fval(4)
    Real (kind=DP) :: ValI

    Real (kind=DP), Parameter :: &
         & s(4) = (/-1.0_DP, 1.0_DP, 1.0_DP, -1.0_DP/), &
         & p(4) = (/-1.0_DP,-1.0_DP, 1.0_DP,  1.0_DP/)
    Real (kind=DP) :: alp(0:2), M(4,4), a(4), ap(4,2), Aux
    Integer :: I, J

    alp(0) = ( (P4(1)-P2(1))*(P1(2)-P3(2)) + (P3(1)-P1(1))*(P4(2)-P2(2)) )/8.0_DP
    alp(1) = ( (P4(1)-P3(1))*(P2(2)-P1(2)) + (P1(1)-P2(1))*(P4(2)-P3(2)) )/8.0_DP
    alp(2) = ( (P4(1)-P1(1))*(P2(2)-P3(2)) + (P3(1)-P2(1))*(P4(2)-P1(2)) )/8.0_DP

    ! First interpolate the function, to obtain the coef.
    ap(1,:) = P1(:)
    ap(2,:) = P2(:)
    ap(3,:) = P3(:)
    ap(4,:) = P4(:)

    M(:,1) = 1.0_DP
    M(:,2) = ap(:,1)
    M(:,3) = ap(:,2)
    M(:,4) = ap(:,1)*ap(:,2)

    a(:) = Fval(:)

!    Do I = 1, 4
!       Write(*,'(4F12.3)')(M(I,J), J=1, 4)
!    End Do
    CALL LUSolve(M, a)

    ValI = 4.0_DP*a(1)*alp(0)
    ValI = ValI + Sum( (a(2)*ap(:,1) + a(3)*ap(:,2)) * &
            & (alp(0) + alp(1)*s(:)/3.0_DP + alp(2)*p(:)/3.0_DP) )

    alp(0) = alp(0)/4.0_DP
    alp(1) = alp(1)/6.0_DP
    alp(2) = alp(2)/6.0_DP
    Aux = 0.0_DP
    Do I = 1, 4
       Do J = 1, 4
          Aux = Aux + ap(I,1)*ap(J,2)*&
               & ( alp(0)*(1.0_DP+s(I)*s(J)/3.0_DP) * &
               &          (1.0_DP+p(I)*p(J)/3.0_DP) + &
               &   alp(1)*(1.0_DP+p(I)*p(J)/3.0_DP) * &
               &          (s(I)+s(J))               + &
               &   alp(2)*(1.0_DP+s(I)*s(J)/3.0_DP) * &
               &          (p(I)+p(J))                 )
       End Do
    End Do
    ValI = ValI + Aux*a(4)

    Return
  End Function IntQuadrilateral_DP

!  *********************************************
!  *                                           *
  Function IntQuadrilateral_SP(P1,P2,P3,P4,Fval) Result (ValI)
!  *                                           *
!  *********************************************
!  * 
!  *********************************************

    Real (kind=SP), Intent (in) :: P1(2), P2(2), P3(2), &
         & P4(2), Fval(4)
    Real (kind=SP) :: ValI

    Real (kind=SP), Parameter :: &
         & s(4) = (/-1.0_SP, 1.0_SP, 1.0_SP, -1.0_SP/), &
         & p(4) = (/-1.0_SP,-1.0_SP, 1.0_SP,  1.0_SP/)
    Real (kind=SP) :: alp(0:2), M(4,4), a(4), ap(4,2), Aux
    Integer :: I, J

    alp(0) = ( (P4(1)-P2(1))*(P1(2)-P3(2)) + (P3(1)-P1(1))*(P4(2)-P2(2)) )/8.0_SP
    alp(1) = ( (P4(1)-P3(1))*(P2(2)-P1(2)) + (P1(1)-P2(1))*(P4(2)-P3(2)) )/8.0_SP
    alp(2) = ( (P4(1)-P1(1))*(P2(2)-P3(2)) + (P3(1)-P2(1))*(P4(2)-P1(2)) )/8.0_SP

    ap(1,:) = P1(:)
    ap(2,:) = P2(:)
    ap(3,:) = P3(:)
    ap(4,:) = P4(:)

    ! First interpolate the function, to obtain the coef.
    M(:,1) = 1.0_SP
    M(:,2) = ap(:,1)
    M(:,3) = ap(:,2)
    M(:,4) = ap(:,1)*ap(:,2)

    a(:) = Fval(:)
    CALL LUSolve(M, a)

    ! Now compute the integral of the bilinear interpoled function
    ValI = 4.0_SP*a(1)*alp(0)
    ValI = ValI + Sum( (a(2)*ap(:,1) + a(3)*ap(:,2)) * &
            & (alp(0) + alp(1)*s(:)/3.0_SP + alp(2)*p(:)/3.0_SP) )

    alp(0) = alp(0)/4.0_SP
    alp(1) = alp(1)/6.0_SP
    alp(2) = alp(2)/6.0_SP
    Aux = 0.0_SP
    Do I = 1, 4
       Do J = 1, 4
          Aux = Aux + ap(I,1)*ap(J,2)*&
               & ( alp(0)*(1.0_SP+s(I)*s(J)/3.0_SP) * &
               &          (1.0_SP+p(I)*p(J)/3.0_SP) + &
               &   alp(1)*(1.0_SP+p(I)*p(J)/3.0_SP) * &
               &          (s(I)+s(J))               + &
               &   alp(2)*(1.0_SP+s(I)*s(J)/3.0_SP) * &
               &          (p(I)+p(J))                 )
       End Do
    End Do
    ValI = ValI + Aux*a(4)

    Return
  End Function IntQuadrilateral_SP

End MODULE Integration

