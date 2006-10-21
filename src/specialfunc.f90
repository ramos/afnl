!
! MODULE to calculate values of special functions
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

! **************************************************
! *
MODULE SpecialFunc
! *
! ***************************************************
! *
! * Documentation available in:
! * 
! * http://lattice.ft.uam.es/perpag/alberto/
! * 
! ***************************************************

  USE NumTypes
  USE Constants
  USE Error

  ! Default tol
  Real (kind=DP), Parameter :: DEFTOL = 0.001_DP


  Interface ThetaChar
     Module Procedure ThetaChar_DPC, ThetaChar_DP
  End Interface

  Interface Hermite
     Module Procedure Hermite_DP, Hermite_SP
  End Interface

  Interface HermiteFunc
     Module Procedure HermiteFunc_DP, HermiteFunc_SP
  End Interface

  Interface Basis
     Module Procedure Basis_DPC, Basis_SPC
  End Interface


  Private DEFTOL, Theta3, ThetaChar_DPC, ThetaChar_SPC, Hermite_DP, &
       &Hermite_SP, HermiteFunc_SP, HermiteFunc_DP, Basis_DPC, Basis_SPC

CONTAINS

! ***********************************
! *
  Real (kind=DP) Function GammaLn(X)
! *
! ***********************************
! * Returns Ln(Gamma(X))
! * 
! * This routine follows the ideas 
! * of "Numerical Recipes in Fortran
! * 77, secon ed. (vol. 1)"
! ***********************************

    Real (kind=DP), Intent (In) :: X
    
    Real (kind=DP) :: Stp, Coef(6), Temp

    Save Coef, Stp
    Data Coef, Stp /76.18009172947146_DP, -86.50532032941677_DP, &
         & 24.01409824083091_DP, -1.231739572450155_DP,&
         & 0.1208650973866179E-2_DP, -0.5395239384953E-5_DP, &
         & 2.5066282746310005_DP/ 

    Temp = 1.000000000190015_DP
    Do I = 1, 6
       Temp = Temp + Coef(I)/(X - 1.D0 + Real(I,kind=DP)) 
    End Do

    GammaLn = Log(Stp * Temp) + (X-0.5D0)*Log(X + 4.5_DP ) - (X + 4.5_DP)
   
    Return
  End Function GammaLn
  
! *************************************
! *
  Recursive Complex (kind=DPC) Function Theta3(z, tau, Tol) Result (value)
! *
! *************************************
! * Calculate the value of the third 
! * Jacobi theta function of parameter 
! * tau in the complex point z, with 
! * precision Tol.
! *
! * We use the series defined in 
! * "A course in modern analysis,  
! * Watson & Wittaker". In particular our 
! * Theta functions have quasi-periods
! * pi and tau*pi
! * 
! * This is a low level function
! ************************************* 


    Complex (kind=DPC), Intent (in) :: z, tau
    Real (kind=DP), Intent (in) :: Tol
    
    Complex (kind=DPC) :: z2, ThetaOld ! Real value, with modulus.
    
    
    Value = (1.D0, 0.D0)
    ThetaOld = Value + 2*Tol
    k = 0

    ! The series that defines the Theta Functions
    ! converges very fast for tau >= 1.
    ! If tau >= 1, then the use the power series
    ! expansion for the Jacobi Function, else use
    ! Jacobi imaginary transformation to make 
    ! tau > 1
    If (Aimag(tau) >= 1.D0) Then
       Do While (Abs(Value - ThetaOld) / Abs(Value) > Tol)
          ThetaOld = Value
          k = k + 1
          
          Value = Value + 2.D0 * Exp((0.D0, 1.D0) * PI_DP * tau * k**2) * & 
               & Cos(2*k*z)
       End Do
    Else
       Value = Sqrt((0.D0, -1.D0)*tau) * &
            & Exp(z**2/((0.D0,1.D0)*PI_DP*tau)) * &
            & Theta3(z/tau, (-1.D0, 0.D0)/tau, Tol)
    End If

    Return
  End Function Theta3  
  
! *************************************
! *
  Complex (kind=DPC) Function ThetaChar_DPC(a, b, z, tau, Prec)
! *
! *************************************
! *
! * Calculates the Theta functions with
! * characteristics a,b of the complex 
! * variable z, with precision Prec.
! *
! *************************************

    Complex (kind=DPC), Intent (in) :: a, b, z, tau
    Real (kind=DP), Intent (in), Optional :: Prec

    Real (kind=DP) :: Tol
    
    If (Present(Prec)) Then
       Tol = Prec
    Else
       Tol = DEFTOL
    End If
    
    If (Aimag(tau) <= 0.D0 ) Then
       CALL Abort("Theta_Char_DPC", &
            & "Theta: tau parameter *must* have Imag(tau) > 0.")
    End If

    ThetaChar_DPC = exp( UnitImag_DPC*a*(a*PI_DP*tau + 2.0_DP*(Z &
         & + b*PI_DP))) * &
         & theta3(Z + PI_DP*(tau*a + Cmplx(b,kind=DPC)), tau, Tol )

    Return
  End Function ThetaChar_DPC

! *************************************
! *
  Complex (kind=DPC) Function ThetaChar_DP(a, b, z, tau, Prec)
! *
! *************************************
! *
! * Calculates the Theta functions with
! * characteristics a,b of the complex 
! * variable z, with precision Prec.
! *
! *************************************

    Complex (kind=DPC), Intent (in) :: z, tau
    Real (Kind=DP), Intent (in) :: a, b
    Real (kind=DP), Intent (in), Optional :: Prec

    Real (kind=DP) :: Tol
    
    If (Present(Prec)) Then
       Tol = Prec
    Else
       Tol = DEFTOL
    End If
    
    If (Aimag(tau) <= 0.D0 ) Then
       CALL Abort("Theta_Char_DP", &
            & "Theta: tau parameter *must* have Imag(tau) > 0.")
    End If

    ThetaChar_DP = exp( UnitImag_DPC*a*(a*PI_DP*tau + 2.0_DP*(Z&
         & +Cmplx(b*PI_DP,kind=DPC)))) * &
         & theta3(Z + PI_DP*(tau*a + Cmplx(b,kind=DPC)), tau, Tol )

    Return
  End Function ThetaChar_DP


! *************************************
! *
  Complex (kind=DPC) Function Theta(i, z, tau, Prec)
! *
! *************************************
! *
! * This is an "interface" to the call of 
! * different theta functions. The calculus
! * of all the theta functions relays 
! * on the calculus of theta3.
! * 
! * We use periodicity properties to 
! * We also try to check some errors,
! * like tau < 1, or a call to the 
! * 56 theta function.
! *************************************
    Integer, Intent (in) :: i
    Complex (kind=DPC), Intent (in) :: z, tau
    Real (kind=DP), Intent (in), Optional :: Prec
    
    Real (kind=DP) :: Tol
    Complex (kind=DPC) :: zmod, Factor ! z with the necassary mod op.
    
    If (Present(Prec)) Then
       Tol = Prec
    Else
       Tol = DEFTOL
    End If
    
    If (Aimag(tau) <= 0.D0 ) Then
       CALL Abort("Theta", "tau parameter *must* have Imag(tau) > 0.")
    End If
    
!    nrest = Int(Aimag(z)/(Pi*Aimag(Tau) ) )
!    Write(*,*)'Nrest: ', nrest
!    Factor = exp(-(0.0_DP,1.0_DP) * (2.0_DP * z + Pi*tau))

!    zmod = zmod - Cmplx(nrest*Pi*Aimag(Tau), kind=DPC)

    If (i == 1) Then
       Theta = -UnitImag_DPC * exp(UnitImag_DPC*(z + PI_DP*tau/4.0_DP)) * &
            & Theta3(z + HALFPI_DP*(tau - Cmplx(1.0_DP, kind=DPC)), tau, Tol)
    Else If (i == 2) Then
       Theta =  exp(UnitImag_DPC*(z + PI_DP*tau/4.0_DP)) * &
            & Theta3(z + HALFPI_DP*tau, tau, Tol)
    Else If (i == 3) Then 
       zmod = (1.0_DP, 0.0_DP) * ( Mod(Real(z), PI_DP) ) + &
            & UnitImag_DPC * Aimag(z)
       Theta = Theta3(zmod, tau, Tol)
    Else If (i == 4) Then
       Theta = Theta3(z - Cmplx(HALFPI_DP, kind=DPC), tau, Tol)
    Else
       CALL Abort("Theta", "I only know 4 theta functions!!")
    End If
    
    Return
  End Function Theta

! *************************************
! *
  Real (kind=DP) Function HermiteFunc_DP(n, X, Dval)
! *
! *************************************
! *
! * Calculates the n-th degree Hermite 
! * Function in the point X. If Dval is 
! * specified, it returns the value of the
! * derivative. DP Version
! *
! *************************************

    Integer, Intent (in) :: n
    Real (kind=DP), Intent (in) :: X
    Real (kind=DP), Intent (out), Optional :: Dval

    Real (kind=DP) :: Dnorm

    Dnorm = 1.0_DP / &
         & Sqrt(exp(GammaLn(Real(n+1,kind=DP))) * 2.0_DP**n*SRPI_DP) 
    If (Present(Dval)) Then
       HermiteFunc_DP = Dnorm * (exp(-X**2/2.0_DP) * Hermite(n,X,Dval))
       Dval = Dnorm * (exp(-X**2/2.0_DP) * Dval - X*HermiteFunc_DP)
    Else
       HermiteFunc_DP = Dnorm * (exp(-X**2/2.0_DP) * Hermite(n,X))
    End If

    Return
  End Function HermiteFunc_DP

! *************************************
! *
  Real (kind=SP) Function HermiteFunc_SP(n, X, Dval)
! *
! *************************************
! *
! * Calculates the n-th degree Hermite 
! * Function in the point X. If Dval is 
! * specified, it returns the value of the
! * derivative. SP Version
! *
! *************************************

    Integer, Intent (in) :: n
    Real (kind=SP), Intent (in) :: X
    Real (kind=SP), Intent (out), Optional :: Dval

    Real (kind=SP) :: Dnorm

    Dnorm = 1.0_SP / &
         & Sqrt(Real(exp(GammaLn(Real(n+1,kind=DP))),kind=SP) &
         & * 2.0_SP**n*SRPI_SP) 
    If (Present(Dval)) Then
       HermiteFunc_SP = Dnorm * (exp(-X**2/2.0_SP) * Hermite(n,X,Dval))
       Dval = Dnorm * (exp(-X**2/2.0_SP) * Dval - X*HermiteFunc_SP)
    Else
       HermiteFunc_SP = Dnorm * (exp(-X**2/2.0_SP) * Hermite(n,X))
    End If

    Return
  End Function HermiteFunc_SP

! *************************************
! *
  Real (kind=DP) Function Hermite_DP(n, X, Dval)
! *
! *************************************
! *
! * Calculates the n-th degree Hermite 
! * polynomial in the point X. If Dval is 
! * specified, it returns the value of the
! * derivative. DP Version
! *************************************

    Integer, Intent (in) :: n
    Real (kind=DP), Intent (in) :: X
    Real (kind=DP), Intent (out), Optional :: Dval

    Real (Kind=DP) :: Y0, Y1

    If (n < 0) Then
       CALL Abort("Hermite_DP", "I dont know what Hn(X) with n<0 is!")
    End If

    If (n == 0) Then
       Hermite_DP = 1.0_DP
       Return
    Else
       Y0 = 1.0_DP
    End If

    If (n == 1) Then
       Hermite_DP = 2.0_DP * X
       Return
    Else
       Y1 = 2.0_DP * X
    End If

    Do I = 2, n
       Hermite_DP = 2.0_DP*X*Y1 - 2.0_DP*(Real(I,kind=DP) - 1.0_DP)*Y0
       Y0 = Y1
       Y1 = Hermite_DP
    End Do

    If (Present(Dval)) Dval = 2.0_DP*Real(n, kind=DP)*Y0

    Return
  End Function Hermite_DP

! *************************************
! *
  Real (kind=SP) Function Hermite_SP(n, X, Dval)
! *
! *************************************
! *
! * Calculates the n-th degree Hermite 
! * polynomial in the point X. If Dval is 
! * specified, it returns the value of the
! * derivative. SP Version
! *************************************

    Integer, Intent (in) :: n
    Real (kind=SP), Intent (in) :: X
    Real (kind=SP), Intent (out), Optional :: Dval

    Real (Kind=SP) :: Y0, Y1

    If (n < 0) Then
       CALL Abort("Hermite_SP", "I dont know what Hn(X) with n<0 is!")
    End If

    If (n == 0) Then
       Hermite_SP = 1.0_SP
       Return
    Else
       Y0 = 1.0_SP
    End If

    If (n == 1) Then
       Hermite_SP = 2.0_SP * X
       Return
    Else
       Y1 = 2.0_SP * X
    End If

    Do I = 2, n
       Hermite_SP = 2.0_SP*X*Y1 - 2.0_SP*(Real(I,kind=SP) - 1.0_SP)*Y0
       Y0 = Y1
       Y1 = Hermite_SP
    End Do

    If (Present(Dval)) Dval = 2.0_SP*Real(n, kind=SP)*Y0

    Return
  End Function Hermite_SP

! *************************************
! *
  Complex (kind=DPC) Function Basis_DPC(X1, X2, n, s, q, itau, Prec)
! *
! *************************************
! *
! * Computes the value of the basis element 
! * |n,s> in the point (X1,X2) With optional
! * precision Prec.
! * 
! * X1 and X2 are assumed to have period 1.
! * 
! *************************************

    Integer, Intent (in) :: n, s, q
    Real (kind=DP), Intent (in) :: X1, X2, itau
    Real (kind=DP), Intent (in), Optional :: Prec

    Complex (Kind=DPC) :: OldVal, w
    Real (Kind=DP) :: u
    Integer :: Kit

    If (Present(Prec)) Then
       Tol = Prec
    Else
       Tol = DEFTOL
    End If
    
    w = exp(TWOPI_IMAG_DPC * q * X1)
    Kit = 0
    u = Sqrt(TWOPI_DP*q*itau) * ( X2 + &
         & Real(s,kind=DP)/Real(q, kind=DP) + Real(Kit,kind=DP))
    Basis_DPC = HermiteFunc(n,u)
    OldVal = Basis_DPC + Cmplx(10.0_DP)
    Do While (Abs(Basis_DPC-OldVal)/Abs(Basis_DPC+OldVal) > Tol)
       OldVal = Basis_DPC
       Kit = Kit + 1

       u = Sqrt(TWOPI_DP*q*itau) * ( X2 + &
            & Real(s,kind=DP)/Real(q, kind=DP) + Real(Kit,kind=DP))
       Basis_DPC = Basis_DPC + HermiteFunc(n,u) * &
            &w**Kit
       u = Sqrt(TWOPI_DP*q*itau) * ( X2 + &
            & Real(s,kind=DP)/Real(q, kind=DP) - Real(Kit,kind=DP))
       Basis_DPC = Basis_DPC + HermiteFunc(n,u) * &
            &w**(-Kit)

    End Do

    Basis_DPC = exp(TWOPI_IMAG_DPC*s*X1) * exp(PI_IMAG_DPC*q*X1*X2) *&
         & Basis_DPC

    Return
  End Function Basis_DPC

! *************************************
! *
  Complex (kind=SPC) Function Basis_SPC(X1, X2, n, s, q, itau, Prec)
! *
! *************************************
! *
! * Computes the value of the basis element 
! * |n,s> in the point (X1,X2) With optional
! * precision Prec.
! * 
! * X1 and X2 are assumed to have period 1.
! * 
! *************************************

    Integer, Intent (in) :: n, s, q
    Real (kind=SP), Intent (in) :: X1, X2, itau
    Real (kind=SP), Intent (in), Optional :: Prec

    Complex (Kind=SPC) :: OldVal, w
    Real (Kind=SP) :: u
    Integer :: Kit

    If (Present(Prec)) Then
       Tol = Prec
    Else
       Tol = DEFTOL
    End If
    
    w = exp(TWOPI_IMAG_SPC * q * X1)
    Kit = 0
    u = Sqrt(TWOPI_SP*q*itau) * ( X2 + &
         & Real(s,kind=SP)/Real(q, kind=SP) + Real(Kit,kind=SP))
    Basis_SPC = HermiteFunc(n,u)
    OldVal = Basis_SPC + Cmplx(10.0_SP)
    Do While (Abs(Basis_SPC-OldVal)/Abs(Basis_SPC+OldVal) > Tol)
       OldVal = Basis_SPC
       Kit = Kit + 1

       u = Sqrt(TWOPI_SP*q*itau) * ( X2 + &
            & Real(s,kind=SP)/Real(q, kind=SP) + Real(Kit,kind=SP))
       Basis_SPC = Basis_SPC + HermiteFunc(n,u) * &
            &w**Kit
       u = Sqrt(TWOPI_SP*q*itau) * ( X2 + &
            & Real(s,kind=SP)/Real(q, kind=SP) - Real(Kit,kind=SP))
       Basis_SPC = Basis_SPC + HermiteFunc(n,u) * &
            &w**(-Kit)

    End Do

    Basis_SPC = exp(TWOPI_IMAG_SPC*s*X1) * exp(PI_IMAG_SPC*q*X1*X2) *&
         & Basis_SPC

    Return
  End Function Basis_SPC


End MODULE SpecialFunc
