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

  IMPLICIT NONE

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

  Interface Factorial
     Module Procedure Factorial_DP
  End Interface

  Interface erf
     Module Procedure erf_DP, erf_SP
  End Interface

  Interface erfc
     Module Procedure erfc_DP, erfc_SP
  End Interface

  Interface Legendre
     Module Procedure Legendre_DP, Legendre_SP
  End Interface

  Interface SphericalHarmonic
     Module Procedure SphericalHarmonic_DP, SphericalHarmonic_SP
  End Interface

  Interface inverf
     Module Procedure inverf_DP, inverf_SP
  End Interface

  Private DEFTOL, Theta3, ThetaChar_DPC, Hermite_DP, &
       &Hermite_SP, HermiteFunc_SP, HermiteFunc_DP, Basis_DPC,&
       & Basis_SPC, Factorial_DP, erf_DP, erf_SP, &
       & erfc_DP, erfc_SP, Legendre_DP, Legendre_SP, &
       & SphericalHarmonic_DP, SphericalHarmonic_SP, &
       & inverf_DP, inverf_SP

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
    Integer :: I

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
    
    Complex (kind=DPC) :: ThetaOld ! Real value, with modulus.
    Integer :: k
    
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
    Complex (kind=DPC) :: zmod ! z with the necassary mod op.
    
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
    Integer :: I

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
    Integer :: I

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
    Real (Kind=DP) :: u, Tol
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
    Real (Kind=SP) :: u, Tol
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

! *************************************
! *
  Real (kind=DP) Function Factorial_DP(N)
! *
! *************************************
! * Compute the factorial of an integer
! * number. DP version.
! *************************************

    Integer, Parameter :: Nbrute = 51
    Integer, Intent (in) :: N

    Real (kind=DP), Save :: SFact(Nbrute)
    Integer, Save :: Ncompt
    
    Integer :: I


    If ( (N == 0).or.(N == 1) ) Then
       Factorial_DP = 1.0_DP
       Return
    Else If (N < 0) Then
       CALL abort('Factorial_DP', &
            & 'I dont know to compute negative factorials')
    End If

    SFact(1) = 1.0_DP
    Ncompt = 1
    If (N < Ncompt) Then
       Factorial_DP = SFact(N)
    Else If (N <= Nbrute) Then
       Do I = Ncompt+1, N
          SFact(I) = Real(I,kind=DP)*SFact(I-1)
       End Do
       Ncompt = N
       Factorial_DP = SFact(N)
    Else
       Factorial_DP = Exp(GammaLn(Real(N,kind=SP)+1.0_DP))
    End If
          
    Return
  End Function Factorial_DP

! *************************************
! *
  Real (kind=DP) Function erf_DP(X)
! *
! *************************************
! * Compute the error function for real
! * arguments. DP version.
! * Based on implementation by W. J. Cody
! * supossed to be exact up to double 
! * precision.
! *************************************

    Real (kind=DP), Intent (in) :: X

    Real (kind=DP) :: Xabs, Rnum, Rden, Xsqr
    Real (kind=DP) :: Dt1 = 0.46875_DP, Dt2 = 4.0_DP, Xmax = 6.0_DP
    Real (kind=DP) :: Cfa1(5), Cfb1(4), Cfa2(9), Cfb2(8), &
         & Cfa3(6), Cfb3(5)

    Integer :: I

    ! Coeficients for first interval computation
    Data Cfa1 /3.16112374387056560_DP, 1.13864154151050156E2_DP, &
         & 3.77485237685302021E2_DP, 3.20937758913846947E3_DP, &
         & 1.85777706184603153E-1_DP/
    Data Cfb1 /2.36012909523441209E1_DP, 2.44024637934444173E2_DP, &
         & 1.28261652607737228E3_DP, 2.84423683343917062E3_DP/

    ! Coeficients for second interval computation
    Data Cfa2 /5.64188496988670089E-1_DP,8.88314979438837594_DP, &
         & 6.61191906371416295E1_DP,2.98635138197400131E2_DP, &
         & 8.81952221241769090D02,1.71204761263407058E3_DP, &
         & 2.05107837782607147E3_DP,1.23033935479799725E3_DP, &
         & 2.15311535474403846E-8_DP/
    Data Cfb2 / 1.57449261107098347E1_DP,1.17693950891312499E2_DP, &
         & 5.37181101862009858E2_DP,1.62138957456669019E3_DP, &
         & 3.29079923573345963E3_DP,4.36261909014324716E3_DP, &
         & 3.43936767414372164E3_DP,1.23033935480374942E3_DP/

    ! Coeficients for last interval computation
    Data Cfa3 /3.05326634961232344E-1_DP,3.60344899949804439E-1_DP, &
         & 1.25781726111229246E-1_DP,1.60837851487422766E-2_DP, &
         & 6.58749161529837803E-4_DP,1.63153871373020978E-2_DP/  

    Data Cfb3/2.56852019228982242_DP,1.87295284992346047_DP, &
         & 5.27905102951428412E-1_DP,6.05183413124413191E-2_DP, &
         & 2.33520497626869185E-3_DP/


    Xabs = Abs(X)
    Xsqr = X**2
    If (Xabs < Dt1) Then
       Rnum = Cfa1(5)*Xsqr
       Rden = Xsqr
       Do I = 1, 3
          Rnum = (Rnum + Cfa1(I)) * Xsqr
          Rden = (Rden + Cfb1(I)) * Xsqr
       End Do
       erf_DP = Sign(X * (Rnum + Cfa1(4)) / (Rden + Cfb1(4)), X)
       Return
    Else If (Xabs < Dt2) Then
       Rnum = Cfa2(9)*Xabs
       Rden = Xabs
       Do I = 1, 7
          Rnum = (Rnum + Cfa2(I)) * Xabs
          Rden = (Rden + Cfb2(I)) * Xabs
       End Do
       erf_DP = Sign( 1.0_DP - Exp(-X**2) * &
            & (Rnum + Cfa2(8)) / (Rden + Cfb2(8)), X)
       Return
    Else If (Xabs < Xmax) Then
       Xsqr = 1.0_DP / X**2
       Rnum = Cfa3(6)*Xsqr
       Rden = Xsqr
       Do I = 1, 4
          Rnum = (Rnum + Cfa3(I)) * Xsqr
          Rden = (Rden + Cfb3(I)) * Xsqr
       End Do
       erf_DP = Xsqr * (Rnum + Cfa3(5)) / (Rden + Cfb3(5))
       erf_DP = ( 1.0_DP/Srpi_DP - erf_DP ) / Xabs
       erf_DP = Sign( 1.0_DP - Exp(-X**2)*erf_DP, X)
       Return
    Else
       erf_DP = Sign(1.0_DP, X)
       Return
    End If


    Return
  End Function erf_DP

! *************************************
! *
  Real (kind=DP) Function erfc_DP(X)
! *
! *************************************
! * Compute the complementary error function 
! * for real arguments. DP version.
! * Based on implementation by W. J. Cody
! * supossed to be exact up to double 
! * precision.
! *************************************

    Real (kind=DP), Intent (in) :: X

    Real (kind=DP) :: Xabs, Rnum, Rden, Xsqr
    Real (kind=DP) :: Dt1 = 0.46875_DP, Dt2 = 4.0_DP, Xmax = 30.0_DP
    Real (kind=DP) :: Cfa1(5), Cfb1(4), Cfa2(9), Cfb2(8), &
         & Cfa3(6), Cfb3(5)

    Integer :: I

    ! Coeficients for first interval computation
    Data Cfa1 /3.16112374387056560_DP, 1.13864154151050156E2_DP, &
         & 3.77485237685302021E2_DP, 3.20937758913846947E3_DP, &
         & 1.85777706184603153E-1_DP/
    Data Cfb1 /2.36012909523441209E1_DP, 2.44024637934444173E2_DP, &
         & 1.28261652607737228E3_DP, 2.84423683343917062E3_DP/

    ! Coeficients for second interval computation
    Data Cfa2 /5.64188496988670089E-1_DP,8.88314979438837594_DP, &
         & 6.61191906371416295E1_DP,2.98635138197400131E2_DP, &
         & 8.81952221241769090D02,1.71204761263407058E3_DP, &
         & 2.05107837782607147E3_DP,1.23033935479799725E3_DP, &
         & 2.15311535474403846E-8_DP/
    Data Cfb2 / 1.57449261107098347E1_DP,1.17693950891312499E2_DP, &
         & 5.37181101862009858E2_DP,1.62138957456669019E3_DP, &
         & 3.29079923573345963E3_DP,4.36261909014324716E3_DP, &
         & 3.43936767414372164E3_DP,1.23033935480374942E3_DP/

    ! Coeficients for last interval computation
    Data Cfa3 /3.05326634961232344E-1_DP,3.60344899949804439E-1_DP, &
         & 1.25781726111229246E-1_DP,1.60837851487422766E-2_DP, &
         & 6.58749161529837803E-4_DP,1.63153871373020978E-2_DP/  

    Data Cfb3/2.56852019228982242_DP,1.87295284992346047_DP, &
         & 5.27905102951428412E-1_DP,6.05183413124413191E-2_DP, &
         & 2.33520497626869185E-3_DP/


    Xabs = Abs(X)
    Xsqr = X**2
    If (Xabs < Dt1) Then
       Rnum = Cfa1(5)*Xsqr
       Rden = Xsqr
       Do I = 1, 3
          Rnum = (Rnum + Cfa1(I)) * Xsqr
          Rden = (Rden + Cfb1(I)) * Xsqr
       End Do
       erfc_DP = 1.0_DP - Sign(X * (Rnum + Cfa1(4)) / (Rden + Cfb1(4)), X)
       Return
    Else If (Xabs < Dt2) Then
       Rnum = Cfa2(9)*Xabs
       Rden = Xabs
       Do I = 1, 7
          Rnum = (Rnum + Cfa2(I)) * Xabs
          Rden = (Rden + Cfb2(I)) * Xabs
       End Do
       erfc_DP = Exp(-X**2) * (Rnum + Cfa2(8)) / (Rden + Cfb2(8))
    Else If (Xabs < Xmax) Then
       Xsqr = 1.0_DP / X**2
       Rnum = Cfa3(6)*Xsqr
       Rden = Xsqr
       Do I = 1, 4
          Rnum = (Rnum + Cfa3(I)) * Xsqr
          Rden = (Rden + Cfb3(I)) * Xsqr
       End Do
       erfc_DP = Xsqr * (Rnum + Cfa3(5)) / (Rden + Cfb3(5))
       erfc_DP = ( 1.0_DP/Srpi_DP - erfc_DP ) / Xabs
       erfc_DP = Exp(-X**2)*erfc_DP
    Else
       erfc_DP = 0.0_DP
    End If

    If (X < 0.0_DP) erfc_DP = 2.0_DP - erfc_DP    

    Return
  End Function erfc_DP


! *************************************
! *
  Real (kind=SP) Function erf_SP(X)
! *
! *************************************
! * Compute the error function for real
! * arguments. SP version.
! * Based on implementation by W. J. Cody
! * supossed to be exact up to double 
! * precision.
! *************************************

    Real (kind=SP), Intent (in) :: X

    Real (kind=SP) :: Xabs, Rnum, Rden, Xsqr
    Real (kind=SP) :: Dt1 = 0.46875_SP, Dt2 = 4.0_SP, Xmax = 6.0_SP
    Real (kind=SP) :: Cfa1(5), Cfb1(4), Cfa2(9), Cfb2(8), &
         & Cfa3(6), Cfb3(5)

    Integer :: I

    ! Coeficients for first interval computation
    Data Cfa1 /3.16112374387056560_SP, 1.13864154151050156E2_SP, &
         & 3.77485237685302021E2_SP, 3.20937758913846947E3_SP, &
         & 1.85777706184603153E-1_SP/
    Data Cfb1 /2.36012909523441209E1_SP, 2.44024637934444173E2_SP, &
         & 1.28261652607737228E3_SP, 2.84423683343917062E3_SP/

    ! Coeficients for second interval computation
    Data Cfa2 /5.64188496988670089E-1_SP,8.88314979438837594_SP, &
         & 6.61191906371416295E1_SP,2.98635138197400131E2_SP, &
         & 8.81952221241769090D02,1.71204761263407058E3_SP, &
         & 2.05107837782607147E3_SP,1.23033935479799725E3_SP, &
         & 2.15311535474403846E-8_SP/
    Data Cfb2 / 1.57449261107098347E1_SP,1.17693950891312499E2_SP, &
         & 5.37181101862009858E2_SP,1.62138957456669019E3_SP, &
         & 3.29079923573345963E3_SP,4.36261909014324716E3_SP, &
         & 3.43936767414372164E3_SP,1.23033935480374942E3_SP/

    ! Coeficients for last interval computation
    Data Cfa3 /3.05326634961232344E-1_SP,3.60344899949804439E-1_SP, &
         & 1.25781726111229246E-1_SP,1.60837851487422766E-2_SP, &
         & 6.58749161529837803E-4_SP,1.63153871373020978E-2_SP/  

    Data Cfb3/2.56852019228982242_SP,1.87295284992346047_SP, &
         & 5.27905102951428412E-1_SP,6.05183413124413191E-2_SP, &
         & 2.33520497626869185E-3_SP/


    Xabs = Abs(X)
    Xsqr = X**2
    If (Xabs < Dt1) Then
       Rnum = Cfa1(5)*Xsqr
       Rden = Xsqr
       Do I = 1, 3
          Rnum = (Rnum + Cfa1(I)) * Xsqr
          Rden = (Rden + Cfb1(I)) * Xsqr
       End Do
       erf_SP = Sign(X * (Rnum + Cfa1(4)) / (Rden + Cfb1(4)), X)
       Return
    Else If (Xabs < Dt2) Then
       Rnum = Cfa2(9)*Xabs
       Rden = Xabs
       Do I = 1, 7
          Rnum = (Rnum + Cfa2(I)) * Xabs
          Rden = (Rden + Cfb2(I)) * Xabs
       End Do
       erf_SP = Sign( 1.0_SP - Exp(-X**2) * &
            & (Rnum + Cfa2(8)) / (Rden + Cfb2(8)), X)
       Return
    Else If (Xabs < Xmax) Then
       Xsqr = 1.0_SP / X**2
       Rnum = Cfa3(6)*Xsqr
       Rden = Xsqr
       Do I = 1, 4
          Rnum = (Rnum + Cfa3(I)) * Xsqr
          Rden = (Rden + Cfb3(I)) * Xsqr
       End Do
       erf_SP = Xsqr * (Rnum + Cfa3(5)) / (Rden + Cfb3(5))
       erf_SP = ( 1.0_SP/Srpi_SP - erf_SP ) / Xabs
       erf_SP = Sign( 1.0_SP - Exp(-X**2)*erf_SP, X)
       Return
    Else
       erf_SP = Sign(1.0_SP, X)
       Return
    End If


    Return
  End Function erf_SP

! *************************************
! *
  Real (kind=SP) Function erfc_SP(X)
! *
! *************************************
! * Compute the complementary error function 
! * for real arguments. SP version.
! * Based on implementation by W. J. Cody
! * supossed to be exact up to double 
! * precision.
! *************************************

    Real (kind=SP), Intent (in) :: X

    Real (kind=SP) :: Xabs, Rnum, Rden, Xsqr
    Real (kind=SP) :: Dt1 = 0.46875_SP, Dt2 = 4.0_SP, Xmax = 30.0_SP
    Real (kind=SP) :: Cfa1(5), Cfb1(4), Cfa2(9), Cfb2(8), &
         & Cfa3(6), Cfb3(5)

    Integer :: I

    ! Coeficients for first interval computation
    Data Cfa1 /3.16112374387056560_SP, 1.13864154151050156E2_SP, &
         & 3.77485237685302021E2_SP, 3.20937758913846947E3_SP, &
         & 1.85777706184603153E-1_SP/
    Data Cfb1 /2.36012909523441209E1_SP, 2.44024637934444173E2_SP, &
         & 1.28261652607737228E3_SP, 2.84423683343917062E3_SP/

    ! Coeficients for second interval computation
    Data Cfa2 /5.64188496988670089E-1_SP,8.88314979438837594_SP, &
         & 6.61191906371416295E1_SP,2.98635138197400131E2_SP, &
         & 8.81952221241769090D02,1.71204761263407058E3_SP, &
         & 2.05107837782607147E3_SP,1.23033935479799725E3_SP, &
         & 2.15311535474403846E-8_SP/
    Data Cfb2 / 1.57449261107098347E1_SP,1.17693950891312499E2_SP, &
         & 5.37181101862009858E2_SP,1.62138957456669019E3_SP, &
         & 3.29079923573345963E3_SP,4.36261909014324716E3_SP, &
         & 3.43936767414372164E3_SP,1.23033935480374942E3_SP/

    ! Coeficients for last interval computation
    Data Cfa3 /3.05326634961232344E-1_SP,3.60344899949804439E-1_SP, &
         & 1.25781726111229246E-1_SP,1.60837851487422766E-2_SP, &
         & 6.58749161529837803E-4_SP,1.63153871373020978E-2_SP/  

    Data Cfb3/2.56852019228982242_SP,1.87295284992346047_SP, &
         & 5.27905102951428412E-1_SP,6.05183413124413191E-2_SP, &
         & 2.33520497626869185E-3_SP/


    Xabs = Abs(X)
    Xsqr = X**2
    If (Xabs < Dt1) Then
       Rnum = Cfa1(5)*Xsqr
       Rden = Xsqr
       Do I = 1, 3
          Rnum = (Rnum + Cfa1(I)) * Xsqr
          Rden = (Rden + Cfb1(I)) * Xsqr
       End Do
       erfc_SP = 1.0_SP - Sign(X * (Rnum + Cfa1(4)) / (Rden + Cfb1(4)), X)
       Return
    Else If (Xabs < Dt2) Then
       Rnum = Cfa2(9)*Xabs
       Rden = Xabs
       Do I = 1, 7
          Rnum = (Rnum + Cfa2(I)) * Xabs
          Rden = (Rden + Cfb2(I)) * Xabs
       End Do
       erfc_SP = Exp(-X**2) * (Rnum + Cfa2(8)) / (Rden + Cfb2(8))
    Else If (Xabs < Xmax) Then
       Xsqr = 1.0_SP / X**2
       Rnum = Cfa3(6)*Xsqr
       Rden = Xsqr
       Do I = 1, 4
          Rnum = (Rnum + Cfa3(I)) * Xsqr
          Rden = (Rden + Cfb3(I)) * Xsqr
       End Do
       erfc_SP = Xsqr * (Rnum + Cfa3(5)) / (Rden + Cfb3(5))
       erfc_SP = ( 1.0_SP/Srpi_SP - erfc_SP ) / Xabs
       erfc_SP = Exp(-X**2)*erfc_SP
    Else
       erfc_SP = 0.0_SP
    End If

    If (X < 0.0_SP) erfc_SP = 2.0_SP - erfc_SP

    Return
  End Function erfc_SP

! *************************************
! *
  Real (kind=DP) Function Legendre_DP(l, m, X) Result (Leg)
! *
! *************************************
! * Compute the value of the associated Legendre 
! * polinomial (l,m) at X.
! *************************************

    Real (kind=DP), Intent (in) :: X
    Integer, Intent (in) :: l, m

    Integer :: I, mabs
    Real (kind=DP) :: Val(0:1), NegFac, Aux

    Leg = 0.0_DP
    mabs = Abs(m)
    If (m < 0) Then
       If (l == mabs) Then
          NegFac = (-1)**mabs / &
               & Factorial(l-m)
       Else
          NegFac = (-1)**mabs * &
               & Exp(GammaLn(Real(l-mabs,kind=DP)) - & 
               & GammaLn(Real(l+mabs,kind=DP)) )
       End If
    Else
       NegFac = 1.0_DP
    End If

    ! Compute Pm,m
    Val(0) = 1.0_DP
    Do I = 3, 2*mabs-1, 2
       Val(0) = Val(0)*Real(I,kind=DP)
    End Do
    Val(0) = Val(0)*(-Sqrt(1.0_DP-X**2))**mabs

    If (l==mabs) Then
       Leg = NegFac * Val(0)
       Return
    End If

    ! Compute Pm,m+1
    Val(1) = X*(2.0_DP*mabs+1.0_DP)*Val(0)
    If (l == mabs+1) Then
       Leg = NegFac * Val(1)
       Return
    End If

    ! For higher l, use the recursion relation
    Do I = mabs+2, l
       Aux = (2.0_DP*(I-1)+1.0_DP) * X * Val(1) - &
            & (I-1+mabs) * Val(0)
       Aux = Aux / Real(I-mabs,kind=DP)
       Val(0) = Val(1)
       Val(1) = Aux
    End Do
    Leg = NegFac * Val(1)

    Return
  End Function Legendre_DP

! *************************************
! *
  Real (kind=SP) Function Legendre_SP(l, m, X) Result (Leg)
! *
! *************************************
! * Compute the value of the associated Legendre 
! * polinomial (l,m) at X.
! *************************************

    Real (kind=SP), Intent (in) :: X
    Integer, Intent (in) :: l, m

    Integer :: I, mabs
    Real (kind=SP) :: Val(0:1), NegFac, Aux

    Leg = 0.0_SP
    mabs = Abs(m)
    If (m < 0) Then
       If (l == mabs) Then
          NegFac = (-1)**mabs / &
               & Factorial(l-m)
       Else
          NegFac = (-1)**mabs * &
               & Real(Exp(GammaLn(Real(l-mabs,kind=DP)) - & 
               & GammaLn(Real(l+mabs,kind=DP)) ),kind=SP)
       End If
    Else
       NegFac = 1.0_SP
    End If

    ! Compute Pm,m
    Val(0) = 1.0_SP
    Do I = 3, 2*mabs-1, 2
       Val(0) = Val(0)*Real(I,kind=SP)
    End Do
    Val(0) = Val(0)*(-Sqrt(1.0_SP-X**2))**mabs

    If (l==mabs) Then
       Leg = NegFac * Val(0)
       Return
    End If

    ! Compute Pm,m+1
    Val(1) = X*(2.0_SP*mabs+1.0_SP)*Val(0)
    If (l == mabs+1) Then
       Leg = NegFac * Val(1)
       Return
    End If

    ! For higher l, use the recursion relation
    Do I = mabs+2, l
       Aux = (2.0_SP*(I-1)+1.0_SP) * X * Val(1) - &
            & (I-1+mabs) * Val(0)
       Aux = Aux / Real(I-mabs,kind=SP)
       Val(0) = Val(1)
       Val(1) = Aux
    End Do
    Leg = NegFac * Val(1)

    Return
  End Function Legendre_SP

! *************************************
! *
  Complex (kind=DPC) Function SphericalHarmonic_DP(l, m, th, ph) Result (sph)
! *
! *************************************
! * Compute the value of the associated Legendre 
! * polinomial (l,m) at X.
! *************************************
    
    Real (kind=DP), Intent (in) :: th, ph
    Integer, Intent (in) :: l, m

    Real (kind=DP) :: Fac

    If (l == m) Then
       Fac = (2.0_DP*l + 1.0_DP)/(4.0_DP*PI_DP) / &
            & Factorial(l+m)
    Else If (l == -m) Then
       Fac = (2.0_DP*l + 1.0_DP)/(4.0_DP*PI_DP) * &
            & Factorial(l-m)
    Else
       Fac = (2.0_DP*l + 1.0_DP)/(4.0_DP*PI_DP) * &
            & Exp(GammaLn(Real(l-m,kind=DP)) - & 
            & GammaLn(Real(l+m,kind=DP)) )
    End If

    Sph = Cmplx(Sqrt(Fac) * Legendre(l,m,Cos(th)),kind=DPC) * &
         & Exp(UnitImag_DPC*m*ph)


    Return
  End Function SphericalHarmonic_DP


! *************************************
! *
  Complex (kind=SPC) Function SphericalHarmonic_SP(l, m, th, ph) Result (sph)
! *
! *************************************
! * Compute the value of the associated Legendre 
! * polinomial (l,m) at X.
! *************************************
    
    Real (kind=SP), Intent (in) :: th, ph
    Integer, Intent (in) :: l, m

    Real (kind=SP) :: Fac

    If (l == m) Then
       Fac = (2.0_SP*l + 1.0_SP)/(4.0_SP*PI_SP) / &
            & Real(Factorial(l+m), kind=SP)
    Else If (l == -m) Then
       Fac = (2.0_SP*l + 1.0_SP)/(4.0_SP*PI_SP) * &
            & Real(Factorial(l-m), kind=SP)
    Else
       Fac = (2.0_SP*l + 1.0_SP)/(4.0_SP*PI_SP) * &
            & Real(Exp(GammaLn(Real(l-m,kind=DP)) - & 
            & GammaLn(Real(l+m,kind=DP)) ), kind=SP)
    End If

    Sph = Cmplx(Sqrt(Fac) * Legendre(l,m,Cos(th)),kind=SPC) * &
         & Exp(UnitImag_SPC*m*ph)


    Return
  End Function SphericalHarmonic_SP

! *************************************
! *
  Real (kind=DP) Function inverf_DP(X)
! *
! *************************************
! * Compute the inverse error function for real
! * arguments. DP version.
! * Based on implementation by Peter John Acklam 
! * <pjacklam@online.no> supossed to be exact up 
! * to double precision.
! *************************************

    Real (kind=DP), Intent (in) :: X

    Real (kind=DP) :: Cfa(6), Cfb(5), Cfc(6), Cfd(4)

    Real (kind=DP) :: pl = 0.02425_DP, ph = 0.97575_DP, q, r, Xx

    ! Coeficients for first interval computation
    Data Cfa / -3.969683028665376e+01_DP, 2.209460984245205e+02_DP, &
         & -2.759285104469687e+02_DP, 1.383577518672690e+02_DP, &
         & -3.066479806614716e+01_DP, 2.506628277459239e+00_DP /
    Data Cfb / -5.447609879822406e+01_DP, 1.615858368580409e+02_DP, &
         & -1.556989798598866e+02_DP, 6.680131188771972e+01_DP, &
         & -1.328068155288572e+01_DP /

    ! Coeficients for second interval computation
    Data Cfc / -7.784894002430293e-03_DP, -3.223964580411365e-01_DP, &
         & -2.400758277161838e+00_DP, -2.549732539343734e+00_DP, &
         & 4.374664141464968e+00_DP, 2.938163982698783e+00_DP /
    Data Cfd / 7.784695709041462e-03_DP, 3.224671290700398e-01_DP, &
         & 2.445134137142996e+00_DP, 3.754408661907416e+00_DP/

    Xx = (X+1.0_DP)/2.0_DP

    If ((Xx > 0.0_DP).and.(Xx < pl)) Then
       q = Sqrt(-2.0_DP*Log(Xx))
       inverf_DP = &
            & (((((Cfc(1)*q + Cfc(2))*q + Cfc(3))*q + Cfc(4))*q + Cfc(5))*q + Cfc(6)) /&
            & ((((Cfd(1)*q + Cfd(2))*q + Cfd(3))*q + Cfd(4))*q + 1.0_DP)
    Else If ((Xx >= pl).and.(Xx <= ph)) Then
       q = Xx - 0.5_DP
       r = q*q
       inverf_DP = &
            & (((((Cfa(1)*r + Cfa(2))*r + Cfa(3))*r + Cfa(4))*r +&
            & Cfa(5))*r + Cfa(6))*q /&
            & (((((Cfb(1)*r + Cfb(2))*r + Cfb(3))*r + &
            & Cfb(4))*r + Cfb(5))*r + 1.0_DP)       
    Else If ((Xx > ph).and.(Xx < 1.0_DP)) Then
       q = Sqrt(-2.0_DP*Log(1.0_DP - Xx))
       inverf_DP = - &
            & (((((Cfc(1)*q + Cfc(2))*q + Cfc(3))*q + Cfc(4))*q + Cfc(5))*q + Cfc(6)) /&
            & ((((Cfd(1)*q + Cfd(2))*q + Cfd(3))*q + Cfd(4))*q + 1.0_DP)
    Else
       CALL Perror("inverf", "Argument out of domain")
       
    End If


    ! Refinement using Haley's rational method
    If ((Xx > 0.0_DP).and.(Xx < 1.0_DP) ) Then
       q = 0.5_DP * erfc_DP(-inverf_DP/SR2_DP) - Xx
       r = q * Sqrt(TWOPI_DP) * exp(inverf_DP**2/2.0_DP)
       inverf_DP = inverf_DP - r/(1.0_DP + r*inverf_DP/2.0_DP)
    End If
 

    inverf_DP = inverf_DP / SR2_DP
    
    Return
  End Function inverf_DP

! *************************************
! *
  Real (kind=SP) Function inverf_SP(X)
! *
! *************************************
! * Compute the inverse error function for real
! * arguments. SP version.
! * Based on implementation by Peter John Acklam 
! * <pjacklam@online.no> supossed to be exact up 
! * to double precision.
! *************************************

    Real (kind=SP), Intent (in) :: X

    Real (kind=SP) :: Cfa(6), Cfb(5), Cfc(6), Cfd(4)

    Real (kind=SP) :: pl = 0.02425_SP, ph = 0.97575_SP, q, r, Xx

    ! Coeficients for first interval computation
    Data Cfa / -3.969683028665376e+01_SP, 2.209460984245205e+02_SP, &
         & -2.759285104469687e+02_SP, 1.383577518672690e+02_SP, &
         & -3.066479806614716e+01_SP, 2.506628277459239e+00_SP /
    Data Cfb / -5.447609879822406e+01_SP, 1.615858368580409e+02_SP, &
         & -1.556989798598866e+02_SP, 6.680131188771972e+01_SP, &
         & -1.328068155288572e+01_SP /

    ! Coeficients for second interval computation
    Data Cfc / -7.784894002430293e-03_SP, -3.223964580411365e-01_SP, &
         & -2.400758277161838e+00_SP, -2.549732539343734e+00_SP, &
         & 4.374664141464968e+00_SP, 2.938163982698783e+00_SP /
    Data Cfd / 7.784695709041462e-03_SP, 3.224671290700398e-01_SP, &
         & 2.445134137142996e+00_SP, 3.754408661907416e+00_SP/

    Xx = (X+1.0_SP)/2.0_SP

    If ((Xx > 0.0_SP).and.(Xx < pl)) Then
       q = Sqrt(-2.0_SP*Log(Xx))
       inverf_SP = &
            & (((((Cfc(1)*q + Cfc(2))*q + Cfc(3))*q + Cfc(4))*q + Cfc(5))*q + Cfc(6)) /&
            & ((((Cfd(1)*q + Cfd(2))*q + Cfd(3))*q + Cfd(4))*q + 1.0_SP)
    Else If ((Xx >= pl).and.(Xx <= ph)) Then
       q = Xx - 0.5_SP
       r = q*q
       inverf_SP = &
            & (((((Cfa(1)*r + Cfa(2))*r + Cfa(3))*r + Cfa(4))*r +&
            & Cfa(5))*r + Cfa(6))*q /&
            & (((((Cfb(1)*r + Cfb(2))*r + Cfb(3))*r + &
            & Cfb(4))*r + Cfb(5))*r + 1.0_SP)       
    Else If ((Xx > ph).and.(Xx < 1.0_SP)) Then
       q = Sqrt(-2.0_SP*Log(1.0_SP - Xx))
       inverf_SP = - &
            & (((((Cfc(1)*q + Cfc(2))*q + Cfc(3))*q + Cfc(4))*q + Cfc(5))*q + Cfc(6)) /&
            & ((((Cfd(1)*q + Cfd(2))*q + Cfd(3))*q + Cfd(4))*q + 1.0_SP)
    Else
       CALL Perror("inverf", "Argument out of domain")
       
    End If


    ! Refinement using Haley's rational method
    If ((Xx > 0.0_SP).and.(Xx < 1.0_SP) ) Then
       q = 0.5_SP * erfc_SP(-inverf_SP/SR2_SP) - Xx
       r = q * Sqrt(TWOPI_SP) * exp(inverf_SP**2/2.0_SP)
       inverf_SP = inverf_SP - r/(1.0_SP + r*inverf_SP/2.0_SP)
    End If
 

    inverf_SP = inverf_SP / SR2_SP
    
    Return
  End Function inverf_SP

End MODULE SpecialFunc
