!
! MODULE with routine to calculate Roots of equations.
!
! Copyright (C) 2005  Alberto Ramos <alberto@martin.ft.uam.es>
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

! $ v. 1.0; Released: 15/09/2005; $

! ***************************************************
! *
MODULE Root
! *
! ***************************************************

  USE NumTypes
  USE Error
  USE Constants
  
  Real (kind=SP), Parameter :: DEFTOL_SP = 2.0E-2_SP
  Real (kind=DP), Parameter :: DEFTOL_DP = 2.0E-2_DP

  Interface RootPol
     Module Procedure RootPol2R_SP, RootPol2C_SP, RootPol3R_SP, RootPol3C_SP,&
          & RootPol4R_SP, RootPol2R_DP, RootPol2C_DP, RootPol3R_DP, &
          & RootPol3C_DP, RootPol4R_DP, RootPol4C_DP, RootPol4C_SP
  End Interface
  
  Interface Newton
     Module Procedure Newton_SP, Newton_DP
  End Interface

  Interface Bisec
     Module Procedure Bisec_SP, Bisec_DP
  End Interface

  Private RootPol2R_SP, RootPol2C_SP, RootPol3R_SP, RootPol3C_SP, &
       & RootPol4R_SP, DEFTOL_SP, DEFTOL_DP, RootPol2R_DP, &
       & RootPol2C_DP, RootPol3R_DP, RootPol3C_DP, RootPol4R_DP, &
       & Newton_SP, Newton_DP, Bisec_SP, Bisec_DP, RootPol4C_DP, &
       & RootPol4C_SP

CONTAINS

! **********************************************
! *                                            *
  Subroutine RootPol2R_DP(a, b, z1, z2) 
! *                                            *
! **********************************************
! * Calculates the Roots of a 2 degree 
! * polinomial. Real DP Version.
! **********************************************

    Real (kind=DP), Intent (in) :: a, b
    Complex (kind=DPC), Intent (out) :: z1, z2
    Complex (kind=DPC) :: Q

    Q = Sqrt(Cmplx(a**2 - 4.0_DP*b,kind=DPC))
    If (a > 0.0_DP) Then
       Q = -0.5_DP * (Cmplx(a,kind=DPC) + Q)
    Else
       Q = -0.5_DP * (Cmplx(a,kind=DPC) - Q)
    End If
    
    z1 = Q
    z2 = Cmplx(b,kind=DPC)/Q

    Return
  End Subroutine RootPol2R_DP

! **********************************************
! *                                            *
  Subroutine RootPol2C_DP(a, b, z1, z2) 
! *                                            *
! **********************************************
! * Calculates the Roots of a 2 degree 
! * polinomial. Complex DP Version.
! **********************************************

    Complex (kind=DPC), Intent (in) :: a, b
    Complex (kind=DPC), Intent (out) :: z1, z2
    Complex (kind=DPC) :: Q
    
    Q = Sqrt(a**2 - 4.0_DP*b)
    If (Real(Conjg(a)*Q,kind=DP) > 0.0_DP) Then
       Q = -0.5_DP*(a + Q)
    Else
       Q = -0.5_DP*(a - Q)
    End If
    
    z1 = Q
    z2 = b/Q

    Return
  End Subroutine RootPol2C_DP

! **********************************************
! *                                            *
  Subroutine RootPol3R_DP(a, b, c, z1, z2, z3)
! *                                            *
! **********************************************
! * Calculates the Roots of a 3 degree 
! * polinomial. Real DP Version.
! **********************************************

    Real (kind=DP), Intent (in) :: a, b, c
    Complex (kind=DPC), Intent (out) :: z1, z2, z3

    Real (kind=DP) :: Q, R, T
    Complex (kind=DPC) :: S, S2

    Q = ( a**2 - 3.0_DP*b  ) / 9.0_DP
    R = (-9.0_DP*a*b + 27.0_DP*c + 2.0_DP*a**3)/54.0_DP

    If (R**2 < Q**3) Then
       T = Acos(R/Sqrt(Q**3))
       
       Q = -2.0_DP*Sqrt(Q)
       z1 = Q*Cos(T/3.0_DP) - a/3.0_DP
       z2 = Q*Cos((T+TWOPI_DP)/3.0_DP) - a/3.0_DP
       z3 = Q*Cos((T-TWOPI_DP)/3.0_DP) - a/3.0_DP
    Else 
       S = Sqrt(Cmplx(R**2-Q**3,kind=DPC))
       If (Real(R*S) > 0) Then
          S = -(Cmplx(R,kind=DPC) + S)**(1.0_DP/3.0_DP)
       Else
          S = -(Cmplx(R,kind=DPC) - S)**(1.0_DP/3.0_DP)
       End If
       If (Abs(S) == 0.0_DP) Then
          S2 = 0.0_DP
       Else
          S2 = Cmplx(Q,kind=DPC)/S
       End If

       z1 = S + S2 - a/3.0_DP
       z2 = -0.5_DP*(S+S2) - a/3.0_DP + &
            & Cmplx(0.0_DP,SR3_DP/2.0_DP,kind=DPC) * (S-S2)
       z3 = -0.5_DP*(S+S2) - a/3.0_DP - &
            & Cmplx(0.0_DP,SR3_DP/2.0_DP,kind=DPC) * (S-S2)
    End If

    Return
  End Subroutine RootPol3R_DP

! **********************************************
! *                                            *
  Subroutine RootPol3C_DP(a, b, c, z1, z2, z3)
! *                                            *
! **********************************************
! * Calculates the Roots of a 3 degree 
! * polinomial. Complex DP Version.
! **********************************************

    Complex (kind=DPC), Intent (in) :: a, b, c
    Complex (kind=DPC), Intent (out) :: z1, z2, z3

    Complex (kind=DPC) :: S, S2, Q, R

    Q = ( a**2 - 3.0_DP*b  ) / 9.0_DP
    R = (-9.0_DP*a*b + 27.0_DP*c + 2.0_DP*a**3)/54.0_DP

    S = Sqrt(Cmplx(R**2-Q**3,kind=DPC))
    If (Real(Conjg(R)*S,kind=DP) > 0.0_DP) Then
       S = -(Cmplx(R,kind=DPC) + S)**(1.0_DP/3.0_DP)
    Else
       S = -(Cmplx(R,kind=DPC) - S)**(1.0_DP/3.0_DP)
    End If
    If (Abs(S) == 0.0_DP) Then
       S2 = 0.0_DP
    Else
       S2 = Q/S
    End If
    
    z1 = S + S2 - a/3.0_DP
    z2 = -0.5_DP*(S+S2) - a/3.0_DP + &
         & Cmplx(0.0_DP,SR3_DP/2.0_DP,kind=DPC) * (S-S2)
    z3 = -0.5_DP*(S+S2) - a/3.0_DP - &
         & Cmplx(0.0_DP,SR3_DP/2.0_DP,kind=DPC) * (S-S2)


    Return
  End Subroutine RootPol3C_DP

! **********************************************
! *                                            *
  Subroutine RootPol4R_DP(a, b, c, d, z1, z2, z3, z4)
! *                                            *
! **********************************************
! * Calculates the Roots of a 4 degree 
! * polinomial. Real DP Version.
! **********************************************

    Real (kind=DP), Intent (in) :: a, b, c, d
    Complex (kind=DPC), Intent (out) :: z1, z2, z3, z4

    Real (kind=DP) :: e, f, g
    Complex (kind=DPC) :: Dumm3, Dumm1, Dumm2, p, q, r

    e = b - 3.0_DP*a**2/8.0_DP
    f = c + a**3/8.0_DP - a*b/2.0_DP
    g = d - 3.0_DP * a**4/256.0_DP + a**2 * b/16.0_DP - a*c/4.0_DP
    CALL RootPol(e/2.0_DP, (e**2-4.0_DP*g)/16.0_DP , -f**2/64.0_DP, &
         & Dumm1, Dumm2, Dumm3)

    p = Sqrt(Dumm1)
    q = Sqrt(Dumm2)
    r = -f/(8.00_DP*p*q)

    z1 =  p + q + r - Cmplx(a,kind=DPC)/4.0_DP
    z2 =  p - q - r - Cmplx(a,kind=DPC)/4.0_DP
    z3 = -p + q - r - Cmplx(a,kind=DPC)/4.0_DP
    z4 = -p - q + r - Cmplx(a,kind=DPC)/4.0_DP

    Return
  End Subroutine RootPol4R_DP

! **********************************************
! *                                            *
  Subroutine RootPol4C_DP(a, b, c, d, z1, z2, z3, z4)
! *                                            *
! **********************************************
! * Calculates the Roots of a 4 degree 
! * polinomial. Cmplex DP Version.
! **********************************************

    Complex (kind=DPC), Intent (in) :: a, b, c, d
    Complex (kind=DPC), Intent (out) :: z1, z2, z3, z4

    Complex (kind=DPC) :: e, f, g
    Complex (kind=DPC) :: Dumm3, Dumm1, Dumm2, p, q, r

    e = b - 3.0_DP*a**2/8.0_DP
    f = c + a**3/8.0_DP - a*b/2.0_DP
    g = d - 3.0_DP * a**4/256.0_DP + a**2 * b/16.0_DP - a*c/4.0_DP
    CALL RootPol(e/2.0_DP, (e**2-4.0_DP*g)/16.0_DP , -f**2/64.0_DP, &
         & Dumm1, Dumm2, Dumm3)

    p = Sqrt(Dumm1)
    q = Sqrt(Dumm2)
    r = -f/(8.00_DP*p*q)

    z1 =  p + q + r - Cmplx(a,kind=DPC)/4.0_DP
    z2 =  p - q - r - Cmplx(a,kind=DPC)/4.0_DP
    z3 = -p + q - r - Cmplx(a,kind=DPC)/4.0_DP
    z4 = -p - q + r - Cmplx(a,kind=DPC)/4.0_DP

    Return
  End Subroutine RootPol4C_DP


! **********************************************
! *                                            *
  Real (kind=DP) Function Newton_DP(Xo, Fnew, Tolerance) 
! *                                            *
! **********************************************
! * Calculates the root of the Function 
! * Fnew(X, F, D), that returns the value of the 
! * function (F) and its Derivative (D). 
! * DP Version.
! **********************************************

    Real (kind=DP), Intent (in) :: Xo
    Real (kind=DP), Intent (in), Optional :: Tolerance

    Real (kind=DP) :: Xnew, Xant, D, Func, Tol

    Interface
       Subroutine FNew(Xo, F, D)

         USE NumTypes

         Real (kind=DP), Intent (in) :: Xo
         Real (kind=DP), Intent (out) :: F, D
       End Subroutine FNew
    End Interface


    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_DP
    Else
       Tol = Tolerance
    End If

    Xant = Xo + 2.0_DP
    Xnew = Xo 
    Do While (Abs(Xnew-Xant).gt.Tol)
       Xant = Xnew
       CALL Fnew(Xant, Func, D)
       Xnew = Xant - Func/D
    End Do
    
    Newton_DP = Xnew
    
    Return
  End Function Newton_DP


! **********************************************
! *                                            *
  Real (kind=DP) Function Bisec_DP(aout, bout, F, Tolerance)
! *                                            *
! **********************************************
! * Calculates the root of the Function 
! * F(X), provided that you give a pair of points
! * a, b such that F(a)*F(b) < 0.
! * DP Version.
! **********************************************

    Real (kind=DP), Intent (in) :: aout, bout
    Real (kind=DP), Intent (in), Optional :: Tolerance

    Real (kind=DP) :: Err, a, b

    Interface
       Function F(X)

         USE NumTypes

         Real (kind=DP), Intent (in) :: X
         Real (kind=DP) :: F
       End Function F
    End Interface

    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_DP
    Else
       Tol = Tolerance
    End If

    If (F(aout)*F(bout) .ge. 0.0_DP) &
         & CALL Abort('Bisec_DP','Bad initial segment')

    a = aout
    b = bout
    Err = Tol + 1.0_DP
    Do While (Err > Tol)
       Bisec_DP = (a+b)/2.0_DP
       If (Abs(F(Bisec_DP)) .eq. 0.0_DP) Return
       If (F(Bisec_DP)*F(b) < 0.0_DP) Then
          Err = Abs(Bisec_DP - a)
          a = Bisec_DP
       Else
          Err = Abs(Bisec_DP - b)
          b = Bisec_DP
       End If
    End Do

    
    Return
  End Function Bisec_DP

! **********************************************
! *                                            *
  Subroutine RootPol2R_SP(a, b, z1, z2) 
! *                                            *
! **********************************************
! * Calculates the Roots of a 2 degree 
! * polinomial. Real SP Version.
! **********************************************

    Real (kind=SP), Intent (in) :: a, b
    Complex (kind=SPC), Intent (out) :: z1, z2
    Complex (kind=SPC) :: Q

    Q = Sqrt(Cmplx(a**2 - 4.0_SP*b,kind=SPC))
    If (a > 0.0_SP) Then
       Q = -0.5_SP * (Cmplx(a,kind=SPC) + Q)
    Else
       Q = -0.5_SP * (Cmplx(a,kind=SPC) - Q)
    End If
    
    z1 = Q
    z2 = Cmplx(b,kind=SPC)/Q

    Return
  End Subroutine RootPol2R_SP

! **********************************************
! *                                            *
  Subroutine RootPol2C_SP(a, b, z1, z2) 
! *                                            *
! **********************************************
! * Calculates the Roots of a 2 degree 
! * polinomial. Complex SP Version.
! **********************************************

    Complex (kind=SPC), Intent (in) :: a, b
    Complex (kind=SPC), Intent (out) :: z1, z2
    Complex (kind=SPC) :: Q
    
    Q = Sqrt(a**2 - 4.0_SP*b)
    If (Real(Conjg(a)*Q,kind=SP) > 0.0_SP) Then
       Q = -0.5_SP*(a + Q)
    Else
       Q = -0.5_SP*(a - Q)
    End If
    
    z1 = Q
    z2 = b/Q

    Return
  End Subroutine RootPol2C_SP

! **********************************************
! *                                            *
  Subroutine RootPol3R_SP(a, b, c, z1, z2, z3)
! *                                            *
! **********************************************
! * Calculates the Roots of a 3 degree 
! * polinomial. Real SP Version.
! **********************************************

    Real (kind=SP), Intent (in) :: a, b, c
    Complex (kind=SPC), Intent (out) :: z1, z2, z3

    Real (kind=SP) :: Q, R, T
    Complex (kind=SPC) :: S, S2

    Q = ( a**2 - 3.0_SP*b  ) / 9.0_SP
    R = (-9.0_SP*a*b + 27.0_SP*c + 2.0_SP*a**3)/54.0_SP

    If (R**2 < Q**3) Then
       T = Acos(R/Sqrt(Q**3))
       
       Q = -2.0_SP*Sqrt(Q)
       z1 = Q*Cos(T/3.0_SP) - a/3.0_SP
       z2 = Q*Cos((T+TWOPI_SP)/3.0_SP) - a/3.0_SP
       z3 = Q*Cos((T-TWOPI_SP)/3.0_SP) - a/3.0_SP
    Else 
       S = Sqrt(Cmplx(R**2-Q**3,kind=SPC))
       If (Real(R*S) > 0) Then
          S = -(Cmplx(R,kind=SPC) + S)**(1.0_SP/3.0_SP)
       Else
          S = -(Cmplx(R,kind=SPC) - S)**(1.0_SP/3.0_SP)
       End If
       If (Abs(S) == 0.0_SP) Then
          S2 = 0.0_SP
       Else
          S2 = Cmplx(Q,kind=SPC)/S
       End If

       z1 = S + S2 - a/3.0_SP
       z2 = -0.5_SP*(S+S2) - a/3.0_SP + &
            & Cmplx(0.0_SP,SR3_SP/2.0_SP,kind=SPC) * (S-S2)
       z3 = -0.5_SP*(S+S2) - a/3.0_SP - &
            & Cmplx(0.0_SP,SR3_SP/2.0_SP,kind=SPC) * (S-S2)
    End If

    Return
  End Subroutine RootPol3R_SP

! **********************************************
! *                                            *
  Subroutine RootPol3C_SP(a, b, c, z1, z2, z3)
! *                                            *
! **********************************************
! * Calculates the Roots of a 3 degree 
! * polinomial. Complex SP Version.
! **********************************************

    Complex (kind=SPC), Intent (in) :: a, b, c
    Complex (kind=SPC), Intent (out) :: z1, z2, z3

    Complex (kind=SPC) :: S, S2, Q, R

    Q = ( a**2 - 3.0_SP*b  ) / 9.0_SP
    R = (-9.0_SP*a*b + 27.0_SP*c + 2.0_SP*a**3)/54.0_SP

    S = Sqrt(Cmplx(R**2-Q**3,kind=SPC))
    If (Real(Conjg(R)*S,kind=SP) > 0.0_SP) Then
       S = -(Cmplx(R,kind=SPC) + S)**(1.0_SP/3.0_SP)
    Else
       S = -(Cmplx(R,kind=SPC) - S)**(1.0_SP/3.0_SP)
    End If
    If (Abs(S) == 0.0_SP) Then
       S2 = 0.0_SP
    Else
       S2 = Q/S
    End If
    
    z1 = S + S2 - a/3.0_SP
    z2 = -0.5_SP*(S+S2) - a/3.0_SP + &
         & Cmplx(0.0_SP,SR3_SP/2.0_SP,kind=SPC) * (S-S2)
    z3 = -0.5_SP*(S+S2) - a/3.0_SP - &
         & Cmplx(0.0_SP,SR3_SP/2.0_SP,kind=SPC) * (S-S2)


    Return
  End Subroutine RootPol3C_SP

! **********************************************
! *                                            *
  Subroutine RootPol4R_SP(a, b, c, d, z1, z2, z3, z4)
! *                                            *
! **********************************************
! * Calculates the Roots of a 4 degree 
! * polinomial. Real SP Version.
! **********************************************

    Real (kind=SP), Intent (in) :: a, b, c, d
    Complex (kind=SPC), Intent (out) :: z1, z2, z3, z4

    Real (kind=SP) :: e, f, g
    Complex (kind=SPC) :: Dumm3, Dumm1, Dumm2, p, q, r

    e = b - 3.0_SP*a**2/8.0_SP
    f = c + a**3/8.0_SP - a*b/2.0_SP
    g = d - 3.0_SP * a**4/256.0_SP + a**2 * b/16.0_SP - a*c/4.0_SP
    CALL RootPol(e/2.0_SP, (e**2-4.0_SP*g)/16.0_SP , -f**2/64.0_SP, &
         & Dumm1, Dumm2, Dumm3)

    p = Sqrt(Dumm1)
    q = Sqrt(Dumm2)
    r = -f/(8.00_SP*p*q)

    z1 =  p + q + r - Cmplx(a,kind=SPC)/4.0_SP
    z2 =  p - q - r - Cmplx(a,kind=SPC)/4.0_SP
    z3 = -p + q - r - Cmplx(a,kind=SPC)/4.0_SP
    z4 = -p - q + r - Cmplx(a,kind=SPC)/4.0_SP

    Return
  End Subroutine RootPol4R_SP

! **********************************************
! *                                            *
  Subroutine RootPol4C_SP(a, b, c, d, z1, z2, z3, z4)
! *                                            *
! **********************************************
! * Calculates the Roots of a 4 degree 
! * polinomial. Cmplex SP Version.
! **********************************************

    Complex (kind=SPC), Intent (in) :: a, b, c, d
    Complex (kind=SPC), Intent (out) :: z1, z2, z3, z4

    Complex (kind=SPC) :: e, f, g
    Complex (kind=SPC) :: Dumm3, Dumm1, Dumm2, p, q, r

    e = b - 3.0_SP*a**2/8.0_SP
    f = c + a**3/8.0_SP - a*b/2.0_SP
    g = d - 3.0_SP * a**4/256.0_SP + a**2 * b/16.0_SP - a*c/4.0_SP
    CALL RootPol(e/2.0_SP, (e**2-4.0_SP*g)/16.0_SP , -f**2/64.0_SP, &
         & Dumm1, Dumm2, Dumm3)

    p = Sqrt(Dumm1)
    q = Sqrt(Dumm2)
    r = -f/(8.00_SP*p*q)

    z1 =  p + q + r - Cmplx(a,kind=SPC)/4.0_SP
    z2 =  p - q - r - Cmplx(a,kind=SPC)/4.0_SP
    z3 = -p + q - r - Cmplx(a,kind=SPC)/4.0_SP
    z4 = -p - q + r - Cmplx(a,kind=SPC)/4.0_SP

    Return
  End Subroutine RootPol4C_SP


! **********************************************
! *                                            *
  Real (kind=SP) Function Newton_SP(Xo, Fnew, Tolerance) 
! *                                            *
! **********************************************
! * Calculates the root of the Function 
! * Fnew(X, F, D), that returns the value of the 
! * function (F) and its Derivative (D). 
! * SP Version.
! **********************************************

    Real (kind=SP), Intent (in) :: Xo
    Real (kind=SP), Intent (in), Optional :: Tolerance

    Real (kind=SP) :: Xnew, Xant, D, Func, Tol

    Interface
       Subroutine FNew(Xo, F, D)

         USE NumTypes

         Real (kind=SP), Intent (in) :: Xo
         Real (kind=SP), Intent (out) :: F, D
       End Subroutine FNew
    End Interface


    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_SP
    Else
       Tol = Tolerance
    End If

    Xant = Xo + 2.0_SP
    Xnew = Xo 
    Do While (Abs(Xnew-Xant).gt.Tol)
       Xant = Xnew
       CALL Fnew(Xant, Func, D)
       Xnew = Xant - Func/D
    End Do
    
    Newton_SP = Xnew
    
    Return
  End Function Newton_SP


! **********************************************
! *                                            *
  Real (kind=SP) Function Bisec_SP(aout, bout, F, Tolerance)
! *                                            *
! **********************************************
! * Calculates the root of the Function 
! * F(X), provided that you give a pair of points
! * a, b such that F(a)*F(b) < 0.
! * SP Version.
! **********************************************

    Real (kind=SP), Intent (in) :: aout, bout
    Real (kind=SP), Intent (in), Optional :: Tolerance

    Real (kind=SP) :: Err, a, b

    Interface
       Function F(X)

         USE NumTypes

         Real (kind=SP), Intent (in) :: X
         Real (kind=SP) :: F
       End Function F
    End Interface

    If (.not. Present(Tolerance)) Then
       Tol = DEFTOL_SP
    Else
       Tol = Tolerance
    End If

    If (F(aout)*F(bout) .ge. 0.0_SP) &
         & CALL Abort('Bisec_SP','Bad initial segment')

    a = aout
    b = bout
    Err = Tol + 1.0_SP
    Do While (Err > Tol)
       Bisec_SP = (a+b)/2.0_SP
       If (Abs(F(Bisec_SP)) .eq. 0.0_SP) Return
       If (F(Bisec_SP)*F(b) < 0.0_SP) Then
          Err = Abs(Bisec_SP - a)
          a = Bisec_SP
       Else
          Err = Abs(Bisec_SP - b)
          b = Bisec_SP
       End If
    End Do

    
    Return
  End Function Bisec_SP

End MODULE Root

