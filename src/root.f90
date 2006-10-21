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

  
  Real (kind=SP), Parameter :: DEFTOL_SP = 2.0E-2_SP
  Real (kind=DP), Parameter :: DEFTOL_DP = 2.0E-2_SP

  Interface RootPol
     Module Procedure RootPol2R_SP, RootPol2C_SP, RootPol3R_SP, RootPol3C_SP,&
          & RootPol4R_SP, RootPol2R_DP, RootPol2C_DP, RootPol3R_DP, &
          & RootPol3C_DP, RootPol4R_DP
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
       & Newton_SP, Newton_DP, Bisec_SP, Bisec_DP

CONTAINS

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

    z1 = ( -Cmplx(a, kind=SPC) + Sqrt(Cmplx(a**2-4*b, kind=SPC)) )/ 2.0_SP
    z2 = ( -Cmplx(a, kind=SPC) - Sqrt(Cmplx(a**2-4*b, kind=SPC)) )/ 2.0_SP
    

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

    z1 = ( -a + Sqrt(a**2-4*b) )/ 2.0_SP
    z2 = ( -a - Sqrt(a**2-4*b) )/ 2.0_SP
    
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

    Real (kind=SP) :: Q, R, Rint
    Complex (kind=SPC) :: S, T, Wint

    Q = ( 3*b - a**2  ) / 9.0_SP
    R = (9.0_SP*a*b-27.0_SP*c-2.0_SP*a**3)/54.0_SP

    Rint = Q**3 + R**2
    S = Cmplx( Sign((Abs(R + Rint))**(1./3.), R+Rint), kind=SPC)
    T = Cmplx( Sign((Abs(R - Rint))**(1./3.), R-Rint), kind=SPC)

    
    z1 = S + T - a/3.0_SP
    z2 = -(S+T)*.5 - Cmplx(a, kind=SPC)/3. + (0., 0.5)*Sqrt(3.0_SP)*(S-T)
    z3 = -(S+T)*.5 - Cmplx(a, kind=SPC)/3. - (0., 0.5)*Sqrt(3.0_SP)*(S-T)

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

    Complex (kind=SPC) :: Q, R, S, T, Wint

    Q = ( 3*b - a**2  ) / 9.0_SP
    R = (9.0_SP*a*b-27.0_SP*c-2.0_SP*a**3)/54.0_SP

    Wint = Sqrt(Q**3 + R**2)
    S = (R + Wint)**(1./3.)
    T = (R - Wint)**(1./3.)
    
    z1 = S + T - a/3.0_SP
    z2 = -(S+T)*.5 - a/3. + (0., 0.5)*Sqrt(3.0_SP)*(S-T)
    z3 = -(S+T)*.5 - a/3. - (0., 0.5)*Sqrt(3.0_SP)*(S-T)

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

    Complex (kind=SPC) :: Wint, Dumm1, Dumm2

    CALL RootPol(-b, a*c-4.*d, 4.*b*d-c**3-a**2*d, Wint, Dumm1, Dumm2)

    Dumm1 = Sqrt(Cmplx(a**2 - 4.*b + 4.*Wint, kind=SPC))
    Dumm2 = Sqrt(Cmplx(Wint**2 - 4.*d, kind=SPC))
    CALL RootPol( 0.5*(a+Dumm1), 0.5*(Wint-Dumm2), z1, z2 )
    CALL RootPol( 0.5*(a-Dumm1), 0.5*(Wint+Dumm2), z3, z4 )

    Return
  End Subroutine RootPol4R_SP


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

    z1 = ( -Cmplx(a, kind=DPC) + Sqrt(Cmplx(a**2-4*b, kind=DPC)) )/ 2.0_DP
    z2 = ( -Cmplx(a, kind=DPC) - Sqrt(Cmplx(a**2-4*b, kind=DPC)) )/ 2.0_DP
    

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

    z1 = ( -a + Sqrt(a**2-4*b) )/ 2.0_DP
    z2 = ( -a - Sqrt(a**2-4*b) )/ 2.0_DP
    
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

    Real (kind=DP) :: Q, R, Rint
    Complex (kind=DPC) :: S, T, Wint

    Q = ( 3*b - a**2  ) / 9.0_DP
    R = (9.0_DP*a*b-27.0_DP*c-2.0_DP*a**3)/54.0_DP

    Rint = Q**3 + R**2
    S = Cmplx( Sign((Abs(R + Rint))**(1./3.), R+Rint), kind=DPC)
    T = Cmplx( Sign((Abs(R - Rint))**(1./3.), R-Rint), kind=DPC)

    
    z1 = S + T - a/3.0_DP
    z2 = -(S+T)*.5 - Cmplx(a, kind=DPC)/3. + (0., 0.5)*Sqrt(3.0_DP)*(S-T)
    z3 = -(S+T)*.5 - Cmplx(a, kind=DPC)/3. - (0., 0.5)*Sqrt(3.0_DP)*(S-T)

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

    Complex (kind=DPC) :: Q, R, S, T, Wint

    Q = ( 3*b - a**2  ) / 9.0_DP
    R = (9.0_DP*a*b-27.0_DP*c-2.0_DP*a**3)/54.0_DP

    Wint = Sqrt(Q**3 + R**2)
    S = (R + Wint)**(1./3.)
    T = (R - Wint)**(1./3.)
    
    z1 = S + T - a/3.0_DP
    z2 = -(S+T)*.5 - a/3. + (0., 0.5)*Sqrt(3.0_DP)*(S-T)
    z3 = -(S+T)*.5 - a/3. - (0., 0.5)*Sqrt(3.0_DP)*(S-T)

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

    Complex (kind=DPC) :: Wint, Dumm1, Dumm2

    CALL RootPol(-b, a*c-4.*d, 4.*b*d-c**3-a**2*d, Wint, Dumm1, Dumm2)

    Dumm1 = Sqrt(Cmplx(a**2 - 4.*b + 4.*Wint, kind=DPC))
    Dumm2 = Sqrt(Cmplx(Wint**2 - 4.*d, kind=DPC))
    CALL RootPol( 0.5*(a+Dumm1), 0.5*(Wint-Dumm2), z1, z2 )
    CALL RootPol( 0.5*(a-Dumm1), 0.5*(Wint+Dumm2), z3, z4 )

    Return
  End Subroutine RootPol4R_DP


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


End MODULE Root

