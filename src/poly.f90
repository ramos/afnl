!
! A MODULE for manipulation and definition of polynomials
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

! $ v. 1.0; Released: 21/09/2005; $
! $ v. 2.0; Released: 26/12/2007; $

! ***********************************************************
! *
MODULE Polynomial
! *
! ***********************************************************

  USE NumTypes
  USE Error

  IMPLICIT NONE

  Type Pol
     Real (kind=DP), Allocatable :: Coef(:)
     Integer :: dg
  End Type Pol

  Type CmplxPol
     Complex (kind=DPC), Allocatable :: Coef(:)
     Integer :: dg
  End Type CmplxPol


  Interface Assignment (=)
     Module Procedure IgualM, IgualP, IgualCM, IgualCP
  End Interface

  Interface Operator (+)
     Module Procedure SumaP, SumaCP
  End Interface

  Interface Operator (-)
     Module Procedure RestaP, RestaCP
  End Interface

  Interface Operator (*)
     Module Procedure MultP, MultCl, MultCr, MultCP, MultCCl, &
          & MultCCr, MultCCCl, MultCCCr
  End Interface


  Interface Init
     Module Procedure Init, InitC
  End Interface

  Interface Degree
     Module Procedure Degree, DegreeC
  End Interface

  Interface Value
     Module Procedure Value, ValueC
  End Interface

  Interface Deriv
     Module Procedure Deriv, DerivC
  End Interface

  Interface Integra
     Module Procedure Integra, IntegraC
  End Interface

  Interface InterPolValue
     Module Procedure InterPolValue, InterPolValueC
  End Interface

  Interface InterPol
     Module Procedure InterPol, InterPolC
  End Interface


  Private IgualM, IgualP, SumaP, RestaP, MultP, MultCl, MultCr, &
       & IgualCM, IgualCP, SumaCP, RestaCP, MultCP, MultCCl, &
       & MultCCr, MultCCCr, MultCCCl, DerivC, IntegraC, InitC, &
       & DegreeC, ValueC, InterPolValueC, InterPolC


CONTAINS

!  *********************************************
!  *                                           *
  Subroutine Init(P, Ngrad)
!  *                                           *
!  *********************************************
!  * Init a polinomial of degree Ngrad
!  *********************************************

    Type (Pol), Intent (inout) :: P
    Integer, Intent (in) :: Ngrad


    ALLOCATE(P%Coef(0:Ngrad))

    P%Coef = 0.D0
    P%dg = Ngrad

    Return
  End Subroutine Init

!  *********************************************
!  *                                           *
  Integer Function Degree(P)
!  *                                           *
!  *********************************************
!  * Returns the degree of a ploy.
!  *********************************************

    Type (pol) :: P

    Degree = P%dg

    Return
  End Function Degree

!  *********************************************
!  *                                           *
  Subroutine IgualP(Pol1, Pol2)
!  *                                           *
!  *********************************************
!  * Equate two polinomials
!  *********************************************

    Type (Pol), Intent(out) :: Pol1
    Type (Pol), Intent(in) :: Pol2

    CALL Init(Pol1, Pol2%dg)
    Pol1%Coef = Pol2%Coef
    
    Return
  End Subroutine IgualP

!  *********************************************
!  *                                           *
  Subroutine IgualM(P, Coef)
!  *                                           *
!  *********************************************
!  * Iguala un polinomio a un Vector de 
!  * coeficientes
!  *********************************************

    Type (Pol), Intent (out) :: P
    Real (kind=DP), Intent (in) :: Coef(:)

    CALL Init(P,Size(Coef)-1)
    P%Coef = Coef
    
    Return
  End Subroutine IgualM

!  *********************************************
!  *                                           *
  Type (Pol) Function SumaP(Pol1, Pol2) Result (Su)
!  *                                           *
!  *********************************************
!  * Suma Dos polinomios
!  * 
!  *********************************************

    Type (Pol), Intent (in) :: Pol1, Pol2

    Integer :: I1, I2

    I1 = Pol1%dg
    I2 = Pol2%dg
    CALL Init(Su, Max(I1,I2))

    If (I1 > I2) Then
       Su%Coef = Pol1%Coef
       Su%Coef(0:I2) = Pol1%Coef(0:I2) + Pol2%Coef(0:I2)
    Else
       Su%Coef = Pol2%Coef
       Su%Coef(0:I1) = Pol1%Coef(0:I1) + Pol2%Coef(0:I1)
    End If


    Return
  End Function SumaP

!  *********************************************
!  *                                           *
  Type (Pol) Function RestaP(Pol1, Pol2) Result (Su)
!  *                                           *
!  *********************************************
!  * Resta Dos polinomios
!  * 
!  *********************************************

    Type (Pol), Intent (in) :: Pol1, Pol2

    Integer :: I1, I2

    I1 = Pol1%dg
    I2 = Pol2%dg
    CALL Init(Su, Max(I1,I2))

    If (I1 > I2) Then
       Su%Coef = Pol1%Coef
       Su%Coef(0:I2) = Pol1%Coef(0:I2) - Pol2%Coef(0:I2)
    Else
       Su%Coef = -Pol2%Coef
       Su%Coef(0:I1) = Pol1%Coef(0:I1) - Pol2%Coef(0:I1)
    End If

    Return
  End Function RestaP

!  *********************************************
!  *                                           *
  Type (Pol) Function MultP(Pol1, Pol2)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (Pol), Intent (in) :: Pol1, Pol2

    Integer :: I, J, Igrado1, Igrado2

    Igrado1 = Pol1%dg
    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado1 + Igrado2)


    MultP%Coef = 0.0_DP
    Do I = 0, Igrado1
       Do J = 0, Igrado2
          MultP%Coef(I+J) = MultP%Coef(I+J) + Pol1%Coef(I)*Pol2%Coef(J)
       End Do
    End Do

    Return
  End Function MultP

!  *********************************************
!  *                                           *
  Type (Pol) Function MultCl(C, Pol2) Result (MultP)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (Pol), Intent (in) :: Pol2
    Real (kind=DP), Intent (in) :: C

    Integer :: Igrado2

    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*Pol2%Coef

    Return
  End Function MultCl

!  *********************************************
!  *                                           *
  Type (Pol) Function MultCr(Pol2, C) Result (MultP)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (Pol), Intent (in) :: Pol2
    Real (kind=DP), Intent (in)  :: C

    Integer :: Igrado2

    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*Pol2%Coef


    Return
  End Function MultCr

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Value(Polin, X)
!  *                                           *
!  *********************************************
!  * Evaluates a Poly in X.
!  *********************************************
    
    Type (Pol), Intent (in) :: Polin
    Real (kind=DP), Intent (in) :: X

    Integer :: I

    Value = Polin%Coef(Polin%dg)*X
    Do I = Polin%dg-1, 1, -1
       Value = (Value + Polin%Coef(I))*X
    End Do
    Value = Value + Polin%Coef(0)

    Return
  End Function Value

!  *********************************************
!  *                                           *
  Type (Pol) Function Deriv(P)
!  *                                           *
!  *********************************************
!  * Calculates de derivative of a poly.
!  *********************************************

    Type (Pol), Intent (in) :: P

    Integer :: I

    CALL Init(Deriv, P%dg-1)
    
    Do I = 0, P%dg - 1
       Deriv%Coef(I) = Real(I+1, kind=DP)*P%Coef(I+1)
    End Do
    
    Return
  End Function Deriv

!  *********************************************
!  *                                           *
  Type (Pol) Function Integra(P, C)
!  *                                           *
!  *********************************************
!  * Integra un polinomio
!  * 
!  *********************************************

    Type (Pol), Intent (in) :: P
    Real (kind=DP), Optional :: C

    Integer :: I

    CALL Init(Integra, P%dg + 1)

    If (Present(C)) Then
       Integra%Coef(0) = C
    Else
       Integra%Coef(0) = 0.D0
    End If
       
    Do I = 1, Integra%dg
       Integra%Coef(I) = P%Coef(I-1)/Real(I,kind=DP)
    End Do
    
    Return
  End Function Integra

!  *********************************************
!  *                                           *
  Real (kind=DP) Function InterpolValue(X, Y, Xo)
!  *                                           *
!  *********************************************
!  * Calculates the value of the Interpolation 
!  * polynomial of data X(:), Y(:) in the point 
!  * Xo.
!  *********************************************

    Real (kind=DP), Intent (in) :: X(:), Y(:), Xo

    Real (kind=DP) :: a(Size(X),Size(X))
    Integer :: I, J, Npoints

    Npoints = Size(X)
    a = 0.0_DP 
    a(:,1) = Y(:)
    Do J = 2, Npoints
       Do I = J, Npoints
          a(I,J) = ( (Xo-X(I-J+1))*a(I,J-1) + (X(I)-Xo)*a(I-1,J-1) ) &
               & / (X(I) - X(I-J+1))
       End Do
    End Do
    
    InterpolValue = a(Npoints, Npoints)

    Return
  End Function InterpolValue

!  *********************************************
!  *                                           *
  Type (Pol) Function Interpol(X, Y)
!  *                                           *
!  *********************************************
!  * Calculates the Interpolation polynomial of
!  * data X(:), Y(:) using divided diference and
!  * the Newton form of the Interpolating 
!  * polynomial.
!  * 
!  * NOTE:
!  * In general this problem has a lot of 
!  * numerical errors, so it is much more 
!  * interesting to calculate the value of the 
!  * interpolating polynomial using the function 
!  * InterpolValue of this Module. 
!  *********************************************

    Real (kind=DP), Intent (in) :: X(:), Y(:)
    
    Type (Pol) :: Base, Mul
    Real (kind=DP) :: a(Size(X), Size(X))
    Integer :: I, J

    CALL Init(Interpol, Size(X)-1)
    CALL Init(Base,0)
    CALL Init(Mul,1)

    ! First we construct the Divided diferences
    a = 0.0_DP 
    a(:,1) = Y(:)
    Do J = 2, Interpol%dg + 1
       Do I = J, Interpol%dg + 1
          a(I,J) = ( a(I,J-1) - a(I-1,J-1) ) / (X(I) - X(I-J+1))
       End Do
    End Do

    ! Now construc the Polynomial in the Newton form.
    Interpol%Coef(0) = Y(1)
    Base%Coef(0) = 1.0_DP
    Do I = 1, Interpol%dg
       Mul%Coef(1) = 1.0_DP
       Mul%Coef(0) = -X(I)
       Base = Base * Mul

       Interpol = Interpol + a(I+1,I+1)*Base
    End Do

    Return
  End Function Interpol

!  *********************************************
!  *                                           *
  Subroutine Spline(X, Y, Ypp0, YppN, Pols)
!  *                                           *
!  *********************************************
!  * Calculates the Interpolation Spline polynomial 
!  * of data X(:), Y(:), and returns N-1 
!  * polynomials of degree 3.
!  *********************************************

    USE Linear

    Real (kind=DP), Intent (in) :: X(:), Y(:), Ypp0, YppN
    Type (Pol), Intent (out) :: Pols(:)

    Real (kind=DP) :: M(2:Size(X)-1,2:Size(X)-1), Ypp(Size(X)), a, &
         & b, c, d, h, hm
    Integer :: Npt, I

    Npt = Size(X)

    Ypp(1)   = Ypp0
    Ypp(Npt) = YppN
    M = 0.0_DP
    
    ! First and last rows
    h  = X(3) - X(2)
    hm = X(2) - X(1)
    M(2,2) = h+3.0_DP*hm
    M(2,3) = h
    Ypp(2) = 6.0_DP * ( &
         & (Y(3)-Y(2))/h - (Y(2)-Y(1))/hm ) - Ypp(1)*hm

    h  = X(Npt) - X(Npt-1)
    hm = X(Npt-1) - X(Npt-2)
    M(Npt-1,Npt-1) = h+3.0_DP*h
    M(Npt-1,Npt-2) = hm
    Ypp(Npt-1) = 6.0_DP * ( &
         & (Y(Npt)-Y(Npt-1))/h - (Y(Npt-1)-Y(Npt-2))/hm ) - Ypp(Npt)*h
    Do I = 3, Npt-2
       h  = X(I+1) - X(I)
       hm = X(I) - X(I-1)
       M(I,I)   = h + 3.0_DP*hm
       M(I,I-1) = hm
       M(I,I+1) = h

       Ypp(I) = 6.0_DP * ( &
         & (Y(I+1)-Y(I))/h - (Y(I)-Y(I-1))/hm )
    End Do

    CALL lusolve(M(2:Npt-1,2:Npt-1),Ypp(2:Npt-1))


    ! Now we have to construct the polynomials
    Do I = 1, Npt-1
       h  = X(I+1) - X(I)
       a = (Ypp(I+1)-Ypp(I))/(6.0_DP*h)
       b = Ypp(I)/2.0_DP
       c = (Y(I+1)-Y(I))/h - (Ypp(I+1)+2.0_DP*Ypp(I))*h/6.0_DP
       d = Y(I)

       CALL Init(Pols(I), 3)
       Pols(I)%Coef(0) = d-c*X(I)+b*X(I)**2-a*X(I)**3
       Pols(I)%Coef(1) = c-2.0_DP*b*X(I)+3.0_DP*a*X(I)**2
       Pols(I)%Coef(2) = b-3.0_DP*a*X(I)
       Pols(I)%Coef(3) = a
    End Do


    Return
  End Subroutine Spline

!  *********************************************
!  *                                           *
  Subroutine InitC(P, Ngrad)
!  *                                           *
!  *********************************************
!  * Init a polinomial of degree Ngrad
!  *********************************************

    Type (CmplxPol), Intent (inout) :: P
    Integer, Intent (in) :: Ngrad


    ALLOCATE(P%Coef(0:Ngrad))

    P%Coef = Cmplx(0.0_DP, kind=DPC)
    P%dg = Ngrad

    Return
  End Subroutine InitC

!  *********************************************
!  *                                           *
  Integer Function DegreeC(P)
!  *                                           *
!  *********************************************
!  * Returns the degree of a ploy.
!  *********************************************

    Type (CmplxPol) :: P

    DegreeC = P%dg

    Return
  End Function DegreeC

!  *********************************************
!  *                                           *
  Subroutine IgualCP(Pol1, Pol2)
!  *                                           *
!  *********************************************
!  * Equate two polinomials
!  *********************************************

    Type (CmplxPol), Intent(out) :: Pol1
    Type (CmplxPol), Intent(in) :: Pol2

    CALL Init(Pol1, Pol2%dg)
    Pol1%Coef = Pol2%Coef
    
    Return
  End Subroutine IgualCP

!  *********************************************
!  *                                           *
  Subroutine IgualCM(P, Coef)
!  *                                           *
!  *********************************************
!  * Iguala un polinomio a un Vector de 
!  * coeficientes
!  *********************************************

    Type (CmplxPol), Intent (out) :: P
    Complex (kind=DPC), Intent (in) :: Coef(:)

    CALL Init(P,Size(Coef)-1)
    P%Coef = Coef
    
    Return
  End Subroutine IgualCM

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function SumaCP(Pol1, Pol2) Result (Su)
!  *                                           *
!  *********************************************
!  * Suma Dos polinomios
!  * 
!  *********************************************

    Type (CmplxPol), Intent (in) :: Pol1, Pol2

    Integer :: I1, I2

    I1 = Pol1%dg
    I2 = Pol2%dg
    CALL Init(Su, Max(I1,I2))

    If (I1 > I2) Then
       Su%Coef = Pol1%Coef
       Su%Coef(0:I2) = Pol1%Coef(0:I2) + Pol2%Coef(0:I2)
    Else
       Su%Coef = Pol2%Coef
       Su%Coef(0:I1) = Pol1%Coef(0:I1) + Pol2%Coef(0:I1)
    End If


    Return
  End Function SumaCP

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function RestaCP(Pol1, Pol2) Result (Su)
!  *                                           *
!  *********************************************
!  * Resta Dos polinomios
!  * 
!  *********************************************

    Type (CmplxPol), Intent (in) :: Pol1, Pol2

    Integer :: I1, I2

    I1 = Pol1%dg
    I2 = Pol2%dg
    CALL Init(Su, Max(I1,I2))

    If (I1 > I2) Then
       Su%Coef = Pol1%Coef
       Su%Coef(0:I2) = Pol1%Coef(0:I2) - Pol2%Coef(0:I2)
    Else
       Su%Coef = -Pol2%Coef
       Su%Coef(0:I1) = Pol1%Coef(0:I1) - Pol2%Coef(0:I1)
    End If

    Return
  End Function RestaCP

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function MultCP(Pol1, Pol2) Result (MultP)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (CmplxPol), Intent (in) :: Pol1, Pol2

    Integer :: I, J, Igrado1, Igrado2

    Igrado1 = Pol1%dg
    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado1 + Igrado2)

    MultP%Coef = 0.0_DP
    Do I = 0, Igrado1
       Do J = 0, Igrado2
          MultP%Coef(I+J) = MultP%Coef(I+J) + Pol1%Coef(I)*Pol2%Coef(J)
       End Do
    End Do

    Return
  End Function MultCP

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function MultCCl(C, Pol2) Result (MultP)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (CmplxPol), Intent (in) :: Pol2
    Real (kind=DP), Intent (in) :: C

    Integer :: Igrado2

    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*Pol2%Coef

    Return
  End Function MultCCl

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function MultCCCl(C, Pol2) Result (MultP)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (CmplxPol), Intent (in) :: Pol2
    Complex (kind=DPC), Intent (in) :: C

    Integer :: Igrado2

    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*Pol2%Coef

    Return
  End Function MultCCCl

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function MultCCr(Pol2, C) Result (MultP)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (CmplxPol), Intent (in) :: Pol2
    Real (kind=DP), Intent (in)  :: C

    Integer :: Igrado2

    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*Pol2%Coef


    Return
  End Function MultCCr

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function MultCCCr(Pol2, C) Result (MultP)
!  *                                           *
!  *********************************************
!  * Product of two poly
!  *********************************************

    Type (CmplxPol), Intent (in) :: Pol2
    Complex (kind=DPC), Intent (in)  :: C

    Integer :: Igrado2

    Igrado2 = Pol2%dg
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*Pol2%Coef


    Return
  End Function MultCCCr

!  *********************************************
!  *                                           *
  Complex (kind=DPC) Function ValueC(Polin, X)
!  *                                           *
!  *********************************************
!  * Evaluates a Poly in X.
!  *********************************************
    
    Type (CmplxPol), Intent (in) :: Polin
    Complex (kind=DPC), Intent (in) :: X

    Integer :: I

    ValueC = Polin%Coef(Polin%dg)*X
    Do I = Polin%dg-1, 1, -1
       ValueC = (ValueC + Polin%Coef(I))*X
    End Do
    ValueC = ValueC + Polin%Coef(0)

    Return
  End Function ValueC

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function DerivC(P) Result (Deriv)
!  *                                           *
!  *********************************************
!  * Calculates de derivative of a poly.
!  *********************************************

    Type (CmplxPol), Intent (in) :: P

    Integer :: I

    CALL Init(Deriv, P%dg-1)
    
    Do I = 0, P%dg - 1
       Deriv%Coef(I) = Real(I+1, kind=DP)*P%Coef(I+1)
    End Do
    
    Return
  End Function DerivC

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function IntegraC(P, C) Result (Integra)
!  *                                           *
!  *********************************************
!  * Integra un polinomio
!  * 
!  *********************************************

    Type (CmplxPol), Intent (in) :: P
    Complex (kind=DPC), Optional :: C

    Integer :: I

    CALL Init(Integra, P%dg + 1)

    If (Present(C)) Then
       Integra%Coef(0) = C
    Else
       Integra%Coef(0) = 0.D0
    End If
       
    Do I = 1, Integra%dg
       Integra%Coef(I) = P%Coef(I-1)/Real(I,kind=DP)
    End Do
    
    Return
  End Function IntegraC

!  *********************************************
!  *                                           *
  Complex (kind=DPC) Function InterpolValueC(X, Y, Xo)
!  *                                           *
!  *********************************************
!  * Calculates the value of the Interpolation 
!  * polynomial of data X(:), Y(:) in the point 
!  * Xo.
!  *********************************************

    Complex (kind=DPC), Intent (in) :: X(:), Y(:), Xo

    Complex (kind=DPC) :: a(Size(X),Size(X))
    Integer :: I, J, Npoints

    Npoints = Size(X)
    a = (0.0_DP, 0.0_DP) 
    a(:,1) = Y(:)
    Do J = 2, Npoints
       Do I = J, Npoints
          a(I,J) = ( (Xo-X(I-J+1))*a(I,J-1) + (X(I)-Xo)*a(I-1,J-1) ) &
               & / (X(I) - X(I-J+1))
       End Do
    End Do
    
    InterpolValueC = a(Npoints, Npoints)

    Return
  End Function InterpolValueC

!  *********************************************
!  *                                           *
  Type (CmplxPol) Function InterpolC(X, Y) Result (InterPol)
!  *                                           *
!  *********************************************
!  * Calculates the Interpolation polynomial of
!  * data X(:), Y(:) using divided diference and
!  * the Newton form of the Interpolating 
!  * polynomial.
!  * 
!  * NOTE:
!  * In general this problem has a lot of 
!  * numerical errors, so it is much more 
!  * interesting to calculate the value of the 
!  * interpolating polynomial using the function 
!  * InterpolValue of this Module. 
!  *********************************************

    Complex (kind=DPC), Intent (in) :: X(:), Y(:)
    
    Type (CmplxPol) :: Base, Mul
    Complex (kind=DPC) :: a(Size(X), Size(X))
    Integer :: I, J

    CALL Init(Interpol, Size(X)-1)
    CALL Init(Base,0)
    CALL Init(Mul,1)

    ! First we construct the Divided diferences
    a = 0.0_DP 
    a(:,1) = Y(:)
    Do J = 2, Interpol%dg + 1
       Do I = J, Interpol%dg + 1
          a(I,J) = ( a(I,J-1) - a(I-1,J-1) ) / (X(I) - X(I-J+1))
       End Do
    End Do

    ! Now construc the Polynomial in the Newton form.
    Interpol%Coef(0) = Y(1)
    Base%Coef(0) = 1.0_DP
    Do I = 1, Interpol%dg
       Mul%Coef(1) = 1.0_DP
       Mul%Coef(0) = -X(I)
       Base = Base * Mul

       Interpol = Interpol + a(I+1,I+1)*Base
    End Do

    Return
  End Function InterpolC


End MODULE Polynomial
