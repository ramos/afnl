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

! ***********************************************************
! *
MODULE Polynomial
! *
! ***********************************************************

  USE NumTypes
  USE Error


  Type Pol
     Real (kind=DP), Pointer :: Coef(:) => Null()
     Integer :: dg
  End Type Pol


  Interface Assignment (=)
     Module Procedure IgualM, IgualP
  End Interface

  Interface Operator (+)
     Module Procedure SumaP
  End Interface

  Interface Operator (-)
     Module Procedure RestaP
  End Interface

  Interface Operator (*)
     Module Procedure MultP, MultCl, MultCr
  End Interface
     

  Private IgualR, IgualM, IgualP, SumaP, RestaP, MultP, MultCl, MultCr


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

    Type (Pol) :: PP1, PP2 ! Local copies of the input poly

    I1 = Pol1%dg
    CALL Init(PP1, I1)
    PP1 = Pol1
    I2 = Pol2%dg
    CALL Init(PP2, I2)
    PP2 = Pol2
    CALL Init(Su, Max(I1,I2))

    If (I1 > I2) Then
       Su%Coef = Pol1%Coef
       Su%Coef(0:I2) = Pol1%Coef(0:I2) + Pol2%Coef(0:I2)
    Else
       Su%Coef = Pol2%Coef
       Su%Coef(0:I1) = Pol1%Coef(0:I1) + Pol2%Coef(0:I1)
    End If

    Deallocate(PP1%Coef, PP2%Coef)

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

    Type (Pol) :: PP1, PP2 ! Local copies of the input poly

    I1 = Pol1%dg
    CALL Init(PP1, I1)
    PP1 = Pol1
    I2 = Pol2%dg
    CALL Init(PP2, I2)
    PP2 = Pol2
    CALL Init(Su, Max(I1,I2))

    If (I1 > I2) Then
       Su%Coef = PP1%Coef
       Su%Coef(0:I2) = PP1%Coef(0:I2) - PP2%Coef(0:I2)
    Else
       Su%Coef = -PP2%Coef
       Su%Coef(0:I1) = PP1%Coef(0:I1) - PP2%Coef(0:I1)
    End If

    Deallocate(PP1%Coef, PP2%Coef)

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

    Type (Pol) :: PP1, PP2 ! Local copies of the input poly

    Igrado1 = Pol1%dg
    CALL Init(PP1, Igrado1)
    PP1 = Pol1
    Igrado2 = Pol2%dg
    CALL Init(PP2, Igrado2)
    PP2 = Pol2
    CALL Init(MultP, Igrado1 + Igrado2)


    MultP%Coef = 0.0_DP
    Do I = 0, Igrado1
       Do J = 0, Igrado2
          MultP%Coef(I+J) = MultP%Coef(I+J) + PP1%Coef(I)*PP2%Coef(J)
       End Do
    End Do

    Deallocate(PP1%Coef, PP2%Coef)

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

    Type (Pol) :: PP2

    Igrado2 = Pol2%dg
    CALL Init(PP2, Igrado2)
    PP2 = Pol2
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*PP2%Coef

    Deallocate(PP2%Coef)

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

    Type (Pol) :: PP2

    Igrado2 = Pol2%dg
    CALL Init(PP2, Igrado2)
    PP2 = Pol2
    CALL Init(MultP, Igrado2)

    MultP%Coef = C*PP2%Coef

    Deallocate(PP2%Coef)

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

    CALL Init(Deriv, P%dg-1)
    
    Do I = 0, P%dg - 1
       Deriv%Coef(I) = (I+1)*P%Coef(I+1)
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
    
    Type (Pol) :: Base, Mul, P
    Real (kind=DP) :: Coef(Size(X)), a(Size(X), Size(X))

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

End MODULE Polynomial
