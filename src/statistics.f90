!
! MODULE with useful statistical-related functions
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

! ***********************************************************
! *
MODULE Statistics
! *
! *********************************************************** 

  USE NumTypes
  USE Constants, ONLY: TWOPI_DP, TWOPI_SP
  USE Error  
  USE Linear

  IMPLICIT NONE

  Interface Mean
     Module Procedure Mean_DP, Mean_SP
  End Interface

  Interface Var
     Module Procedure Var_DP, Var_SP
  End Interface

  Interface Stddev
     Module Procedure Stddev_DP, Stddev_SP
  End Interface

  Interface Moment
     Module Procedure Moment_DP, Moment_SP
  End Interface

  Interface Normal 
     Module Procedure NormalS, NormalV, NormalS2, NormalV2, &
          & NormalS_SP, NormalV_SP, NormalS2_SP, NormalV2_SP
  End Interface

  Interface ChiSqr
     Module Procedure ChiSqr_DP, ChiSqr_SP
  End Interface

  Interface Laplace
     Module Procedure Laplace_DP, Laplace_SP,Laplace2_DP, Laplace2_SP
  End Interface

  Interface FishTipp
     Module Procedure FishTipp_DP, FishTipp_SP,FishTipp2_DP, FishTipp2_SP
  End Interface

  Interface Histogram
     Module Procedure Histogram_SP, Histogram_DP
  End Interface

  Interface Irand
     Module Procedure Irand_S, Irand_V
  End Interface

  Interface EstBstrp
     Module Procedure EstBstrp_H, EstBstrp
  End Interface


  Interface LinearReg
     Module Procedure LinearReg_DP, LinearReg_SP, &
          & LinearReg_Pol_DP, LinearReg_Pol_SP, &
          & MultiLinearReg_DP, MultiLinearReg_SP
  End Interface


  Private NormalS, NormalV, NormalS2, NormalV2, &
       & NormalS_SP, NormalV_SP, NormalS2_SP, NormalV2_SP, &
       & Stddev_DP, Var_DP, Mean_DP, Moment_DP, &
       & Stddev_SP, Var_SP, Mean_SP, Moment_SP, ChiSqr_SP, &
       & ChiSqr_DP, LinearReg_DP, LinearReg_SP, LinearReg_Pol_DP, &
       & LinearReg_Pol_SP, MultiLinearReg_DP, MultiLinearReg_SP, &
       & Laplace_DP, Laplace_SP,Laplace2_DP, Laplace2_SP, &
       & FishTipp_DP, FishTipp_SP,FishTipp2_DP, FishTipp2_SP, &
       & Irand_S, Irand_V, EstBstrp_H

CONTAINS

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Mean_DP(X)
!  *                                           *
!  *********************************************
!  * Returns the mean value of the elements in 
!  * the vector X(:)
!  *********************************************

    Real (kind=DP), Intent (in) :: X(:)

    Mean_DP = Sum(X) / Size(X)

    Return
  End Function Mean_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Var_DP(X)
!  *                                           *
!  *********************************************
!  * Returns the Variance of the elements in 
!  * the vector X(:)
!  *********************************************

    Real (kind=DP), Intent (in) :: X(:)
    Real (kind=DP) :: Nmed

    Nmed = Mean(X)
    Var_DP = Sum((X(:)-Nmed)**2) / Size(X)

    Return
  End Function Var_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Stddev_DP(X)
!  *                                           *
!  *********************************************
!  * Returns the standard deviantion of the 
!  * elements in the vector X(:)
!  *********************************************

    Real (kind=DP), Intent (in) :: X(:)

    Stddev_DP = Sqrt(Var(X))

    Return
  End Function Stddev_DP

!  *********************************************
!  *                                           *
  Real (kind=DP) Function Moment_DP(X, k)
!  *                                           *
!  *********************************************
!  * Returns the k moment of the data in the 
!  * vector X(:)
!  *********************************************

    Real (kind=DP), Intent (in) :: X(:)
    Integer, Intent (in) :: k
    Real (kind=DP) :: Nmed

    Nmed = Mean(X)
    Moment_DP = Sum((X(:)-Nmed)**k) / Size(X)

    Return
  End Function Moment_DP

!  *********************************************
!  *                                           *
  Subroutine NormalS(X)
!  *                                           *
!  *********************************************
!  * Returns a number with N(0,1) distribution in
!  * the intent (out) Real DP variable X
!  *********************************************

    Real (kind=DP), Intent(out) :: X
    Real (kind=DP) :: U1, U2
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Sqrt(-2.0_DP*Log(U1)) * Cos(TWOPI_DP*U2)


    Return
  End Subroutine  NormalS

!  *********************************************
!  *                                           *
  Subroutine NormalV(X)
!  *                                           *
!  *********************************************
!  * Returns a vector of numbers with N(0,1) 
!  * distribution in the intent (out) Real DP 
!  * variable X
!  *********************************************

    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP) :: U1(Size(X)), U2(Size(X))
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Sqrt(-2.0_DP*Log(U1)) * Cos(TWOPI_DP*U2)


    Return
  End Subroutine  NormalV

!  *********************************************
!  *                                           *
  Subroutine NormalS2(X, Rmed, Rsig)
!  *                                           *
!  *********************************************
!  * Returns a number with N(Rmed, Rsig) 
!  * distribution in the intent (out) Real DP 
!  * variable X
!  *********************************************

    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent (in) :: Rmed, Rsig
    Real (kind=DP) :: U1, U2
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Rsig * Sqrt(-2.0_DP*Log(U1)) * Cos(TWOPI_DP*U2) + Rmed 


    Return
  End Subroutine  NormalS2

!  *********************************************
!  *                                           *
  Subroutine NormalV2(X, Rmed, Rsig)
!  *                                           *
!  *********************************************
!  * Returns a vector of numbers with N(Rmed,Rsig) 
!  * distribution in the intent (out) Real DP 
!  * variable X
!  *********************************************

    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent(in) :: Rmed, Rsig
    Real (kind=DP) :: U1(Size(X)), U2(Size(X))
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Rsig * Sqrt(-2.0_DP*Log(U1)) * Cos(TWOPI_DP*U2) + Rmed

    Return
  End Subroutine  NormalV2

!  *********************************************
!  *                                           *
  Real (kind=SP) Function Mean_SP(X)
!  *                                           *
!  *********************************************
!  * Returns the mean value of the elements in 
!  * the vector X(:)
!  *********************************************

    Real (kind=SP), Intent (in) :: X(:)

    Mean_SP = Sum(X) / Size(X)

    Return
  End Function Mean_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function Var_SP(X)
!  *                                           *
!  *********************************************
!  * Returns the Variance of the elements in 
!  * the vector X(:)
!  *********************************************

    Real (kind=SP), Intent (in) :: X(:)
    Real (kind=SP) :: Nmed

    Nmed = Mean(X)
    Var_SP = Sum((X(:)-Nmed)**2) / Size(X)

    Return
  End Function Var_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function Stddev_SP(X)
!  *                                           *
!  *********************************************
!  * Returns the standard deviantion of the 
!  * elements in the vector X(:)
!  *********************************************

    Real (kind=SP), Intent (in) :: X(:)

    Stddev_SP = Sqrt(Var(X))

    Return
  End Function Stddev_SP

!  *********************************************
!  *                                           *
  Real (kind=SP) Function Moment_SP(X, k)
!  *                                           *
!  *********************************************
!  * Returns the k moment of the data in the 
!  * vector X(:)
!  *********************************************

    Real (kind=SP), Intent (in) :: X(:)
    Integer, Intent (in) :: k
    Real (kind=SP) :: Nmed

    Nmed = Mean(X)
    Moment_SP = Sum((X(:)-Nmed)**k) / Size(X)

    Return
  End Function Moment_SP

!  *********************************************
!  *                                           *
  Subroutine NormalS_SP(X)
!  *                                           *
!  *********************************************
!  * Returns a number with N(0,1) distribution in
!  * the intent (out) Real SP variable X
!  *********************************************

    Real (kind=SP), Intent(out) :: X
    Real (kind=SP) :: U1, U2
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Sqrt(-2.0_SP*Log(U1)) * Cos(TWOPI_SP*U2)


    Return
  End Subroutine  NormalS_SP

!  *********************************************
!  *                                           *
  Subroutine NormalV_SP(X)
!  *                                           *
!  *********************************************
!  * Returns a vector of numbers with N(0,1) 
!  * distribution in the intent (out) Real SP 
!  * variable X
!  *********************************************

    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP) :: U1(Size(X)), U2(Size(X))
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Sqrt(-2.0_SP*Log(U1)) * Cos(TWOPI_SP*U2)


    Return
  End Subroutine  NormalV_SP

!  *********************************************
!  *                                           *
  Subroutine NormalS2_SP(X, Rmed, Rsig)
!  *                                           *
!  *********************************************
!  * Returns a number with N(Rmed, Rsig) 
!  * distribution in the intent (out) Real SP 
!  * variable X
!  *********************************************

    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent (in) :: Rmed, Rsig
    Real (kind=SP) :: U1, U2
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Rsig * Sqrt(-2.0_SP*Log(U1)) * Cos(TWOPI_SP*U2) + Rmed 


    Return
  End Subroutine  NormalS2_SP

!  *********************************************
!  *                                           *
  Subroutine NormalV2_SP(X, Rmed, Rsig)
!  *                                           *
!  *********************************************
!  * Returns a vector of numbers with N(Rmed,Rsig) 
!  * distribution in the intent (out) Real SP 
!  * variable X
!  *********************************************

    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent(in) :: Rmed, Rsig
    Real (kind=SP) :: U1(Size(X)), U2(Size(X))
    
    CALL Random_Number(U1)
    CALL Random_Number(U2)

    X = Rsig * Sqrt(-2.0_SP*Log(U1)) * Cos(TWOPI_SP*U2) + Rmed

    Return
  End Subroutine  NormalV2_SP

!  *********************************************
!  *                                           *
  Function ChiSqr_SP(V, D)
!  *                                           *
!  *********************************************
!  * Returns the Maximum likehood value for a set
!  * of points V(:) with (optional) errors D(:)
!  *********************************************

    Real (kind=SP), Intent (in) :: V(:)
    Real (kind=SP), Intent (in), Optional :: D(:)
    Real (kind=SP) :: ChiSqr_SP

    Real (kind=SP) :: Serr, Sval
    Integer :: Nval

    Nval = Size(V)
    Serr = 0.0_SP
    Sval = 0.0_SP
    If (Present(D)) Then
       Serr = Sum(1.0_SP/D(:)**2)
       Sval = Sum(V(:)/D(:)**2)
       ChiSqr_SP = 1.0_SP/Serr * Sval
       Return
    Else
       ChiSqr_SP = Sum(V(:))/Nval
    End If


    Return
  End Function ChiSqr_SP

!  *********************************************
!  *                                           *
  Function ChiSqr_DP(V, D)
!  *                                           *
!  *********************************************
!  * Returns the Maximum likehood value for a set
!  * of points V(:) with (optional) errors D(:)
!  *********************************************

    Real (kind=DP), Intent (in) :: V(:)
    Real (kind=DP), Intent (in), Optional :: D(:)
    Real (kind=DP) :: ChiSqr_DP

    Real (kind=DP) :: Serr, Sval
    Integer :: Nval

    Nval = Size(V)
    Serr = 0.0_DP
    Sval = 0.0_DP
    If (Present(D)) Then
       Serr = Sum(1.0_DP/D(:)**2)
       Sval = Sum(V(:)/D(:)**2)
       ChiSqr_DP = 1.0_DP/Serr * Sval
       Return
    Else
       ChiSqr_DP = Sum(V(:))/Nval
    End If


    Return
  End Function ChiSqr_DP


!  *********************************************
!  *                                           *
  Subroutine Histogram_SP(Val, Ndiv, Ntics, Vmin, Vmax, h)
!  *                                           *
!  *********************************************
!  * 
!  * 
!  *********************************************

    Real (kind=SP), Intent (in) :: Val(:)
    Integer, Intent (in) :: Ndiv
    Real (kind=SP), Intent (out) :: Vmin, Vmax, h
    Integer, Intent (out) :: Ntics(Ndiv)

    Integer :: I, K
    
    Vmin = MinVal(Val)
    Vmax = MaxVal(Val)
    h = (Vmax - Vmin)/Real(Ndiv,kind=SP)

    Ntics = 0
    Do I = 1, Size(Val)
       K = Int((Val(I) - Vmin)/h) + 1
       Ntics(K) = Ntics(K) + 1
    End Do

    Return
  End Subroutine Histogram_SP

!  *********************************************
!  *                                           *
  Subroutine Histogram_DP(Val, Ndiv, Ntics, Vmin, Vmax, h)
!  *                                           *
!  *********************************************
!  * 
!  * 
!  *********************************************

    Real (kind=DP), Intent (in) :: Val(:)
    Integer, Intent (in) :: Ndiv
    Real (kind=DP), Intent (out) :: Vmin, Vmax, h
    Integer, Intent (out) :: Ntics(Ndiv)

    Integer :: I, K
    
    Vmin = MinVal(Val)
    Vmax = MaxVal(Val)
    h = (Vmax - Vmin)/Real(Ndiv,kind=DP)

    Ntics = 0
    Do I = 1, Size(Val)
       K = Int((Val(I) - Vmin)/h) + 1
       Ntics(K) = Ntics(K) + 1
    End Do

    Return
  End Subroutine Histogram_DP

!  *********************************************
!  *                                           *
  Subroutine LinearReg_DP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a function \sum_i Coef(i)*Func(X,i).
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr is returned
!  *********************************************
    
    Real (kind=DP), Intent(in) :: X(:), Y(:), Yerr(:)
    Real (kind=DP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=DP) :: Sm(Size(Coef), Size(Coef)), Kv(Size(Coef)), Fval
    Integer :: Ipiv(Size(Coef)), Nparm, Npoints, I, J, K, Idet

    Interface
       Function Func(Xx, i)
         
         USE NumTypes

         Real (kind=DP), Intent (in) :: Xx
         Integer, Intent (in) :: i
         Real (kind=DP) :: Func

       End Function Func
    End Interface


    Nparm = Size(Coef)
    Npoints = Size(X)
    Sm = 0.0_DP
    Kv = 0.0_DP
    Do I = 1, Nparm
       Do K = 1, Npoints
          Kv(I) = Kv(I) + Y(K)*Func(X(K), I)/(Yerr(K))**2
       End Do
       Do J = 1, Nparm
          Do K = 1, Npoints
             Sm(I,J) = Sm(I,J) + Func(X(K), I)*Func(X(K), J)/(Yerr(K)**2)
          End Do
       End Do
    End Do


    CALL LU(Sm, Ipiv, Idet)
    Do I = 1, Nparm
       Coef(I) = Kv(Ipiv(I))
    End Do
    Do I = 2, Nparm
       Coef(I) = Coef(I) - Dot_Product(Sm(I,1:I-1),Coef(1:I-1))
    End Do
    Coef(Nparm) = Coef(Nparm) / Sm(Nparm, Nparm)
    Do I = Nparm - 1, 1, -1
       Coef(I) = Coef(I) - Dot_Product(Sm(I, I+1:Nparm),Coef(I+1:Nparm))
       Coef(I) = Coef(I) / Sm(I, I)
    End Do

    Do J = 1, Nparm
       Do I = 1, Nparm
          If (Ipiv(I) == J) Then
             Kv(I) = 1.0_DP
          Else
             Kv(I) = 0.0_DP
          End If
       End Do

       Do I = 2, Nparm
          Kv(I) = Kv(I) - Dot_Product(Sm(I,1:I-1),Kv(1:I-1))
       End Do
       Kv(Nparm) = Kv(Nparm) / Sm(Nparm, Nparm)
       Do I = Nparm - 1, J, -1
          Kv(I) = Kv(I) - Dot_Product(Sm(I, I+1:Nparm),Kv(I+1:Nparm))
          Kv(I) = Kv(I) / Sm(I, I)
       End Do
       Cerr(J) = Sqrt(Kv(J))
    End Do

    ChisqrV = 0.0_DP
    Do I = 1, Npoints
       Fval = 0.0_DP
       Do J = 1, Nparm
          Fval = Fval + Coef(J)*Func(X(I),J)
       End Do
       ChisqrV = ChisqrV + ((Y(I) - Fval)/(Yerr(I)))**2
    End Do
    
    ChisqrV = ChisqrV / Real(Npoints - Nparm,kind=DP)


    Return
  End Subroutine LinearReg_DP

!  *********************************************
!  *                                           *
  Subroutine LinearReg_Pol_DP(X, Y, Yerr, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a polynomial of degree Size(Coef).
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr is returned
!  *********************************************
    
    Real (kind=DP), Intent(in) :: X(:), Y(:), Yerr(:)
    Real (kind=DP), Intent(out) :: Coef(:), Cerr(:), ChisqrV
    
    Real (kind=DP) :: Sm(Size(Coef), Size(Coef)), Kv(Size(Coef)), Fval
    Integer :: Ipiv(Size(Coef)), Nparm, Npoints, I, J, K, Idet


    Nparm = Size(Coef)
    Npoints = Size(X)
    Sm = 0.0_DP
    Kv = 0.0_DP
    Do I = 1, Nparm
       Do K = 1, Npoints
          Kv(I) = Kv(I) + Y(K)*(X(K)**I)/(Yerr(K))**2
       End Do
       Do J = 1, Nparm
          Do K = 1, Npoints
             Sm(I,J) = Sm(I,J) + (X(K)**I)*(X(K)**J)/(Yerr(K)**2)
          End Do
       End Do
    End Do

    CALL LU(Sm, Ipiv, Idet)

    Do I = 1, Nparm
       Coef(I) = Kv(Ipiv(I))
    End Do
    Do I = 2, Nparm
       Coef(I) = Coef(I) - Dot_Product(Sm(I,1:I-1),Coef(1:I-1))
    End Do
    Coef(Nparm) = Coef(Nparm) / Sm(Nparm, Nparm)
    Do I = Nparm - 1, 1, -1
       Coef(I) = Coef(I) - Dot_Product(Sm(I, I+1:Nparm),Coef(I+1:Nparm))
       Coef(I) = Coef(I) / Sm(I, I)
    End Do

    Do J = 1, Nparm
       Do I = 1, Nparm
          If (Ipiv(I) == J) Then
             Kv(I) = 1.0_DP
          Else
             Kv(I) = 0.0_DP
          End If
       End Do

       Do I = 2, Nparm
          Kv(I) = Kv(I) - Dot_Product(Sm(I,1:I-1),Kv(1:I-1))
       End Do
       Kv(Nparm) = Kv(Nparm) / Sm(Nparm, Nparm)
       Do I = Nparm - 1, J, -1
          Kv(I) = Kv(I) - Dot_Product(Sm(I, I+1:Nparm),Kv(I+1:Nparm))
          Kv(I) = Kv(I) / Sm(I, I)
       End Do
       Cerr(J) = Sqrt(Kv(J))
    End Do

    ChisqrV = 0.0_DP
    Do I = 1, Npoints
       Fval = 0.0_DP
       Do J = 1, Nparm
          Fval = Fval + Coef(J)*(X(I)**J)
       End Do
       ChisqrV = ChisqrV + ((Y(I) - Fval)/(Yerr(I)))**2
    End Do
    
    ChisqrV = ChisqrV / Real(Npoints - Nparm,kind=DP)

    
    Return
  End Subroutine LinearReg_Pol_DP

!  *********************************************
!  *                                           *
  Subroutine LinearReg_SP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a function \sum_i Coef(i)*Func(X,i).
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr is returned
!  *********************************************
    
    Real (kind=SP), Intent(in) :: X(:), Y(:), Yerr(:)
    Real (kind=SP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=SP) :: Sm(Size(Coef), Size(Coef)), Kv(Size(Coef)), Fval
    Integer :: Ipiv(Size(Coef)), Nparm, Npoints, I, J, K, Idet

    Interface
       Function Func(Xx, i)
         
         USE NumTypes

         Real (kind=SP), Intent (in) :: Xx
         Integer, Intent (in) :: i
         Real (kind=SP) :: Func

       End Function Func
    End Interface


    Nparm = Size(Coef)
    Npoints = Size(X)
    Sm = 0.0_SP
    Kv = 0.0_SP
    Do I = 1, Nparm
       Do K = 1, Npoints
          Kv(I) = Kv(I) + Y(K)*Func(X(K), I)/(Yerr(K))**2
       End Do
       Do J = 1, Nparm
          Do K = 1, Npoints
             Sm(I,J) = Sm(I,J) + Func(X(K), I)*Func(X(K), J)/(Yerr(K)**2)
          End Do
       End Do
    End Do


    CALL LU(Sm, Ipiv, Idet)
    Do I = 1, Nparm
       Coef(I) = Kv(Ipiv(I))
    End Do
    Do I = 2, Nparm
       Coef(I) = Coef(I) - Dot_Product(Sm(I,1:I-1),Coef(1:I-1))
    End Do
    Coef(Nparm) = Coef(Nparm) / Sm(Nparm, Nparm)
    Do I = Nparm - 1, 1, -1
       Coef(I) = Coef(I) - Dot_Product(Sm(I, I+1:Nparm),Coef(I+1:Nparm))
       Coef(I) = Coef(I) / Sm(I, I)
    End Do

    Do J = 1, Nparm
       Do I = 1, Nparm
          If (Ipiv(I) == J) Then
             Kv(I) = 1.0_SP
          Else
             Kv(I) = 0.0_SP
          End If
       End Do

       Do I = 2, Nparm
          Kv(I) = Kv(I) - Dot_Product(Sm(I,1:I-1),Kv(1:I-1))
       End Do
       Kv(Nparm) = Kv(Nparm) / Sm(Nparm, Nparm)
       Do I = Nparm - 1, J, -1
          Kv(I) = Kv(I) - Dot_Product(Sm(I, I+1:Nparm),Kv(I+1:Nparm))
          Kv(I) = Kv(I) / Sm(I, I)
       End Do
       Cerr(J) = Sqrt(Kv(J))
    End Do

    ChisqrV = 0.0_SP
    Do I = 1, Npoints
       Fval = 0.0_SP
       Do J = 1, Nparm
          Fval = Fval + Coef(J)*Func(X(I),J)
       End Do
       ChisqrV = ChisqrV + ((Y(I) - Fval)/(Yerr(I)))**2
    End Do
    
    ChisqrV = ChisqrV / Real(Npoints - Nparm,kind=SP)


    Return
  End Subroutine LinearReg_SP

!  *********************************************
!  *                                           *
  Subroutine LinearReg_Pol_SP(X, Y, Yerr, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a polynomial of degree Size(Coef).
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr is returned
!  *********************************************
    
    Real (kind=SP), Intent(in) :: X(:), Y(:), Yerr(:)
    Real (kind=SP), Intent(out) :: Coef(:), Cerr(:), ChisqrV
    
    Real (kind=SP) :: Sm(Size(Coef), Size(Coef)), Kv(Size(Coef)), Fval
    Integer :: Ipiv(Size(Coef)), Nparm, Npoints, I, J, K, Idet


    Nparm = Size(Coef)
    Npoints = Size(X)
    Sm = 0.0_SP
    Kv = 0.0_SP
    Do I = 1, Nparm
       Do K = 1, Npoints
          Kv(I) = Kv(I) + Y(K)*(X(K)**I)/(Yerr(K))**2
       End Do
       Do J = 1, Nparm
          Do K = 1, Npoints
             Sm(I,J) = Sm(I,J) + (X(K)**I)*(X(K)**J)/(Yerr(K)**2)
          End Do
       End Do
    End Do

    CALL LU(Sm, Ipiv, Idet)

    Do I = 1, Nparm
       Coef(I) = Kv(Ipiv(I))
    End Do
    Do I = 2, Nparm
       Coef(I) = Coef(I) - Dot_Product(Sm(I,1:I-1),Coef(1:I-1))
    End Do
    Coef(Nparm) = Coef(Nparm) / Sm(Nparm, Nparm)
    Do I = Nparm - 1, 1, -1
       Coef(I) = Coef(I) - Dot_Product(Sm(I, I+1:Nparm),Coef(I+1:Nparm))
       Coef(I) = Coef(I) / Sm(I, I)
    End Do

    Do J = 1, Nparm
       Do I = 1, Nparm
          If (Ipiv(I) == J) Then
             Kv(I) = 1.0_SP
          Else
             Kv(I) = 0.0_SP
          End If
       End Do

       Do I = 2, Nparm
          Kv(I) = Kv(I) - Dot_Product(Sm(I,1:I-1),Kv(1:I-1))
       End Do
       Kv(Nparm) = Kv(Nparm) / Sm(Nparm, Nparm)
       Do I = Nparm - 1, J, -1
          Kv(I) = Kv(I) - Dot_Product(Sm(I, I+1:Nparm),Kv(I+1:Nparm))
          Kv(I) = Kv(I) / Sm(I, I)
       End Do
       Cerr(J) = Sqrt(Kv(J))
    End Do

    ChisqrV = 0.0_SP
    Do I = 1, Npoints
       Fval = 0.0_SP
       Do J = 1, Nparm
          Fval = Fval + Coef(J)*(X(I)**J)
       End Do
       ChisqrV = ChisqrV + ((Y(I) - Fval)/(Yerr(I)))**2
    End Do
    
    ChisqrV = ChisqrV / Real(Npoints - Nparm,kind=SP)
    
    Return
  End Subroutine LinearReg_Pol_SP


!  *********************************************
!  *                                           *
  Subroutine MultiLinearReg_DP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:,:), Y(:)), this routine 
!  * fit the points to a function \sum_i Coef(i)*Func(X(:),i).
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr is returned
!  *********************************************
    
    Real (kind=DP), Intent(in) :: X(:,:), Y(:), Yerr(:)
    Real (kind=DP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=DP) :: Sm(Size(Coef), Size(Coef)), Kv(Size(Coef)), Fval
    Integer :: Ipiv(Size(Coef)), Nparm, Npoints, I, J, K, Idet

    Interface
       Function Func(Xx, i)
         
         USE NumTypes

         Real (kind=DP), Intent (in) :: Xx(:)
         Integer, Intent (in) :: i
         Real (kind=DP) :: Func

       End Function Func
    End Interface


    Nparm = Size(Coef)
    Npoints = Size(X,1)

    Sm = 0.0_DP
    Kv = 0.0_DP
    Do I = 1, Nparm
       Do K = 1, Npoints
          Kv(I) = Kv(I) + Y(K)*Func(X(K,:), I)/(Yerr(K))**2
       End Do
       Do J = 1, Nparm
          Do K = 1, Npoints
             Sm(I,J) = Sm(I,J) + Func(X(K,:), I)*Func(X(K,:), J)/(Yerr(K)**2)
          End Do
       End Do
    End Do


    CALL LU(Sm, Ipiv, Idet)
    Do I = 1, Nparm
       Coef(I) = Kv(Ipiv(I))
    End Do
    Do I = 2, Nparm
       Coef(I) = Coef(I) - Dot_Product(Sm(I,1:I-1),Coef(1:I-1))
    End Do
    Coef(Nparm) = Coef(Nparm) / Sm(Nparm, Nparm)
    Do I = Nparm - 1, 1, -1
       Coef(I) = Coef(I) - Dot_Product(Sm(I, I+1:Nparm),Coef(I+1:Nparm))
       Coef(I) = Coef(I) / Sm(I, I)
    End Do

    Do J = 1, Nparm
       Do I = 1, Nparm
          If (Ipiv(I) == J) Then
             Kv(I) = 1.0_DP
          Else
             Kv(I) = 0.0_DP
          End If
       End Do

       Do I = 2, Nparm
          Kv(I) = Kv(I) - Dot_Product(Sm(I,1:I-1),Kv(1:I-1))
       End Do
       Kv(Nparm) = Kv(Nparm) / Sm(Nparm, Nparm)
       Do I = Nparm - 1, J, -1
          Kv(I) = Kv(I) - Dot_Product(Sm(I, I+1:Nparm),Kv(I+1:Nparm))
          Kv(I) = Kv(I) / Sm(I, I)
       End Do
       Cerr(J) = Sqrt(Kv(J))
    End Do

    ChisqrV = 0.0_DP
    Do I = 1, Npoints
       Fval = 0.0_DP
       Do J = 1, Nparm
          Fval = Fval + Coef(J)*Func(X(I,:),J)
       End Do
       ChisqrV = ChisqrV + ((Y(I) - Fval)/(Yerr(I)))**2
    End Do
    
    ChisqrV = ChisqrV / Real(Npoints - Nparm,kind=DP)


    Return
  End Subroutine MultiLinearReg_DP

!  *********************************************
!  *                                           *
  Subroutine MultiLinearReg_SP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:,:), Y(:)), this routine 
!  * fit the points to a function \sum_i Coef(i)*Func(X(:),i).
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr is returned
!  *********************************************
    
    Real (kind=SP), Intent(in) :: X(:,:), Y(:), Yerr(:)
    Real (kind=SP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=SP) :: Sm(Size(Coef), Size(Coef)), Kv(Size(Coef)), Fval
    Integer :: Ipiv(Size(Coef)), Nparm, Npoints, I, J, K, Idet

    Interface
       Function Func(Xx, i)
         
         USE NumTypes

         Real (kind=SP), Intent (in) :: Xx(:)
         Integer, Intent (in) :: i
         Real (kind=SP) :: Func

       End Function Func
    End Interface


    Nparm = Size(Coef)
    Npoints = Size(X,1)

    Sm = 0.0_SP
    Kv = 0.0_SP
    Do I = 1, Nparm
       Do K = 1, Npoints
          Kv(I) = Kv(I) + Y(K)*Func(X(K,:), I)/(Yerr(K))**2
       End Do
       Do J = 1, Nparm
          Do K = 1, Npoints
             Sm(I,J) = Sm(I,J) + Func(X(K,:), I)*Func(X(K,:), J)/(Yerr(K)**2)
          End Do
       End Do
    End Do


    CALL LU(Sm, Ipiv, Idet)
    Do I = 1, Nparm
       Coef(I) = Kv(Ipiv(I))
    End Do
    Do I = 2, Nparm
       Coef(I) = Coef(I) - Dot_Product(Sm(I,1:I-1),Coef(1:I-1))
    End Do
    Coef(Nparm) = Coef(Nparm) / Sm(Nparm, Nparm)
    Do I = Nparm - 1, 1, -1
       Coef(I) = Coef(I) - Dot_Product(Sm(I, I+1:Nparm),Coef(I+1:Nparm))
       Coef(I) = Coef(I) / Sm(I, I)
    End Do

    Do J = 1, Nparm
       Do I = 1, Nparm
          If (Ipiv(I) == J) Then
             Kv(I) = 1.0_SP
          Else
             Kv(I) = 0.0_SP
          End If
       End Do

       Do I = 2, Nparm
          Kv(I) = Kv(I) - Dot_Product(Sm(I,1:I-1),Kv(1:I-1))
       End Do
       Kv(Nparm) = Kv(Nparm) / Sm(Nparm, Nparm)
       Do I = Nparm - 1, J, -1
          Kv(I) = Kv(I) - Dot_Product(Sm(I, I+1:Nparm),Kv(I+1:Nparm))
          Kv(I) = Kv(I) / Sm(I, I)
       End Do
       Cerr(J) = Sqrt(Kv(J))
    End Do

    ChisqrV = 0.0_SP
    Do I = 1, Npoints
       Fval = 0.0_SP
       Do J = 1, Nparm
          Fval = Fval + Coef(J)*Func(X(I,:),J)
       End Do
       ChisqrV = ChisqrV + ((Y(I) - Fval)/(Yerr(I)))**2
    End Do
    
    ChisqrV = ChisqrV / Real(Npoints - Nparm,kind=SP)


    Return
  End Subroutine MultiLinearReg_SP

!  *********************************************
!  *                                           *
  Subroutine Laplace_SP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with Laplace(Rmu, Rb) 
!  * distribution in the intent (out) Real SP 
!  * variable X.
!  *********************************************

    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent (in) :: Rmu, Rb
    Real (kind=SP) :: U
    
    CALL Random_Number(U)
    U = U - 0.5_SP

    X = Rmu - Rb * Sign(Log(1-2.0_SP*Abs(U)), U)

    Return
  End Subroutine  Laplace_SP

!  *********************************************
!  *                                           *
  Subroutine Laplace_DP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with Laplace(Rmu, Rb) 
!  * distribution in the intent (out) Real DP 
!  * variable X.
!  *********************************************

    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent (in) :: Rmu, Rb
    Real (kind=DP) :: U
    
    CALL Random_Number(U)
    U = U - 0.5_DP

    X = Rmu - Rb * Sign(Log(1-2.0_DP*Abs(U)), U)

    Return
  End Subroutine  Laplace_DP

!  *********************************************
!  *                                           *
  Subroutine Laplace2_SP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with Laplace(Rmu, Rb) 
!  * distribution in the intent (out) Real SP 
!  * variable X.
!  *********************************************

    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent (in) :: Rmu, Rb
    Real (kind=SP) :: U(Size(X))
    
    CALL Random_Number(U)
    U = U - 0.5_SP

    X(:) = Rmu - Rb * Sign(Log(1-2.0_SP*Abs(U(:))), U(:))

    Return
  End Subroutine  Laplace2_SP

!  *********************************************
!  *                                           *
  Subroutine Laplace2_DP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with Laplace(Rmu, Rb) 
!  * distribution in the intent (out) Real DP 
!  * variable X.
!  *********************************************

    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent (in) :: Rmu, Rb
    Real (kind=DP) :: U(Size(X))
    
    CALL Random_Number(U)
    U = U - 0.5_DP

    X(:) = Rmu - Rb * Sign(Log(1-2.0_DP*Abs(U(:))), U(:))

    Return
  End Subroutine  Laplace2_DP

!  *********************************************
!  *                                           *
  Subroutine FishTipp_SP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with FishTipp(Rmu, Rb) 
!  * distribution in the intent (out) Real SP 
!  * variable X.
!  *********************************************

    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent (in) :: Rmu, Rb
    Real (kind=SP) :: U
    
    CALL Random_Number(U)

    X = Rmu - Rb * Log(-Log(U))

    Return
  End Subroutine  FishTipp_SP

!  *********************************************
!  *                                           *
  Subroutine FishTipp_DP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with FishTipp(Rmu, Rb) 
!  * distribution in the intent (out) Real DP 
!  * variable X.
!  *********************************************

    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent (in) :: Rmu, Rb
    Real (kind=DP) :: U
    
    CALL Random_Number(U)

    X = Rmu - Rb * Log(-Log(U))

    Return
  End Subroutine  FishTipp_DP

!  *********************************************
!  *                                           *
  Subroutine FishTipp2_SP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with FishTipp(Rmu, Rb) 
!  * distribution in the intent (out) Real SP 
!  * variable X.
!  *********************************************

    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent (in) :: Rmu, Rb
    Real (kind=SP) :: U(Size(X))
    
    CALL Random_Number(U)

    X = Rmu - Rb * Log(-Log(U(:)))

    Return
  End Subroutine  FishTipp2_SP

!  *********************************************
!  *                                           *
  Subroutine FishTipp2_DP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
!  * Returns a number with FishTipp(Rmu, Rb) 
!  * distribution in the intent (out) Real DP 
!  * variable X.
!  *********************************************

    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent (in) :: Rmu, Rb
    Real (kind=DP) :: U(Size(X))
    
    CALL Random_Number(U)

    X = Rmu - Rb * Log(-Log(U(:)))

    Return
  End Subroutine  FishTipp2_DP

! ********************************************
! *
  Function Irand_S(I, J)
! *
! ********************************************
! * Generates a Random integer number between
! * I and J.
! ********************************************

    Integer, Intent (in) :: I, J
    Integer :: Irand_s
    
    Real (kind=DP) :: U
    
    CALL Random_Number(U)
    Irand_S = Int(U*(J-I+1)) + 1

    Return
  End Function Irand_S

! ********************************************
! *
  Subroutine Irand_V(Ir, I, J)
! *
! ********************************************
! * Fills Ir(:) with Random integers numbers
! * between I and J.
! ********************************************

    Integer, Intent (in) :: I, J
    Integer, Intent (out) :: Ir(:)
    
    Real (kind=DP) :: U(Size(Ir))
    
    CALL Random_Number(U)
    Ir = Int(U*(J-I+1)) + 1

    Return
  End Subroutine Irand_V

! ********************************************
! *
  Subroutine Bootstrap(Ibt)
! *
! ********************************************
! * Generates Nb Bootstrap sequence of N numbers 
! * each. Ibt is an integer two dimensional 
! * array of sizes N x Nb
! ********************************************

    Integer, Intent (out) :: Ibt(:,:)

    Integer :: Nb, N, I

    Nb = Size(Ibt,2)
    N  = Size(Ibt,1)

    Do I = 1, N
       CALL Irand(Ibt(I,:), 1, N) 
    End Do

    Return
  End Subroutine Bootstrap

! ********************************************
! *
  Subroutine ReadBstrp(Ibt, Filename)
! *
! ********************************************
! * Reads Nb Bootstrap sequence of N numbers 
! * each. Ibt is an integer two dimensional 
! * array of sizes N x Nb in the file Filename.
! ********************************************

    Integer, Intent (out) :: Ibt(:,:)
    Character (len=*), Intent (in) :: Filename

    Integer :: I, J
    Character (len=40) :: Fmt 

    Open (Unit=55, File = Trim(Filename), ACTION="READ")
    Write(Fmt,'(1A1,1I12,1A4)')'(',Size(Ibt,1),'I12)'
    Do I = 1, Size(Ibt,2)
       Read(55,Fmt)(Ibt(J,I), J=1, Size(Ibt,1)) 
    End Do

    Return
  End Subroutine ReadBstrp

! ********************************************
! *
  Subroutine SaveBstrp(Ibt, Filename)
! *
! ********************************************
! * Saves Nb Bootstrap sequence of N numbers 
! * each. Ibt is an integer two dimensional 
! * array of sizes N x Nb in the file Filename.
! ********************************************
  
    Integer, Intent (in) :: Ibt(:,:)
    Character (len=*), Intent (in) :: Filename

    Integer :: I, J
    Character (len=40) :: Fmt 

    Open (Unit=55, File = Trim(Filename), ACTION="WRITE")
    Write(Fmt,'(1A1,1I12,1A4)')'(',Size(Ibt,1),'I12)'
    Do I = 1, Size(Ibt,2)
       Write(55,Fmt)(Ibt(J,I), J=1, Size(Ibt,1)) 
    End Do

    Return
  End Subroutine SaveBstrp

! ********************************************
! *
  Subroutine EstBstrp(Data, Ibt, Func, Val, Err)
! *
! ********************************************
! * Estimates using the Bottstrap method the 
! * average and error of an estimator given as
! * a user suplied function
! ********************************************
  
    Real (kind=DP), Intent (in) :: Data(:)
    Integer, Intent (in) :: Ibt(:,:)
    Real (kind=DP), Intent (out) :: Val, Err

    Real (kind=DP) :: Rest(Size(Ibt,2))
    Integer :: I, Nb

    Interface 
       Function Func(X)
         USE NumTypes
         
         Real (kind=DP), Intent (in) :: X(:)
         Real (kind=DP) :: Func

       End Function Func
    End Interface

    Val = Func(Data)
    
    Nb = Size(Ibt,2)
    Do I = 1, Nb
       Rest(I) = Func(Data(Ibt(:,I)))     
    End Do
    Err = Stddev(Rest)

    Return
  End Subroutine EstBstrp

! ********************************************
! *
  Subroutine EstBstrp_H(Data, Ibt, Func, Val, Err, Rest)
! *
! ********************************************
! * Estimates using the Bottstrap method the 
! * average and error of an estimator given as
! * a user suplied function. Returns info for writting
! * the bootstrp histogram.
! ********************************************
  
    Real (kind=DP), Intent (in) :: Data(:)
    Integer, Intent (in) :: Ibt(:,:)
    Real (kind=DP), Intent (out) :: Val, Err
    Real (kind=DP), Intent (out) :: Rest(Size(Ibt,2))

    Integer :: I, Nb

    Interface 
       Function Func(X)
         USE NumTypes
         
         Real (kind=DP), Intent (in) :: X(:)
         Real (kind=DP) :: Func

       End Function Func
    End Interface

    Val = Func(Data)
    
    Nb = Size(Ibt,2)
    Do I = 1, Nb
       Rest(I) = Func(Data(Ibt(:,I)))     
    End Do
    Err = Stddev(Rest)

    Return
  End Subroutine EstBstrp_H

End MODULE Statistics
