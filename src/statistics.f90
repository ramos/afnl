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
  USE NonNumeric
  USE Linear

  IMPLICIT NONE

  Interface Mean
     Module Procedure Mean_DP, Mean_SP
  End Interface

  Interface Median
     Module Procedure Median_DP, Median_SP
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

  Interface NonLinearReg
     Module Procedure NonLinearReg_DP, NonLinearReg_SP, &
          & MultiNonLinearReg_DP, MultiNonLinearReg_SP
  End Interface

  Interface LinearReg
     Module Procedure LinearReg_DP, LinearReg_SP, &
          & LinearReg_Pol_DP, LinearReg_Pol_SP, &
          & MultiLinearReg_DP, MultiLinearReg_SP
  End Interface

  Interface Random_Number
     Module Procedure My_Random_Number_SP, My_Random_Number_DP, &
          & My_Random_Number_SP_S, My_Random_Number_DP_S
  End Interface
  
  ! Parameters to tune the Nonlinear fitting routines:
  !
  !  - Rlambda_DEF: Default value used in the Levenberg-Marquardt
  !             algorithm.
  !  - RFACTOR_DEF: Factor of increasing/decreasing Rlambda
  !                 after iterations.
  !  - TOL: Amount that the ChiSqr must. Used with 
  !         Nconv to stop the iteration.
  !  - Nconv: The number of times the Chisqr must decrease
  !           less than TOL to stop the iteration.
  Real (kind=DP) :: TOL = 1.0E-2_DP, Rlambda_DEF = 1.0E-3_DP, &
       & RFACTOR_DEF = 10.0_DP
  Integer :: Nconv = 2, MaxIter = 10000

  ! Parameters for Luxury random number generator
  Integer, Parameter :: LUX_jump(6) = (/0, 24, 73, 194, 380, 770/)
  Integer, Parameter :: LUX_base = 24
  Integer :: LUX_RndLevel = 5 ! Default a cheap and bad 
                          ! random number generator:
                          ! the intrinsic F90 routine
  Real (kind=SP) :: LUX_Seed(LUX_base)
  Real (kind=SP) :: LUX_carry, LUX_lastbit
  Logical :: LUX_FirstCall = .True.
  Integer :: LUX_ir, LUX_is, LUX_next(LUX_base), LUX_Pos
  


  Private NormalS, NormalV, NormalS2, NormalV2, &
       & NormalS_SP, NormalV_SP, NormalS2_SP, NormalV2_SP, &
       & Stddev_DP, Var_DP, Mean_DP, Moment_DP, &
       & Stddev_SP, Var_SP, Mean_SP, Moment_SP, ChiSqr_SP, &
       & ChiSqr_DP, LinearReg_DP, LinearReg_SP, LinearReg_Pol_DP, &
       & LinearReg_Pol_SP, MultiLinearReg_DP, MultiLinearReg_SP, &
       & Laplace_DP, Laplace_SP,Laplace2_DP, Laplace2_SP, &
       & FishTipp_DP, FishTipp_SP,FishTipp2_DP, FishTipp2_SP, &
       & Irand_S, Irand_V, EstBstrp_H, NonLinearReg_DP, TOL, &
       & Rlambda_DEF, Nconv, RFACTOR_DEF, NonLinearReg_SP, &
       & MultiNonLinearReg_SP, Median_SP, Median_DP, LUX_base, &
       & LUX_seed, LUX_carry, LUX_lastbit, LUX_FirstCall, &
       & LUX_RndLevel, LUX_Init, LUX_jump, My_Random_Number_SP, &
       & My_Random_Number_DP, RanLux_SP, RanLux_DP, RanLux_SP_S,&
       & My_Random_Number_SP_S, My_Random_Number_DP_S, RanLux_DP_S

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
  Real (kind=DP) Function Median_DP(X)
!  *                                           *
!  *********************************************
!  * Returns the median value of the elements in 
!  * the vector X(:)
!  *********************************************

    Real (kind=DP), Intent (in) :: X(:)

    Real (kind=DP) :: Xcp(Size(X))
    Integer :: Ns, Nsd2

    Xcp = X
    Ns = Size(X)
    Nsd2 = Int(Ns/2)

    CALL Qsort(Xcp)
    If (Mod(Ns,2) == 1) Then
       Median_DP = Xcp(Nsd2+1)
    Else
       Median_DP = (Xcp(Nsd2) + Xcp(Nsd2+1))/2.0_DP
    End If


    Return
  End Function Median_DP

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
  Real (kind=SP) Function Median_SP(X)
!  *                                           *
!  *********************************************
!  * Returns the median value of the elements in 
!  * the vector X(:)
!  *********************************************

    Real (kind=SP), Intent (in) :: X(:)

    Real (kind=SP) :: Xcp(Size(X))
    Integer :: Ns, Nsd2

    Xcp = X
    Ns = Size(X)
    Nsd2 = Int(Ns/2)

    CALL Qsort(Xcp)
    If (Mod(Ns,2) == 1) Then
       Median_SP = Xcp(Nsd2+1)
    Else
       Median_SP = (Xcp(Nsd2) + Xcp(Nsd2+1))/2.0_SP
    End If


    Return
  End Function Median_SP

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
    Irand_S = Int((J-I+1)*U + I)

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
    Ir = Int((J-I+1)*U + I)

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
       CALL Irand_V(Ibt(I,:), 1, N) 
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

! ********************************************
! *
  Subroutine BstrpConfInt(Data, Ibt, alpha, Func, dmin, dpls)
! *
! ********************************************
! * Estimates using the Bootstrap method a confidence
! * interval of alpha for the the estimator
! * given as a user suplied function. Returns info 
! * for writting the bootstrp histogram. Returns 
! * interval as (dmin, dpls).
! ********************************************

    Real (kind=DP), Intent (in) :: Data(:), alpha
    Integer, Intent (in) :: Ibt(:,:)
    Real (kind=DP), Intent (out) :: dmin, dpls

    Real (kind=DP) :: Rest(Size(Ibt,2))
    Integer :: Nb, Nm, Np, I

    Interface 
       Function Func(X)
         USE NumTypes
         
         Real (kind=DP), Intent (in) :: X(:)
         Real (kind=DP) :: Func

       End Function Func
    End Interface
    
    Nb = Size(Ibt,2)
    Do I = 1, Nb
       Rest(I) = Func(Data(Ibt(:,I)))     
    End Do
    CALL Qsort(Rest)
    
    Nm = Int(alpha*Nb)!/200.0_DP
    dmin = (Rest(Nm-1) + Rest(Nm) + Rest(Nm+1))/3.0_DP
    Np = Nb - Nm
    dpls = (Rest(Np-1) + Rest(Np) + Rest(Np+1))/3.0_DP

    Return
  End Subroutine BstrpConfInt

!  *********************************************
!  *                                           *
  Subroutine NonLinearReg_DP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a function given by the routine
!  * Func that also returns the first derivatives. 
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr per degree of freedom is 
!  * returned. This routine uses the Levenberg-Marquardt 
!  * algorithm
!  *********************************************

    Real (kind=DP), Intent(in) :: X(:), Y(:), Yerr(:)
    Real (kind=DP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=DP) :: alphap(Size(Coef), Size(Coef)), Rlambda, &
         & Vf, Vd(Size(Coef)), Dcoef(Size(Coef)), Cold, Cnew, Cdif

    Integer :: Npoints, Npar, I, J, K, Kont = 0, Ipiv(Size(Coef)),&
         & Idet, Ktot = 0

    Logical :: Cnt = .True.

    Interface
       Subroutine Func(X, Cf, Valf, ValD)
         
         USE NumTypes

         Real (kind=DP), Intent (in) :: X, Cf(:)
         Real (kind=DP), Intent (out) :: Valf, ValD(Size(Cf))
         
       End Subroutine Func
    End Interface

    Npoints = Size(X)
    Npar = Size(Coef)
    Rlambda = Rlambda_DEF
 
    Ktot = 0
    Cnt = .True.
    Do While (Cnt)
       Ktot = Ktot + 1
       alphap = 0.0_DP
       DCoef = 0.0_DP
       Cold = CalcCh2()
       Do K = 1, Npoints
          CALL Func(X(K), Coef, Vf, VD)
          Do I = 1, Npar          
             DCoef(I) = DCoef(I) + (Y(K) - Vf)/Yerr(K)**2 * VD(I)
             Do J = 1, Npar
                alphap(I,J) = alphap(I,J) + &
                     & ( VD(I)*VD(J) )/Yerr(K)**2
             End Do
          End Do
       End Do
       
       Do I = 1, Npar
          alphap(I,I) = alphap(I,I)*(Rlambda + 1.0_DP)
       End Do

       Cnew = 0.0_DP
       
       CALL LuSolve(alphap, DCoef)
       Coef = Coef + DCoef
       Cnew = CalcCh2()
       Cdif = CNew - Cold
       If (Cdif > 0.0_DP) Then
          RLambda = RFACTOR_DEF*RLambda
          Coef = Coef - DCoef
          kont = 0
       Else
          RLambda = RLambda/RFACTOR_DEF
          If (Abs(Cdif) < TOL) Then
             kont = kont + 1
          Else 
             kont = 0
          End If
       End If
       
       If ( (kont >= Nconv).or.(Ktot > MaxIter) ) Then
!          Write(*,*)'Dentro: ', Kont, Ktot
          ChiSqrV = Cnew/Real(Npoints - Npar, kind=DP)
          Do K = 1, Npoints
             CALL Func(X(K), Coef, Vf, VD)
             Do I = 1, Npar          
                Do J = 1, Npar
                   alphap(I,J) = alphap(I,J) + &
                        & ( VD(I)*VD(J) )/Yerr(K)**2
                End Do
             End Do
          End Do
          
          CALL LU(alphap, Ipiv, Idet)
          Do J = 1, Npar
             Do I = 1, Npar
                If (Ipiv(I) == J) Then
                   DCoef(I) = 1.0_DP
                Else
                   DCoef(I) = 0.0_DP
                End If
             End Do
             
             Do I = 2, Npar
                DCoef(I) = DCoef(I) - & 
                     & Dot_Product(alphap(I,1:I-1),DCoef(1:I-1))
             End Do
             DCoef(Npar) = DCoef(Npar) / alphap(Npar, Npar)
             Do I = Npar - 1, J, -1
                DCoef(I) = DCoef(I) - &
                     & Dot_Product(alphap(I, I+1:Npar),DCoef(I+1:Npar))
                DCoef(I) = DCoef(I) / alphap(I, I)
             End Do
             Cerr(J) = Sqrt(Abs(DCoef(J)))
          End Do
          
          Cnt = .False.
       End If
    End Do

    Return

  CONTAINS

    Function CalcCh2()

      Real (kind=DP) :: VVF, VVD(Size(Coef))
      Real (kind=DP) :: CalcCh2

      Integer :: Ii

      CalcCh2 = 0.0_DP
      Do Ii = 1, Npoints
         CALL Func(X(Ii), Coef, VVf, VVd)
         CalcCh2 = CalcCh2 + ( (Y(Ii) - VVf)/Yerr(Ii) )**2
      End Do

      Return
    End Function CalcCh2

  End Subroutine NonLinearReg_DP

!  *********************************************
!  *                                           *
  Subroutine NonLinearReg_SP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a function given by the routine
!  * Func that also returns the first derivatives. 
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr per degree of freedom is 
!  * returned. This routine uses the Levenberg-Marquardt 
!  * algorithm
!  *********************************************

    Real (kind=SP), Intent(in) :: X(:), Y(:), Yerr(:)
    Real (kind=SP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=SP) :: alphap(Size(Coef), Size(Coef)), Rlambda, &
         & Vf, Vd(Size(Coef)), Dcoef(Size(Coef)), Cold, Cnew, Cdif

    Integer :: Npoints, Npar, I, J, K, Kont = 0, Ipiv(Size(Coef)),&
         & Idet, Ktot = 0

    Logical :: Cnt = .True.

    Interface
       Subroutine Func(X, Cf, Valf, ValD)
         
         USE NumTypes

         Real (kind=SP), Intent (in) :: X, Cf(:)
         Real (kind=SP), Intent (out) :: Valf, ValD(Size(Cf))
         
       End Subroutine Func
    End Interface

    Npoints = Size(X)
    Npar = Size(Coef)
    Rlambda = Rlambda_DEF
    
    Ktot = 0
    Cnt = .True.
    Do While (Cnt)
       Ktot = Ktot + 1
       alphap = 0.0_SP
       DCoef = 0.0_SP
       Cold = CalcCh2()
       Do K = 1, Npoints
          CALL Func(X(K), Coef, Vf, VD)
          Do I = 1, Npar          
             DCoef(I) = DCoef(I) + (Y(K) - Vf)/Yerr(K)**2 * VD(I)
             Do J = 1, Npar
                alphap(I,J) = alphap(I,J) + &
                     & ( VD(I)*VD(J) )/Yerr(K)**2
             End Do
          End Do
       End Do
       
       Do I = 1, Npar
          alphap(I,I) = alphap(I,I)*(Rlambda + 1.0_SP)
       End Do

       Cnew = 0.0_SP
       
       CALL LuSolve(alphap, DCoef)
       Coef = Coef + DCoef
       Cnew = CalcCh2()
       Cdif = CNew - Cold
       If (Cdif > 0.0_SP) Then
          RLambda = RFACTOR_DEF*RLambda
          Coef = Coef - DCoef
          kont = 0
       Else
          RLambda = RLambda/RFACTOR_DEF
          If (Abs(Cdif) < TOL) Then
             kont = kont + 1
          Else 
             kont = 0
          End If
       End If
       
       If ( (kont >= Nconv).or.(Ktot > MaxIter) ) Then 
          ChiSqrV = Cnew/Real(Npoints - Npar, kind=SP)
          Do K = 1, Npoints
             CALL Func(X(K), Coef, Vf, VD)
             Do I = 1, Npar          
                Do J = 1, Npar
                   alphap(I,J) = alphap(I,J) + &
                        & ( VD(I)*VD(J) )/Yerr(K)**2
                End Do
             End Do
          End Do
          
          CALL LU(alphap, Ipiv, Idet)
          Do J = 1, Npar
             Do I = 1, Npar
                If (Ipiv(I) == J) Then
                   DCoef(I) = 1.0_SP
                Else
                   DCoef(I) = 0.0_SP
                End If
             End Do
             
             Do I = 2, Npar
                DCoef(I) = DCoef(I) - & 
                     & Dot_Product(alphap(I,1:I-1),DCoef(1:I-1))
             End Do
             DCoef(Npar) = DCoef(Npar) / alphap(Npar, Npar)
             Do I = Npar - 1, J, -1
                DCoef(I) = DCoef(I) - &
                     & Dot_Product(alphap(I, I+1:Npar),DCoef(I+1:Npar))
                DCoef(I) = DCoef(I) / alphap(I, I)
             End Do
             Cerr(J) = Sqrt(Abs(DCoef(J)))
          End Do
          
          Cnt = .False.
       End If
    End Do

    Return

  CONTAINS

    Function CalcCh2()

      Real (kind=SP) :: VVF, VVD(Size(Coef))
      Real (kind=SP) :: CalcCh2

      Integer :: Ii

      CalcCh2 = 0.0_SP
      Do Ii = 1, Npoints
         CALL Func(X(Ii), Coef, VVf, VVd)
         CalcCh2 = CalcCh2 + ( (Y(Ii) - VVf)/Yerr(Ii) )**2
      End Do

      Return
    End Function CalcCh2

  End Subroutine NonLinearReg_SP

!  *********************************************
!  *                                           *
  Subroutine MultiNonLinearReg_DP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a function given by the routine
!  * Func that also returns the first derivatives. 
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr per degree of freedom is 
!  * returned. This routine uses the Levenberg-Marquardt 
!  * algorithm
!  *********************************************

    Real (kind=DP), Intent(in) :: X(:,:), Y(:), Yerr(:)
    Real (kind=DP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=DP) :: alphap(Size(Coef), Size(Coef)), Rlambda, &
         & Vf, Vd(Size(Coef)), Dcoef(Size(Coef)), Cold, Cnew, Cdif

    Integer :: Npoints, Npar, I, J, K, Kont = 0, Ipiv(Size(Coef)),&
         & Idet, Ndim, Ktot = 0

    Logical :: Cnt = .True.

    Interface
       Subroutine Func(X, Cf, Valf, ValD)
         
         USE NumTypes

         Real (kind=DP), Intent (in) :: X(:), Cf(:)
         Real (kind=DP), Intent (out) :: Valf, ValD(Size(Cf))
         
       End Subroutine Func
    End Interface

    Npoints = Size(X,1)
    Ndim = Size(X,2)
    Npar = Size(Coef)
    Rlambda = Rlambda_DEF
    
    Ktot = 0
    Cnt = .True.
    Do While (Cnt)
       Ktot = Ktot + 1
       alphap = 0.0_DP
       DCoef = 0.0_DP
       Cold = CalcCh2()
       Do K = 1, Npoints
          CALL Func(X(K,:), Coef, Vf, VD)
          Do I = 1, Npar          
             DCoef(I) = DCoef(I) + (Y(K) - Vf)/Yerr(K)**2 * VD(I)
             Do J = 1, Npar
                alphap(I,J) = alphap(I,J) + &
                     & ( VD(I)*VD(J) )/Yerr(K)**2
             End Do
          End Do
       End Do
       
       Do I = 1, Npar
          alphap(I,I) = alphap(I,I)*(Rlambda + 1.0_DP)
       End Do

       Cnew = 0.0_DP
       
       CALL LuSolve(alphap, DCoef)
       Coef = Coef + DCoef
       Cnew = CalcCh2()
       Cdif = CNew - Cold
       If (Cdif > 0.0_DP) Then
          RLambda = RFACTOR_DEF*RLambda
          Coef = Coef - DCoef
          kont = 0
       Else
          RLambda = RLambda/RFACTOR_DEF
          If (Abs(Cdif) < TOL) Then
             kont = kont + 1
          Else 
             kont = 0
          End If
       End If
       
       If ( (kont >= Nconv).or.(Ktot > MaxIter) ) Then 
          ChiSqrV = Cnew/Real(Npoints - Npar, kind=DP)
          Do K = 1, Npoints
             CALL Func(X(K,:), Coef, Vf, VD)
             Do I = 1, Npar          
                Do J = 1, Npar
                   alphap(I,J) = alphap(I,J) + &
                        & ( VD(I)*VD(J) )/Yerr(K)**2
                End Do
             End Do
          End Do
          
          CALL LU(alphap, Ipiv, Idet)
          Do J = 1, Npar
             Do I = 1, Npar
                If (Ipiv(I) == J) Then
                   DCoef(I) = 1.0_DP
                Else
                   DCoef(I) = 0.0_DP
                End If
             End Do
             
             Do I = 2, Npar
                DCoef(I) = DCoef(I) - & 
                     & Dot_Product(alphap(I,1:I-1),DCoef(1:I-1))
             End Do
             DCoef(Npar) = DCoef(Npar) / alphap(Npar, Npar)
             Do I = Npar - 1, J, -1
                DCoef(I) = DCoef(I) - &
                     & Dot_Product(alphap(I, I+1:Npar),DCoef(I+1:Npar))
                DCoef(I) = DCoef(I) / alphap(I, I)
             End Do
             Cerr(J) = Sqrt(Abs(DCoef(J)))
          End Do
          
          Cnt = .False.
       End If
    End Do

    Return

  CONTAINS

    Function CalcCh2()

      Real (kind=DP) :: VVF, VVD(Size(Coef))
      Real (kind=DP) :: CalcCh2

      Integer :: Ii

      CalcCh2 = 0.0_DP
      Do Ii = 1, Npoints
         CALL Func(X(Ii,:), Coef, VVf, VVd)
         CalcCh2 = CalcCh2 + ( (Y(Ii) - VVf)/Yerr(Ii) )**2
      End Do

      Return
    End Function CalcCh2

  End Subroutine MultiNonLinearReg_DP


!  *********************************************
!  *                                           *
  Subroutine MultiNonLinearReg_SP(X, Y, Yerr, Func, Coef, Cerr, ChisqrV)
!  *                                           *
!  *********************************************
!  * Given a set of points (X(:), Y(:)), this routine 
!  * fit the points to a function given by the routine
!  * Func that also returns the first derivatives. 
!  * The errors in the coefficients are returned in 
!  * Cerr(:), and the ChiSqr per degree of freedom is 
!  * returned. This routine uses the Levenberg-Marquardt 
!  * algorithm
!  *********************************************

    Real (kind=SP), Intent(in) :: X(:,:), Y(:), Yerr(:)
    Real (kind=SP), Intent(out) :: Coef(:), Cerr(:), ChisqrV

    Real (kind=SP) :: alphap(Size(Coef), Size(Coef)), Rlambda, &
         & Vf, Vd(Size(Coef)), Dcoef(Size(Coef)), Cold, Cnew, Cdif

    Integer :: Npoints, Npar, I, J, K, Kont = 0, Ipiv(Size(Coef)),&
         & Idet, Ndim, Ktot = 0

    Logical :: Cnt = .True.

    Interface
       Subroutine Func(X, Cf, Valf, ValD)
         
         USE NumTypes

         Real (kind=SP), Intent (in) :: X(:), Cf(:)
         Real (kind=SP), Intent (out) :: Valf, ValD(Size(Cf))
         
       End Subroutine Func
    End Interface

    Npoints = Size(X,1)
    Ndim = Size(X,2)
    Npar = Size(Coef)
    Rlambda = Rlambda_DEF
    
    Ktot = 0
    Cnt = .True.
    Do While (Cnt)
       Ktot = Ktot + 1
       alphap = 0.0_SP
       DCoef = 0.0_SP
       Cold = CalcCh2()
       Do K = 1, Npoints
          CALL Func(X(K,:), Coef, Vf, VD)
          Do I = 1, Npar          
             DCoef(I) = DCoef(I) + (Y(K) - Vf)/Yerr(K)**2 * VD(I)
             Do J = 1, Npar
                alphap(I,J) = alphap(I,J) + &
                     & ( VD(I)*VD(J) )/Yerr(K)**2
             End Do
          End Do
       End Do
       
       Do I = 1, Npar
          alphap(I,I) = alphap(I,I)*(Rlambda + 1.0_SP)
       End Do

       Cnew = 0.0_SP
       
       CALL LuSolve(alphap, DCoef)
       Coef = Coef + DCoef
       Cnew = CalcCh2()
       Cdif = CNew - Cold
       If (Cdif > 0.0_SP) Then
          RLambda = RFACTOR_DEF*RLambda
          Coef = Coef - DCoef
          kont = 0
       Else
          RLambda = RLambda/RFACTOR_DEF
          If (Abs(Cdif) < TOL) Then
             kont = kont + 1
          Else 
             kont = 0
          End If
       End If
       
       If ( (kont >= Nconv).or.(Ktot > MaxIter) ) Then 
          ChiSqrV = Cnew/Real(Npoints - Npar, kind=SP)
          Do K = 1, Npoints
             CALL Func(X(K,:), Coef, Vf, VD)
             Do I = 1, Npar          
                Do J = 1, Npar
                   alphap(I,J) = alphap(I,J) + &
                        & ( VD(I)*VD(J) )/Yerr(K)**2
                End Do
             End Do
          End Do
          
          CALL LU(alphap, Ipiv, Idet)
          Do J = 1, Npar
             Do I = 1, Npar
                If (Ipiv(I) == J) Then
                   DCoef(I) = 1.0_SP
                Else
                   DCoef(I) = 0.0_SP
                End If
             End Do
             
             Do I = 2, Npar
                DCoef(I) = DCoef(I) - & 
                     & Dot_Product(alphap(I,1:I-1),DCoef(1:I-1))
             End Do
             DCoef(Npar) = DCoef(Npar) / alphap(Npar, Npar)
             Do I = Npar - 1, J, -1
                DCoef(I) = DCoef(I) - &
                     & Dot_Product(alphap(I, I+1:Npar),DCoef(I+1:Npar))
                DCoef(I) = DCoef(I) / alphap(I, I)
             End Do
             Cerr(J) = Sqrt(Abs(DCoef(J)))
          End Do
          
          Cnt = .False.
       End If
    End Do

    Return

  CONTAINS

    Function CalcCh2()

      Real (kind=SP) :: VVF, VVD(Size(Coef))
      Real (kind=SP) :: CalcCh2

      Integer :: Ii

      CalcCh2 = 0.0_SP
      Do Ii = 1, Npoints
         CALL Func(X(Ii,:), Coef, VVf, VVd)
         CalcCh2 = CalcCh2 + ( (Y(Ii) - VVf)/Yerr(Ii) )**2
      End Do

      Return
    End Function CalcCh2

  End Subroutine MultiNonLinearReg_SP

! ********************************************
! *
  Subroutine Permutation(Idx)
! *
! ********************************************
! * Generates a random permutation.
! ********************************************

    Integer, Intent (out) :: Idx(:)

    Integer :: I, J, Ns 

    Ns = Size(Idx)

    Forall (I=1:Ns) Idx(I) = I
    Do I = 1, Ns
       J = Irand(1, Ns)
       If (I /= J) CALL Swap(Idx, I, J)
    End Do

    Return
  End Subroutine Permutation

! ********************************************
! *
  Subroutine PutLuxSeed(ISeed)
! *
! ********************************************
! * 
! ********************************************

    Integer, Intent (in), Optional :: ISeed(LUX_base+1)

    Integer :: I, jseed, it24, k 


    If (LUX_FirstCall) Then
       LUX_FirstCall = .False.
       CALL LUX_Init()
    End If
       
    If (Present(Iseed)) Then
       If (ISeed(LUX_base+1) < 0) Then
          LUX_carry = LUX_lastbit
       Else
          LUX_carry = 0.0_SP
       End If
       LUX_ir = Mod(Abs(ISeed(LUX_base+1)), 100)
       LUX_is = Mod(Int(Abs(ISeed(LUX_base+1))/100), 100)
       LUX_RndLevel = Mod(Int(Abs(ISeed(LUX_base+1))/10000), 100)
       LUX_Pos = Mod(Int(Abs(ISeed(LUX_base+1))/1000000), 100)
       
       LUX_Seed(1:LUX_base) = Real(ISeed(1:LUX_base),kind=SP) * LUX_lastbit
    Else
       ! Switch to default initialization. Try to make it equal 
       ! as in the reference code by Allan miller.
       CALL LUX_Init()
       jseed = 314159265
       it24 = 2**24
       Do I = 1, LUX_base
          k = Int(jseed/53668)
          jseed = 40014 * (jseed-k*53668) - k * 12211
          If (jseed < 0) jseed = jseed + 2147483563
          LUX_Seed(I) = Real(Mod(jseed,it24),kind=SP) * LUX_lastbit
       End Do
    End If

    Return
  End Subroutine PutLuxSeed

! ********************************************
! *
  Subroutine GetLuxSeed(ISeed)
! *
! ********************************************
! * Output the seed, and the position on the 
! * generator as a vector of 25 integers.
! ********************************************

    Integer, Intent (out) :: ISeed(LUX_base+1)
     
    If (LUX_firstcall) Then
       CALL Perror('GetLuxSeed: ', 'Seed of the luxury Random Number generator not initialised.')
       Return
    End If

    ISeed(1:LUX_base) = Int(Lux_Seed(1:LUX_base)*Scale(1.0_SP, 12)**2)
    
    ISeed(LUX_base+1) = LUX_ir + 100*LUX_is + 10000*LUX_RndLevel + &
         & 1000000*LUX_Pos
    If (LUX_carry > 0.0_SP) ISeed(LUX_base+1) = -ISeed(LUX_base+1) 

    Return
  End Subroutine GetLuxSeed

! ********************************************
! *
  Subroutine RanLux_SP(Rnd)
! *
! ********************************************
! * 
! ********************************************
    

  Real (kind=SP), Intent (out) :: Rnd(:)

  Real (kind=SP) :: Delta
  Integer :: RndLen, I, J, Njump  

  If (LUX_FirstCall) Then
     LUX_FirstCall = .False.
     CALL LUX_Init()
     CALL PutLuxSeed()
  End If
  
  Rndlen = Size(Rnd, 1)
  If (LUX_RndLevel <= 6) Then
     Njump = LUX_jump(LUX_RndLevel)
  Else If (LUX_RndLevel < 24) Then
     CALL Perror('RanLux_SP', 'This Random Level has no sense for LUXURY. Swithching to default Random Level')
     Njump = LUX_jump(5)
  Else
     Njump = LUX_RndLevel - 24
  End If

  I = 0
  Do While (.True.)
     Do J = LUX_Pos, 24
        Delta = LUX_Seed(LUX_is) - LUX_Seed(LUX_ir) - LUX_carry
        If (Delta < 0.0_SP) Then
           Delta = Delta + 1.0_SP
           LUX_carry = LUX_lastbit
        Else
           LUX_carry = 0.0_SP
        End If
        LUX_Seed(LUX_ir) = Delta
        LUX_ir = LUX_next(LUX_ir)
        LUX_is = LUX_next(LUX_is)
        I = I + 1
        Rnd(I) = Delta
        If (I == RndLen) Then
           LUX_Pos = J + 1
           Return
        End If
     End Do
     LUX_Pos = 1

     Do J = 1, Njump
        Delta = LUX_Seed(LUX_is) - LUX_Seed(LUX_ir) - LUX_carry
        If (Delta < 0.0_SP) Then
           Delta = Delta + 1.0_SP
           LUX_carry = LUX_lastbit
        Else
           LUX_carry = 0.0_SP
        End If
        LUX_Seed(LUX_ir) = Delta
        LUX_ir = LUX_next(LUX_ir)
        LUX_is = LUX_next(LUX_is)
     End Do
  End Do

  Return
End Subroutine RanLux_SP


! ********************************************
! *
Subroutine RanLux_SP_S(Rnd)
! *
! ********************************************
! * 
! ********************************************
    
  Real (kind=SP), Intent (out) :: Rnd
  
  Real (kind=SP) :: R(1)

  CALL RanLux_SP(R)
  Rnd = R(1)

  Return
End Subroutine RanLux_SP_S

! ***********************************
! *
Subroutine LUX_init()
! *
! ************************************

  Integer :: I

  Forall (I=2:LUX_base) LUX_next(I) = I-1
  LUX_next(1) = LUX_base
  LUX_ir = LUX_base
  LUX_is = 10
  
  LUX_lastbit = Scale(1.0_SP, -LUX_base) ! = 2^{-24}
  LUX_carry = 0.0_SP
  LUX_Pos = 1

  Return
End Subroutine LUX_init

! ********************************************
! *
  Subroutine RanLux_DP(Rnd)
! *
! ********************************************
! * 
! ********************************************
    

  Real (kind=DP), Intent (out) :: Rnd(:)

  Real (kind=SP) :: Raux(2*Size(Rnd))
  Integer :: I


  CALL RanLux_SP(Raux)
  Do I = 1, Size(Rnd)
     Rnd(I) = Real(Raux(2*I-1),kind=DP) + Real(Raux(2*I))*LUX_lastbit
  End Do

  Return
End Subroutine RanLux_DP

! ********************************************
! *
  Subroutine RanLux_DP_S(Rnd)
! *
! ********************************************
! * 
! ********************************************
    

  Real (kind=DP), Intent (out) :: Rnd

  Real (kind=SP) :: Raux(2)


  CALL RanLux_SP(Raux)
  Rnd = Real(Raux(1),kind=DP) + Real(Raux(2))*LUX_lastbit


  Return
End Subroutine RanLux_DP_S

! ********************************************
! *
  Subroutine SetLuxLevel(Ilevel)
! *
! ********************************************
! * 
! ********************************************

    Integer, Intent (in) :: Ilevel

    LUX_RndLevel = Ilevel

    Return
  End Subroutine SetLuxLevel

! ********************************************
! *
  Subroutine My_Random_Number_SP(Rnd)
! *
! ********************************************
! * 
! ********************************************

    Real (kind=SP), Intent (out) :: Rnd(:)

    If (LUX_RndLevel == 0) Then
       CALL Random_Number(Rnd)
    Else
       CALL RanLux_SP(Rnd)
    End If
       
    Return
  End Subroutine My_Random_Number_SP


! ********************************************
! *
  Subroutine My_Random_Number_DP(Rnd)
! *
! ********************************************
! * 
! ********************************************

    Real (kind=DP), Intent (out) :: Rnd(:)

    If (LUX_RndLevel == 0) Then
       CALL Random_Number(Rnd)
    Else
       CALL RanLux_DP(Rnd)
    End If
       
    Return
  End Subroutine My_Random_Number_DP

! ********************************************
! *
  Subroutine My_Random_Number_SP_S(Rnd)
! *
! ********************************************
! * 
! ********************************************

    Real (kind=SP), Intent (out) :: Rnd

    If (LUX_RndLevel == 0) Then
       CALL Random_Number(Rnd)
    Else
       CALL RanLux_SP_S(Rnd)
    End If
       
    Return
  End Subroutine My_Random_Number_SP_S

! ********************************************
! *
  Subroutine My_Random_Number_DP_S(Rnd)
! *
! ********************************************
! * 
! ********************************************

    Real (kind=DP), Intent (out) :: Rnd

    If (LUX_RndLevel == 0) Then
       CALL Random_Number(Rnd)
    Else
       CALL RanLux_DP_S(Rnd)
    End If
       
    Return
  End Subroutine My_Random_Number_DP_S

End MODULE Statistics
