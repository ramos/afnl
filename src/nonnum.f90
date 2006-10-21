!
! A MODULE for non numerical routines (sorting and locate)
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

! ***************************************************
! *
MODULE NonNumeric
! *
! ***************************************************
! *
! * NonNumeric routines for sorting and locating
! * data.
! *
! ***************************************************


  USE NumTypes
  USE Error


  Interface Qsort
     Module Procedure Qsort_IN, Qsort_SP, Qsort_DP
  End Interface

  Interface Locate
     Module Procedure Locate_IN, Locate_SP, Locate_DP
  End Interface


  Private Qsort_IN, Qsort_SP, Qsort_DP, Locate_IN, Locate_SP, &
       & Locate_DP

CONTAINS


! ***********************************
! *
  Integer Function Locate_IN(X, Xo, Iin)
! *
! ***********************************
! * From an ascendent ordered set 
! * of points X(:), locate will return
! * an integer K, such that 
! * X(K) < Xo < X(K+1). If 0 is returned, 
! * then Xo is smaller than all X(:). If 
! * Size(X) is returned then is bigger
! * than all X(:). Iin (optional) is an
! * initial guess of the position.
! * Integer version.
! ***********************************

    Integer, Intent (in) :: X(:), Xo
    Integer, Intent (in), Optional :: Iin

    Integer :: Nmin, Nmax, Nnew


    ISize = Size(X)

    If (Xo < X(1)) Then
       Locate_IN = 0
       Return
    Else If (Xo > X(Isize)) Then
       Locate_IN = Isize
       Return
    End If

    Nmin = 0
    Nmax = ISize
    Do While (Nmax - Nmin > 1)
       Nnew = Int((Nmax + Nmin)/2)
       If (X(Nnew) < Xo) Then
          Nmin = Nnew
       Else
          Nmax = Nnew
       End If
    End Do

    Locate_IN = Nmin

    Return
  End Function Locate_IN

! ***********************************
! *
  Integer Function Locate_SP(X, Xo, Iin)
! *
! ***********************************
! * From an ascendent ordered set 
! * of points X(:), locate will return
! * an integer K, such that 
! * X(K) < Xo < X(K+1). If 0 is returned, 
! * then Xo is smaller than all X(:). If 
! * Size(X) is returned then is bigger
! * than all X(:). Integer version.
! ***********************************

    Real (kind=SP), Intent (in) :: X(:), Xo
    Integer, Intent (in), Optional :: Iin

    Integer :: Nmin, Nmax, Nnew


    ISize = Size(X)

    If (Xo < X(1)) Then
       Locate_SP = 0
       Return
    Else If (Xo > X(Isize)) Then
       Locate_SP = Isize
       Return
    End If

    Nmin = 0
    Nmax = ISize
    Do While (Nmax - Nmin > 1)
       Nnew = Int((Nmax + Nmin)/2)
       If (X(Nnew) < Xo) Then
          Nmin = Nnew
       Else
          Nmax = Nnew
       End If
    End Do

    Locate_SP = Nmin

    Return
  End Function Locate_SP

! ***********************************
! *
  Integer Function Locate_DP(X, Xo, Iin)
! *
! ***********************************
! * From an ascendent ordered set 
! * of points X(:), locate will return
! * an integer K, such that 
! * X(K) < Xo < X(K+1). If 0 is returned, 
! * then Xo is smaller than all X(:). If 
! * Size(X) is returned then is bigger
! * than all X(:). Integer version.
! ***********************************

    Real (kind=DP), Intent (in) :: X(:), Xo
    Integer, Intent (in), Optional :: Iin

    Integer :: Nmin, Nmax, Nnew


    ISize = Size(X)

    If (Xo < X(1)) Then
       Locate_DP = 0
       Return
    Else If (Xo > X(Isize)) Then
       Locate_DP = Isize
       Return
    End If

    Nmin = 0
    Nmax = ISize
    Do While (Nmax - Nmin > 1)
       Nnew = Int((Nmax + Nmin)/2)
       If (X(Nnew) < Xo) Then
          Nmin = Nnew
       Else
          Nmax = Nnew
       End If
    End Do

    Locate_DP = Nmin

    Return
  End Function Locate_DP


! ***********************************
! *
  Recursive Subroutine Qsort_IN(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * (default, or Idir=+1) or descendent 
! * order (Idir=-1). If present Ipt, a
! * pointer with the changes is returned
! * in Ipt. Integer version.
! ***********************************

    Integer, Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Integer :: Xcp(Size(X))
    Integer :: Temp, Ip(Size(X)), Ipp(Size(X)), Itemp(Size(X))
    


    Isize = Size(X)
    Xcp = X

    If (Isize == 1) Then
       Return
    Else If (Isize == 2) Then
       If (X(1) < X(2)) Then
          Return
       Else
          Temp = X(1)
          X(1) = X(2)
          X(2) = Temp
          Return
       End If
    Else
       Temp = X(1)
       Kmin = 1
       Kmax = Isize
       Do I = 2, Isize
          If (X(I) < Temp) Then
             Ip(Kmin) = I
             Kmin = Kmin + 1
          Else
             Ip(Kmax) = I
             Kmax = Kmax - 1
          End If
       End Do

       Ip(Kmin) = 1
       Do I = 1, Isize
          X(I) = Xcp(Ip(I))
       End Do

       If (Kmin > 2) CALL &
            & Qsort(X(1:Kmin-1), Ip(1:Kmin-1))
       If (Kmin < Isize-1) CALL &
            & Qsort(X(Kmin+1:Isize), Ip(Kmin+1:Isize))
    End If

    If (Present(Ipt)) Ipt = Ip
    
    Return
  End Subroutine Qsort_IN

! ***********************************
! *
  Recursive Subroutine Qsort_SP(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * (default, or Idir=+1) or descendent 
! * order (Idir=-1). If present Ipt, a
! * pointer with the changes is returned
! * in Ipt. Integer version.
! ***********************************

    Real (kind=SP), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Real (kind=SP) :: Xcp(Size(X)), Temp
    Integer :: Ip(Size(X)), Ipp(Size(X)), Itemp(Size(X))
    

    Isize = Size(X)
    Xcp = X

    Forall (I=1:Isize) Ip(I) = I
    If (Isize == 1) Then
       If (Present(Ipt)) Ipt = Ip
       Return
    Else If (Isize == 2) Then
       If (X(1) < X(2)) Then
          If (Present(Ipt)) Ipt = Ip
          Return
       Else
          Temp = X(1)
          X(1) = X(2)
          X(2) = Temp
          Forall (I=1:2) Ip(I) = 3 - I
          If (Present(Ipt)) Ipt = Ip
          Return
       End If
    Else
       Temp = X(1)
       Kmin = 1
       Kmax = Isize
       Do I = 2, Isize
          If (X(I) < Temp) Then
             Ip(Kmin) = I
             Kmin = Kmin + 1
          Else
             Ip(Kmax) = I
             Kmax = Kmax - 1
          End If
       End Do

       Ip(Kmin) = 1
       Do I = 1, Isize
          X(I) = Xcp(Ip(I))
       End Do

       If (Kmin > 2) Then
          CALL  Qsort(X(1:Kmin-1), Ipp(1:Kmin-1))
          Itemp = Ip
          Do I = 1, Kmin-1
             Ip(I) = Itemp(Ipp(I))
          End Do
       End If

       If (Kmin < Isize-1) Then
          CALL Qsort(X(Kmin+1:Isize), Ipp(Kmin+1:Isize))
          Itemp = Ip
          Do I = Kmin+1, Isize
             Ip(I) = Itemp(Ipp(I)+Kmin)
          End Do
       End If
    End If

    If (Present(Ipt)) Ipt = Ip
    
    Return
  End Subroutine Qsort_SP

! ***********************************
! *
  Recursive Subroutine Qsort_DP(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * (default, or Idir=+1) or descendent 
! * order (Idir=-1). If present Ipt, a
! * pointer with the changes is returned
! * in Ipt. Integer version.
! ***********************************

    Real (kind=DP), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Real (kind=DP) :: Xcp(Size(X)), Temp
    Integer :: Ip(Size(X)), Ipp(Size(X)), Itemp(Size(X))
    

    Isize = Size(X)
    Xcp = X

    Forall (I=1:Isize) Ip(I) = I
    If (Isize == 1) Then
       If (Present(Ipt)) Ipt = Ip
       Return
    Else If (Isize == 2) Then
       If (X(1) < X(2)) Then
          If (Present(Ipt)) Ipt = Ip
          Return
       Else
          Temp = X(1)
          X(1) = X(2)
          X(2) = Temp
          Forall (I=1:2) Ip(I) = 3 - I
          If (Present(Ipt)) Ipt = Ip
          Return
       End If
    Else
       Temp = X(1)
       Kmin = 1
       Kmax = Isize
       Do I = 2, Isize
          If (X(I) < Temp) Then
             Ip(Kmin) = I
             Kmin = Kmin + 1
          Else
             Ip(Kmax) = I
             Kmax = Kmax - 1
          End If
       End Do

       Ip(Kmin) = 1
       Do I = 1, Isize
          X(I) = Xcp(Ip(I))
       End Do

       If (Kmin > 2) Then
          CALL  Qsort(X(1:Kmin-1), Ipp(1:Kmin-1))
          Itemp = Ip
          Do I = 1, Kmin-1
             Ip(I) = Itemp(Ipp(I))
          End Do
       End If

       If (Kmin < Isize-1) Then
          CALL Qsort(X(Kmin+1:Isize), Ipp(Kmin+1:Isize))
          Itemp = Ip
          Do I = Kmin+1, Isize
             Ip(I) = Itemp(Ipp(I)+Kmin)
          End Do
       End If
    End If

    If (Present(Ipt)) Ipt = Ip

    Return
  End Subroutine Qsort_DP


End MODULE NonNumeric

