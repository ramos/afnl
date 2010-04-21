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

  IMPLICIT NONE

  Interface Swap
     Module Procedure Swap_IN, Swap_SP, Swap_DP
  End Interface

  Interface Insrt
     Module Procedure Insrt_IN, Insrt_SP, Insrt_DP
  End Interface

  Interface Qsort
     Module Procedure NewQsort_IN, NewQsort_SP, NewQsort_DP
  End Interface

  Interface Locate
     Module Procedure Locate_IN, Locate_SP, Locate_DP
  End Interface

  Interface Partition
     Module Procedure Partition_SP, Partition_IN, Partition_DP
  End Interface


  Private Locate_IN, Locate_SP, &
       & Locate_DP, Insrt_IN, Insrt_SP, Insrt_DP, Swap_IN, &
       & Swap_SP, Swap_DP, NewQsort_IN, NewQsort_SP, NewQsort_DP, &
       & Partition_SP, Partition_IN, Partition_DP, Partition

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

    Integer :: Nmin, Nmax, Nnew, Isize


    ISize = Size(X)

    If (Present(Iin)) Then
       Locate_IN = IIn
    End If

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

    Integer :: Nmin, Nmax, Nnew, Isize


    ISize = Size(X)

    If (Present(Iin)) Then
       Locate_SP = IIn
    End If

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

    Integer :: Nmin, Nmax, Nnew, Isize


    ISize = Size(X)

    If (Present(Iin)) Then
       Locate_DP = IIn
    End If

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
  Subroutine Insrt_IN(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order.
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. Integer 
! * version.
! ***********************************

    Integer, Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Integer :: Rtmp
    Integer :: I, J


    If (Present(Ipt)) Then
       Forall (I=1:Size(X)) Ipt(I) = I
       
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
                CALL Swap(Ipt, J, J+1)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    Else
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    End If

    Return
  End Subroutine Insrt_IN


! ***********************************
! *
  Subroutine Insrt_SP(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order.
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. Integer 
! * version.
! ***********************************

    Real (kind=SP), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Real (kind=SP) :: Rtmp
    Integer :: I, J


    If (Present(Ipt)) Then
       Forall (I=1:Size(X)) Ipt(I) = I
       
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
                CALL Swap(Ipt, J, J+1)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    Else
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    End If

    Return
  End Subroutine Insrt_SP

! ***********************************
! *
  Subroutine Insrt_DP(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order.
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. Integer 
! * version.
! ***********************************

    Real (kind=DP), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Real (kind=DP) :: Rtmp
    Integer :: I, J


    If (Present(Ipt)) Then
       Forall (I=1:Size(X)) Ipt(I) = I
       
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
                CALL Swap(Ipt, J, J+1)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    Else
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    End If

    Return
  End Subroutine Insrt_DP

! ***********************************
! *
  Subroutine Swap_SP(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************

    Real (kind=SP), Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J

    Real (kind=SP) :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap_SP

! ***********************************
! *
  Subroutine Swap_DP(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************

    Real (kind=DP), Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J

    Real (kind=DP) :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap_DP

! ***********************************
! *
  Subroutine Swap_IN(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************

    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J

    Integer :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap_IN

! ***********************************
! *
  Subroutine NewQsort_SP(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. SP version
! ***********************************

    Type Limits
       Integer :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer, Parameter :: Isw = 10

    Real (kind=SP), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Integer :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Type (Limits), Allocatable :: Stack(:)
    
    
    Allocate(Stack(2*Size(X)))

    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)
!          Write(*,*)Ispos, ISmax

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrt(X(Ileft:Iright))
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    End If

    Deallocate(Stack)

    Return
    
  CONTAINS

    ! ***********************************
    Integer Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright)
    ! * for Qsort.
    ! ***********************************
      
      Real (kind=SP), Intent (in) :: XX(:)
      Integer, Intent (in) :: IIleft, IIright
      
      Real (kind=SP) :: XXcp(3)
      Integer :: IIpt(3), IImd
      
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IIright)
      XXcp(3) = XX(IImd)
      
      CALL Insrt(XXcp, IIpt)
      
      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return
    End Function ChoosePiv

    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
    ! ***********************************

      Real (kind=SP), Intent (inout) :: XX(:)
      Integer, Intent (inout) :: IIpt(:)
      Integer, Intent (in) :: IIl, IIr
      
      Real (kind=SP) :: RRtmp
      Integer :: II, JJ
      

      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do
      
      Return
    End Subroutine InsrtLC


  End Subroutine NewQsort_SP

! ***********************************
! *
  Subroutine NewQsort_DP(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. DP version
! ***********************************

    Type Limits
       Integer :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer, Parameter :: Isw = 10

    Real (kind=DP), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Integer :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Type (Limits), Allocatable :: Stack(:)
    
    
    Allocate(Stack(2*Size(X)))



    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)
!          Write(*,*)Ispos, ISmax

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrt(X(Ileft:Iright))
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    End If

    Deallocate(Stack)

    Return
    
  CONTAINS

    ! ***********************************
    Integer Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright)
    ! * for Qsort.
    ! ***********************************
      
      Real (kind=DP), Intent (in) :: XX(:)
      Integer, Intent (in) :: IIleft, IIright
      
      Real (kind=DP) :: XXcp(3)
      Integer :: IIpt(3), IImd
      
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IIright)
      XXcp(3) = XX(IImd)
      
      CALL Insrt(XXcp, IIpt)
      
      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return
    End Function ChoosePiv

    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
    ! ***********************************

      Real (kind=DP), Intent (inout) :: XX(:)
      Integer, Intent (inout) :: IIpt(:)
      Integer, Intent (in) :: IIl, IIr
      
      Real (kind=DP) :: RRtmp
      Integer :: II, JJ
      

      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do
      
      Return
    End Subroutine InsrtLC


  End Subroutine NewQsort_DP

! ***********************************
! *
  Subroutine NewQsort_IN(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. DP version
! ***********************************

    Type Limits
       Integer :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer, Parameter :: Isw = 10

    Integer, Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
    
    Integer :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Type (Limits), Allocatable :: Stack(:)
    
    
    Allocate(Stack(2*Size(X)))



    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)
!          Write(*,*)Ispos, ISmax

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrt(X(Ileft:Iright))
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    End If

    Deallocate(Stack)

    Return
    
  CONTAINS

    ! ***********************************
    Integer Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright)
    ! * for Qsort.
    ! ***********************************
      
      Integer, Intent (in) :: XX(:)
      Integer, Intent (in) :: IIleft, IIright
      
      Real (kind=DP) :: XXcp(3)
      Integer :: IIpt(3), IImd
      
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IIright)
      XXcp(3) = XX(IImd)
      
      CALL Insrt(XXcp, IIpt)
      
      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return
    End Function ChoosePiv

    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
    ! ***********************************

      Integer, Intent (inout) :: XX(:)
      Integer, Intent (inout) :: IIpt(:)
      Integer, Intent (in) :: IIl, IIr
      
      Integer :: RRtmp
      Integer :: II, JJ
      

      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do
      
      Return
    End Subroutine InsrtLC


  End Subroutine NewQsort_IN

! ***********************************
! *
  Integer Function Partition_SP(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by Qsort.
! ***********************************

    Real (kind=SP), Intent (inout) :: X(:)
    Integer, Intent (in) :: Ileft, Iright, Ipv
    Integer, Intent (inout), Optional :: Ipt(:)
    
    Real (kind=SP) :: Rpv
    Integer :: I

    Rpv = X(Ipv)
    CALL Swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             CALL Swap(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL Swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap(Ipt, Ipvfn, Iright)

    Return
  End Function Partition_SP

! ***********************************
! *
  Integer Function Partition_DP(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by Qsort.
! ***********************************

    Real (kind=DP), Intent (inout) :: X(:)
    Integer, Intent (in) :: Ileft, Iright, Ipv
    Integer, Intent (inout), Optional :: Ipt(:)
    
    Real (kind=DP) :: Rpv
    Integer :: I

    Rpv = X(Ipv)
    CALL Swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             CALL Swap(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL Swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap(Ipt, Ipvfn, Iright)

    Return
  End Function Partition_DP

! ***********************************
! *
  Integer Function Partition_IN(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by Qsort.
! ***********************************

    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: Ileft, Iright, Ipv
    Integer, Intent (inout), Optional :: Ipt(:)
    
    Integer :: Rpv
    Integer :: I

    Rpv = X(Ipv)
    CALL Swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             CALL Swap(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL Swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap(Ipt, Ipvfn, Iright)

    Return
  End Function Partition_IN

End MODULE NonNumeric

