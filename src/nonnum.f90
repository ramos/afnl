! MODULE for non numerical routines (sorting and locate)
!
! "THE BEER-WARE LICENSE":
! Alberto Ramos wrote this file. As long as you retain this 
! notice you can do whatever you want with this stuff. If we meet some 
! day, and you think this stuff is worth it, you can buy me a beer in 
! return. <alberto.ramos@cern.ch>

MODULE NonNumeric

  USE ISO_FORTRAN_ENV, Only : error_unit, output_unit, iostat_end
  USE NumTypes
  USE Time

  IMPLICIT NONE

  interface realloc
     module procedure realloc_D, realloc_I
  end interface realloc

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

  Interface WriteBuffer
     Module Procedure WriteBuffer, WriteBuffer2
  End Interface WriteBuffer

  Interface ReadBuffer
     Module Procedure ReadBuffer, ReadBuffer2
  End Interface ReadBuffer

  Interface AddBuffer
     Module Procedure AddBuffer, AddBufferS, AddbufferR, &
          & AddbufferD, AddbufferRV, AddbufferDV
  End Interface AddBuffer

  Interface GetOpt
     Module Procedure GetOptCh, GetOptInt, GetOptDP, GetOptSP, GetOptTest
  End Interface GetOpt

  Interface Hash
     Module Procedure Hash_i32, Hash_i64, Hash_f32, Hash_f64, &
          & Hash_z32, Hash_z64
  End Interface Hash

  Character :: HDR=':'
  Integer, Parameter :: MAXLN=1000

  Integer :: Start_Hash = 314159265
  Integer :: table(0:255) = &
       & (/ -520310904, -23640796 ,-1755350438,88677750,&
       & -858319660,-2101220352 ,-1124865239,-2095760803,&
       & 44142808,1428115983,-1051390829, 1708352218,&
       & 1722642428,774853257,-983812907,580229275, & 
       & 458811929,1385753749,-889423576,-202374012,&
       & -512184309 ,-1833420199,879783905,1064866508,87813449,&
       & -1317557734,1275168835,-1012938719,132108151,&
       & -1397571578,-980931160,-367568929,-941835371,&
       & -1659849734,1428861935,-1437799479 ,-779465134,&
       & 1667768843,1953079905,-1145752466,-12411307, & 
       & -64593167,1886342641,1930672847,530346494,&
       & 1723067341,491614680,1435901097,96193509,&
       & -2091922483 ,-1489648432 ,-468587754,1986973157,&
       & -913563874,-517392711,-1549896364,91708779,&
       & -598424656,1417853772,-36573354,1224210515,& 
       & 1798300502,-1606047633,184606139,798926046,&
       & -1452817147 ,-1536717465,1384522629,-387404411,&
       & 1331304314,810407370 ,-1909117988,-1244393878,&
       & 1791454215,1928670958,180738758, 456692149,&
       & -873748860,885539777,1337483815,1005448863, &
       & -756459815,65899348,554233642,-317322574,&
       & -1790463615,340348919,-488050614,1256839619,&
       & 1944139941,2069847744, 1979897372,-1110137825,&
       & -577503759,1607741816,2091098967,-2127124658,&
       & -1679820319,-1175785395,1917773833,-1773654459 ,& 
       & 553064530,1685687797,-402586688,547875482,1535445404,& 
       & -1876593190,283394127,455278762,462872717,-1736499462,& 
       & -431196092,990322979,1943716543,457632225,-1162460171,& 
       & 7224510,-316660299,1653121126,-356898455,-1520267046,& 
       & 1294113126,-702268712,1076591054,-559072099,952824385,& 
       & 188814664,-626504961,1064196292,-2034554489,724966511,& 
       & -2013666511,717684207,220007275,1686086786,2138089290,& 
       & -293936682,-1400968910,-1397103180,236648461,803609507,& 
       & -1953190615,-557468531,-1542736622,1638685862,1837384809,& 
       & 126242359,-1998084886,1926153284,1425372941,-1117635688,& 
       & 1476236195,1424983543,-342458768,1460677252,494183745, & 
       & -2011045553,-955133617 ,-2064459243 ,-1304826384 ,267700582, & 
       & -605258983,552131509,-1260164126,-104399828,1941920290, & 
       & 1610044994,-668369221,441807305,-482645958,1988645233, & 
       & 2072623063,704294290,685324416,831514048,-1566531683, & 
       & 1499992087, -770471118, 1941194843, -1422644807,-1519837591, & 
       & 64336491,-656418956,-977203574,2071121624,665990048, & 
       & 1231935040,-896111611,1221186953,1249395953,-1614887967, & 
       & -937293250, -1811954021, -1035737928, -534402710, -1391628134,& 
       & 1145936422,919172998,1652242472,-1012298421,-1067982730, & 
       & 66134802,-499962716,853827516,-28119185,-404107906, & 
       & 278680114,-841625883,1516954392,-1255512391,&
       & -1852979883 ,-185941740,500081968,261370080,-302237855, &
       & -40718178,-861784371,-298112727,-713517878,-605484984,&
       & 510758174,1516754003,-1720521400,-726037970,-450268245,&
       & -798457399,-378130769,739370423,-1556413841,-1526980292,&
       & -699204581,-916568279,-2083403346,-1942818214,1912912161,&
       & -1770934318,363338731,262615298,-1005128728,-1987768399 ,&
       & 1803066496,2106827396,2135399440,-1140689457,98958985,&
       & 1567478145,71616901,-1948064131,1369392574,1754779077,&
       & -1023731795,123565065,-1089699658,1777593670,101116517,&
       & 1621514258 /)


  Private Locate_IN, Locate_SP, &
       & Locate_DP, Insrt_IN, Insrt_SP, Insrt_DP, Swap_IN, &
       & Swap_SP, Swap_DP, NewQsort_IN, NewQsort_SP, NewQsort_DP, &
       & Partition_SP, Partition_IN, Partition_DP, Partition, HDR, &
       & MAXLN, Start_Hash, table, AddBufferS, AddbufferR, &
       & AddbufferD, AddbufferRV, AddbufferDV, aabort, pperror

CONTAINS

! ****************************************
! *
  subroutine realloc_Z(d, n)
! *
! ****************************************
    complex (kind=DPC), intent (inout), allocatable :: d(:)
    integer, intent (in) :: n
    integer :: nold
    complex (kind=DPC), allocatable :: dcp(:)
    
    nold = size(d)
    if (allocated(d)) then
       allocate(dcp(nold))
       dcp = d
       
       deallocate(d)
       allocate(d(n))
       if (n > nold) then
          d(1:nold)  = dcp(1:nold)
          d(nold+1:) = (0.0_DP, 0.0_DP)
       else
          d(1:n)  = dcp(1:n)
       end if
    else
       allocate(d(n))
    end if

    return
  end subroutine realloc_Z

! ****************************************
! *
  subroutine realloc_D(d, n)
! *
! ****************************************
    real (kind=DP), intent (inout), allocatable :: d(:)
    integer, intent (in) :: n
    integer :: nold
    real (kind=DP), allocatable :: dcp(:)
    
    nold = size(d)
    if (allocated(d)) then
       allocate(dcp(nold))
       dcp = d
       
       deallocate(d)
       allocate(d(n))
       if (n > nold) then
          d(1:nold)  = dcp(1:nold)
          d(nold+1:) = 0.0_DP
       else
          d(1:n)  = dcp(1:n)
       end if
    else
       allocate(d(n))
    end if

    return
  end subroutine realloc_D

! ****************************************
! *
  subroutine realloc_I(d, n)
! *
! ****************************************
    integer, intent (inout), allocatable :: d(:)
    integer, intent (in) :: n
    integer :: nold
    integer, allocatable :: dcp(:)
    
    nold = size(d)
    if (allocated(d)) then
       allocate(dcp(nold))
       dcp = d
       
       deallocate(d)
       allocate(d(n))
       if (n > nold) then
          d(1:nold)  = dcp(1:nold)
          d(nold+1:) = 0
       else
          d(1:n)  = dcp(1:n)
       end if
    else
       allocate(d(n))
    end if

    return
  end subroutine realloc_I

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

       ! Initialize the stack
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

! ***************************************
! *
  Integer Function Pt2N(P, Npt)
! *
! ***************************************
! * Tranforms a Ndimensional point in an
! * integer index.
! ***************************************  

    Integer, Intent (in) :: P(:), Npt(:)
    Integer :: Nfac, I, Ndim

    Ndim = Size(Npt)

    Pt2N = P(1)-1
    Nfac = 1
    Do I = 2, Ndim
       Nfac = Nfac * Npt(I-1)
       Pt2N = Pt2N + Nfac*(P(I) - 1)
    End Do

    Return
  End Function Pt2N

! ***************************************
! *
  Function N2Pt(Npos, Npt) Result (P)
! *
! ***************************************
! * Tranforms a Ndimensional point in an
! * integer index.
! ***************************************  

    Integer, Intent (in) :: Npos, Npt(:)
    Integer :: Ncp, Ndim, I

    Integer :: P(Size(Npt))

    Ndim = Size(Npt)

    Ncp = Npos
    Do I = 1, Ndim
       P(I) = Mod(Ncp, Npt(I)) + 1
       Ncp = Int(Ncp/Npt(I))
    End Do

    Return
  End Function N2Pt

! ***************************************
! *
  Integer Function NumberOfLines(fname) Result (Nlines)
! *
! ***************************************

    Character (len=*), Intent (in) :: fname

    Open (Unit=77, File=Trim(fname), ACTION="READ")
    
    Nlines = 0
    Do While (.True.)
       Read(77,*,END=10)
       Nlines = Nlines + 1
    End Do
10  Continue

    Close(77)

    Return
  End Function NumberOfLines

! ***************************************
! *
  Integer Function NumberOfColumns(fname) Result (Ncol)
! *
! ***************************************

    Character (len=*), Intent (in) :: fname

    Integer :: J
    Character (len=100000) :: foo
    Real (kind=DP) :: A

    Open (Unit=77, File=Trim(fname), ACTION="READ")
    
    Read(77,'(1A)')foo
    Close(77)
    
    Ncol=1
    Do While (.True.)
       Read(foo,*,ERR=20, END=20)(A,J=1,Ncol)
       Ncol = Ncol + 1
    End Do
20  Continue
    Ncol = Ncol - 1
    
    Return
  End Function NumberOfColumns

! ***************************************
! *
  Integer Function Hash_i32(Iarr, Ist) Result (Hash)
! *
! ***************************************
! * Computes a HASH for an integer array
! * in a system independent way. 
! * With the help of the transfer intrinsic
! * this routine generates a Checksum for 
! * any data
! ***************************************
  
    Integer (kind=4), Intent (in) :: Iarr(:)
    Integer, Optional, Intent (in) :: Ist

    Integer :: I, J1, J2, J3, J4
    
    If (Present(Ist)) Then
       Hash = Ist
    Else
       Hash = Start_Hash
    End If

    J1 = 0
    J2 = 0
    J3 = 0
    J4 = 0  
    Do I = 1, Size(Iarr)
       CALL mvbits(Iarr(I), 24,  8, J1, 0)
       CALL mvbits(Iarr(I), 16,  8, J2, 0)
       CALL mvbits(Iarr(I), 8,  8, J3, 0)
       CALL mvbits(Iarr(I), 0,  8, J4, 0)
       
       Hash = Ieor(Hash,table(J1))
       Hash = Ishftc(Hash, 3)
       Hash = Ieor(Hash,table(J2))
       Hash = Ieor(Hash,table(J3))
       Hash = Ishftc(Hash, 5)
       Hash = Ieor(Hash,table(J4))
    End Do

    Return
  End Function Hash_i32

! ***************************************
! *
  Integer Function Hash_i64(Iarr, Ist) Result (Hash)
! *
! ***************************************
! * Computes a HASH for an integer array
! * in a system independent way. 
! * With the help of the transfer intrinsic
! * this routine generates a Checksum for 
! * any data
! ***************************************
  
    Integer (kind=8), Intent (in) :: Iarr(:)
    Integer, Optional, Intent (in) :: Ist

    Integer :: I, Itmp(2)
    
    If (Present(Ist)) Then
       Hash = Ist
    Else
       Hash = Start_Hash
    End If

    Do I = 1, Size(Iarr)
       Itmp = Transfer(Iarr(I), Itmp)
       Hash = Hash_i32(Itmp,Hash)
    End Do

    Return
  End Function Hash_i64

! ***************************************
! *
  Integer Function Hash_f32(arr, Ist) Result (Hash)
! *
! ***************************************
! * Computes a HASH for an integer array
! * in a system independent way. 
! * With the help of the transfer intrinsic
! * this routine generates a Checksum for 
! * any data
! ***************************************
  
    Real (kind=4), Intent (in) :: arr(:)
    Integer, Optional, Intent (in) :: Ist

    Integer :: I, Itmp(1)
    
    If (Present(Ist)) Then
       Hash = Ist
    Else
       Hash = Start_Hash
    End If

    Do I = 1, Size(arr)
       Itmp = Transfer(arr(I), Itmp)
       Hash = Hash_i32(Itmp,Hash)
    End Do

    Return
  End Function Hash_f32

! ***************************************
! *
  Integer Function Hash_f64(arr, Ist) Result (Hash)
! *
! ***************************************
! * Computes a HASH for an integer array
! * in a system independent way. 
! * With the help of the transfer intrinsic
! * this routine generates a Checksum for 
! * any data
! ***************************************
  
    Real (kind=8), Intent (in) :: arr(:)
    Integer, Optional, Intent (in) :: Ist

    Integer :: I, Itmp(2)
    
    If (Present(Ist)) Then
       Hash = Ist
    Else
       Hash = Start_Hash
    End If

    Do I = 1, Size(arr)
       Itmp = Transfer(arr(I), Itmp)
       Hash = Hash_i32(Itmp,Hash)
    End Do

    Return
  End Function Hash_f64

! ***************************************
! *
  Integer Function Hash_z64(arr, Ist) Result (Hash)
! *
! ***************************************
! * Computes a HASH for an integer array
! * in a system independent way. 
! * With the help of the transfer intrinsic
! * this routine generates a Checksum for 
! * any data
! ***************************************
  
    Complex (kind=8), Intent (in) :: arr(:)
    Integer, Optional, Intent (in) :: Ist

    Integer :: I, Itmp(4)
    
    If (Present(Ist)) Then
       Hash = Ist
    Else
       Hash = Start_Hash
    End If

    Do I = 1, Size(arr)
       Itmp = Transfer(arr(I), Itmp)
       Hash = Hash_i32(Itmp,Hash)
    End Do

    Return
  End Function Hash_z64

! ***************************************
! *
  Integer Function Hash_z32(arr, Ist) Result (Hash)
! *
! ***************************************
! * Computes a HASH for an integer array
! * in a system independent way. 
! * With the help of the transfer intrinsic
! * this routine generates a Checksum for 
! * any data
! ***************************************
  
    Complex (kind=4), Intent (in) :: arr(:)
    Integer, Optional, Intent (in) :: Ist

    Integer :: I, Itmp(2)
    
    If (Present(Ist)) Then
       Hash = Ist
    Else
       Hash = Start_Hash
    End If

    Do I = 1, Size(arr)
       Itmp = Transfer(arr(I), Itmp)
       Hash = Hash_i32(Itmp,Hash)
    End Do

    Return
  End Function Hash_z32

! ***************************************
! *
  Subroutine WriteBuffer(buf, fn, Comment) 
! *
! ***************************************

    Integer, Intent (in) :: buf(:)
    Character (len=*), Intent (in) :: fn
    Character (len=*), Intent (in), Optional :: Comment
    Integer :: ifn
    Logical :: is_used

    ifn = 69
    Do 
       Inquire(ifn, Opened=is_used)
       If (is_used) Then
          ifn = ifn + 1
       Else
          Exit
       End If
    End Do


    Open (File=Trim(fn), Unit=ifn, Form='FORMATTED', &
         & Access='STREAM', Action="WRITE")
    
    Write(ifn,'(1A)')'####'
    Write(ifn,'(1X,1A,1A)'     )"# Date and Time: ", &
         & asctime(gettime())    
    If (Present(Comment)) Then
       Write(ifn,'(1X,1A,1A)'  )"# Comment:       ", Trim(Comment)
    End If
    Write(ifn,'(1A)')'####'
    Write(ifn,'(2A)')HDR, 'DATA'
    Close(ifn)

    Open (File=Trim(fn), Unit=ifn, Form='UNFORMATTED', &
         & Access='STREAM', Position='APPEND')

    Write(ifn)Size(buf)
    Write(ifn)buf
    Write(ifn)Hash(buf)

    Close(ifn)

    Return
  End Subroutine WriteBuffer

! ***************************************
! *
  Subroutine WriteBuffer2(buf, fn, Comment) 
! *
! ***************************************

    Integer, Intent (in) :: buf(:)
    Character (len=*), Intent (in) :: fn
    Character (len=*), Intent (in) :: Comment(:)

    Integer :: I
    Integer :: ifn
    Logical :: is_used

    ifn = 69
    Do 
       Inquire(ifn, Opened=is_used)
       If (is_used) Then
          ifn = ifn + 1
       Else
          Exit
       End If
    End Do

    Open (File=Trim(fn), Unit=ifn, Form='FORMATTED', &
         & Access='STREAM', Action="WRITE")
    
    Write(ifn,'(1A)')'####'
    Write(ifn,'(1X,1A,1A)'     )"# Date and Time: ", &
         & asctime(gettime())    
    Write(ifn,'(1X,1A,1A)'  )"# Comment:       ", Trim(Comment(1))
    Do I = 2, Size(Comment)
       Write(ifn,'(1X,1A,1A)'  )"#                ", Trim(Comment(I))
    End Do
    Write(ifn,'(1A)')'####'
    Write(ifn,'(2A)')HDR, 'DATA'
    Close(ifn)

    Open (File=Trim(fn), Unit=ifn, Form='UNFORMATTED', &
         & Access='STREAM', Position='APPEND')

    Write(ifn)Size(buf)
    Write(ifn)buf
    Write(ifn)Hash(buf)

    Close(ifn)

    Return
  End Subroutine WriteBuffer2

! ***************************************
! *
  Subroutine ReadBuffer(buf, fn, chk) 
! *
! ***************************************

    Integer, Intent (out), allocatable :: buf(:)
    Character (len=*), Intent (in) :: fn
    Logical, Intent (out), Optional :: chk

    Logical :: is_ok
    Integer :: Ihs

    CALL ReadBuffer2(buf, fn, is_ok, Ihs)
    If (Present(chk)) chk = is_ok

    Return
  End Subroutine ReadBuffer

! ***************************************
! *
  Subroutine ReadBuffer2(buf, fn, chk, Ihsh) 
! *
! ***************************************

    Integer, Intent (out), allocatable :: buf(:)
    Character (len=*), Intent (in) :: fn
    Logical, Intent (out) :: chk
    Integer, Intent (out) :: Ihsh

    Character :: RD
    Integer :: Npos, Nsz, Is, CHKsum

    Open (File=Trim(fn), Unit=69, Form='FORMATTED', &
         & Access='STREAM', Action="READ")
    
    Read(69,*, IOSTAT=Is)
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')

    Do
       Read(69,'(1A1)', IOSTAT=Is)RD
       If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
       If (RD == HDR) Exit
    End Do
    Inquire(69,POS=Npos)
    Close(69)

    Open (File=Trim(fn), Unit=69, Form='UNFORMATTED', &
         & Access='STREAM', Action="READ")

    Read(69, POS=NPOS, IOSTAT=Is)Nsz
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')

    If (Allocated(buf)) Deallocate(buf)
    Allocate(buf(Nsz))
    Read(69, IOSTAT=Is)buf
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
    
    Read(69)CHKSum
    Close(69)

    chk = .False.
    Ihsh = Hash(buf)
    If (CHKsum == Ihsh) chk = .True.

    Return
  End Subroutine ReadBuffer2

! ***************************************
! *
  Function OpenBuffer(fn, nnsz) 
! *
! ***************************************

    Character (len=*), Intent (in) :: fn
    Integer, Intent (out), Optional :: nnsz

    Integer :: OpenBuffer

    Character :: RD
    Integer :: Npos, Nsz, Is
    Logical :: is_used

    OpenBuffer = 56
    Do 
       Inquire(OpenBuffer, Opened=is_used)
       If (is_used) Then
          OpenBuffer = OpenBuffer + 1
       Else
          Exit
       End If
    End Do

    Open (File=Trim(fn), Unit=OpenBuffer, Form='FORMATTED', &
         & Access='STREAM', Action="READ")
    
    Do
       Read(OpenBuffer,'(1A1)', IOSTAT=Is)RD
       If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
       If (RD == HDR) Exit
    End Do
    Inquire(OpenBuffer,POS=Npos)
    Close(OpenBuffer)

    Open (File=Trim(fn), Unit=OpenBuffer, Form='UNFORMATTED', &
         & Access='STREAM', Action="READ")

    Read(OpenBuffer, POS=NPOS, IOSTAT=Is)Nsz
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')

    If (Present(nnsz)) nnsz = Nsz

    Return
  End Function OpenBuffer

! ***************************************
! *
  Subroutine AddBuffer(buf, fn) 
! *
! ***************************************

    Integer, Intent (in) :: buf(:)
    Character (len=*), Intent (in) :: fn

    Character :: RD
    Integer :: Npos, Nsz, Is, CHKsum, N, I
    Integer :: ifn
    Logical :: is_used

    ifn = 69
    Do 
       Inquire(ifn, Opened=is_used)
       If (is_used) Then
          ifn = ifn + 1
       Else
          Exit
       End If
    End Do

    CHKsum = GetHash(fn)
    Open (File=Trim(fn), Unit=ifn, Form='FORMATTED', &
         & Access='STREAM', Action="READ")
    
    Do
       Read(ifn,'(1A1)', IOSTAT=Is)RD
       If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
       If (RD == HDR) Exit
    End Do
    Inquire(ifn,POS=Npos)
    Close(ifn)

    Open (File=Trim(fn), Unit=ifn, Form='UNFORMATTED', &
         & Access='STREAM', Action="READ")

    Read(ifn, POS=NPOS, IOSTAT=Is)Nsz
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
    Close(ifn)

    Open (File=Trim(fn), Unit=ifn, Form='UNFORMATTED', &
         & Access='STREAM', Action="READWRITE")
    Write(ifn, POS=NPOS, IOSTAT=Is)Nsz+Size(buf)


    Do I = 1, Nsz
       Read(ifn)N
    End Do
    Write(ifn)buf
    Write(ifn)Hash(buf, CHKsum)


    Close(ifn)

    Return
  End Subroutine AddBuffer

! ***************************************
! *
  Subroutine AddBufferD(D, fn) 
! *
! ***************************************

    Real (kind=DP), Intent (in) :: D
    Character (len=*), Intent (in) :: fn

    Integer :: Nib
    Integer, Allocatable :: ib(:)

    Nib = Size(Transfer(D, ib))
    Allocate(ib(Nib))
    ib = Transfer(D, ib)

    CALL AddBuffer(ib, fn)
    Deallocate(ib)

    Return
  End Subroutine AddBufferD

! ***************************************
! *
  Subroutine AddBufferDV(D, fn) 
! *
! ***************************************

    Real (kind=DP), Intent (in) :: D(:)
    Character (len=*), Intent (in) :: fn

    Integer, Allocatable :: ib(:)
    Integer :: Nib

    Nib = Size(Transfer(D, ib))
    Allocate(ib(Nib))
    ib = Transfer(D, ib)

    CALL AddBuffer(ib, fn)
    Deallocate(ib)

    Return
  End Subroutine AddBufferDV

! ***************************************
! *
  Subroutine AddBufferR(D, fn) 
! *
! ***************************************

    Real (kind=SP), Intent (in) :: D
    Character (len=*), Intent (in) :: fn

    Integer, Allocatable :: ib(:)
    Integer :: Nib

    Nib = Size(Transfer(D, ib))
    Allocate(ib(Nib))
    ib = Transfer(D, ib)

    CALL AddBuffer(ib, fn)
    Deallocate(ib)

    Return
  End Subroutine AddBufferR

! ***************************************
! *
  Subroutine AddBufferRV(D, fn) 
! *
! ***************************************

    Real (kind=SP), Intent (in) :: D(:)
    Character (len=*), Intent (in) :: fn

    Integer, Allocatable :: ib(:)
    Integer :: Nib

    Nib = Size(Transfer(D, ib))
    Allocate(ib(Nib))
    ib = Transfer(D, ib)

    CALL AddBuffer(ib, fn)
    Deallocate(ib)

    Return
  End Subroutine AddBufferRV

! ***************************************
! *
  Subroutine AddBufferS(buf, fn) 
! *
! ***************************************

    Integer, Intent (in) :: buf
    Character (len=*), Intent (in) :: fn

    Character :: RD
    Integer :: Npos, Nsz, Is, CHKsum, N, I
    Integer :: ifn
    Logical :: is_used

    ifn = 69
    Do 
       Inquire(ifn, Opened=is_used)
       If (is_used) Then
          ifn = ifn + 1
       Else
          Exit
       End If
    End Do

    CHKsum = GetHash(fn)
    Open (File=Trim(fn), Unit=ifn, Form='FORMATTED', &
         & Access='STREAM', Action="READ")
    
    Do
       Read(ifn,'(1A1)', IOSTAT=Is)RD
       If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
       If (RD == HDR) Exit
    End Do
    Inquire(ifn,POS=Npos)
    Close(ifn)

    Open (File=Trim(fn), Unit=ifn, Form='UNFORMATTED', &
         & Access='STREAM', Action="READ")

    Read(ifn, POS=NPOS, IOSTAT=Is)Nsz
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
    Close(ifn)

    Open (File=Trim(fn), Unit=ifn, Form='UNFORMATTED', &
         & Access='STREAM', Action="READWRITE")
    Write(ifn, POS=NPOS, IOSTAT=Is)Nsz+1


    Do I = 1, Nsz
       Read(ifn)N
    End Do
    Write(ifn)buf
    Write(ifn)Hash((/buf/), CHKsum)


    Close(ifn)

    Return
  End Subroutine AddBufferS

! ***************************************
! *
  Function GetHash(fn) 
! *
! ***************************************

    Character (len=*), Intent (in) :: fn
    Integer :: GetHash

    Character :: RD
    Integer :: Is, Nsz, I, Npos

    Open (File=Trim(fn), Unit=69, Form='FORMATTED', &
         & Access='STREAM', Action="READ")
    
    Do
       Read(69,'(1A1)', IOSTAT=Is)RD
       If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
       If (RD == HDR) Exit
    End Do
    Inquire(69,POS=Npos)
    Close(69)

    Open (File=Trim(fn), Unit=69, Form='UNFORMATTED', &
         & Access='STREAM', Action="READ")

    Read(69, POS=NPOS, IOSTAT=Is)Nsz
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
    Do I = 1, Nsz
       Read(69)Is
    End Do
    Read(69)GetHash
    Close(69)

    Return
  End Function GetHash

! ***************************************
! *
  Subroutine SetTable(tb) 
! *
! ***************************************
    
    Integer, Intent (in) :: tb(0:255)

    table(0:255) = tb(0:255)

    Return
  End Subroutine SetTable

! ***************************************
! *
  Subroutine SetStartHash(Ih) 
! *
! ***************************************
    
    Integer, Intent (in) :: Ih

    Start_Hash = Abs(Ih)

    Return
  End Subroutine SetStartHash

! ***************************************
! *
  Function GetDT(fn) 
! *
! ***************************************

    Character (len=*), Intent (in) :: fn
    Character (len=25) :: GetDT

    Character (len=17) :: foo
    Integer :: Is

    Open (File=Trim(fn), Unit=69, Form='FORMATTED', &
         & Access='STREAM', Action="READ")
    
    Read(69,*, IOSTAT=Is)
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')

    Read(69,'(1X,1A,1A24)', IOSTAT=Is)foo, GetDT
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
    Close(69)

    Return
  End Function GetDT

! ***************************************
! *
  Subroutine GetComment(fn, cm) 
! *
! ***************************************

    Character (len=*), Intent (in) :: fn
    Character (len=*), Intent (out) :: cm

    Character (len=17) :: foo
    Integer :: Is

    Open (File=Trim(fn), Unit=69, Form='FORMATTED', &
         & Access='STREAM', Action="READ")
    
    Read(69,*, IOSTAT=Is)
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
    Read(69,*, IOSTAT=Is)
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')

    Read(69,'(1X,1A,1A)', IOSTAT=Is)foo, cm
    If (Is /= 0) CALL Aabort('[ReadBuffer]', 'Error reading buffer')
    Close(69)

    Return
  End Subroutine GetComment

! ***************************************
! *
  Function GetOptCh(opt, val) Result (GetOpt)
! *
! ***************************************

    Character (len=*), Intent (in)  :: opt
    Character (len=*), Intent (out) :: val
    Logical :: GetOpt

    Integer :: I, Narg, Ist
    Character (len=MAXLN) :: arg

    Narg   = Command_argument_count()

    GetOpt = .False.
    If (Narg == 0) Return

    I = 0
    Do 
       I = I + 1
       Call Get_command_argument(I, arg)
       
       If ( (arg == opt).and.(len(Trim(arg))==len(Trim(opt)))) Then
          I = I+1
          Call Get_command_argument(I, val, STATUS=Ist)
          If (Ist /= 0) CALL Aabort('GetOptCh', &
               & 'Fail retriving the argument '//Trim(opt))
          GetOpt = .True.
          Return
       End If
          
       If (I == Narg) Return
    End Do
    

    Return
  End Function GetOptCh

! ***************************************
! *
  Function GetOptInt(opt, Ival) Result (GetOpt)
! *
! ***************************************

    Character (len=*), Intent (in)  :: opt
    Integer, Intent (out) :: Ival
    Logical :: GetOpt

    Integer :: I, Narg, Ist
    Character (len=MAXLN) :: arg
    Character (len=10000) :: foo

    Narg   = Command_argument_count()

    GetOpt = .False.
    If (Narg == 0) Return

    I = 0
    Do 
       I = I + 1
       Call Get_command_argument(I, arg)
       
       If ( (arg == opt).and.(len(Trim(arg))==len(Trim(opt)))) Then
          I = I+1
          Call Get_command_argument(I, foo, STATUS=Ist)
          If (Ist /= 0) CALL Aabort('GetOptCh', &
               & 'Fail retriving the argument '//Trim(opt))

          Read(foo,*,IOSTAT=Ist)Ival
          If (Ist /= 0) CALL Aabort('GetOptCh', &
               & 'Argument error '//Trim(opt))

          GetOpt = .True.
          Return
       End If
          
       If (I == Narg) Return
    End Do
    

    Return
  End Function GetOptInt

! ***************************************
! *
  Function GetOptSP(opt, val) Result (GetOpt)
! *
! ***************************************

    Character (len=*), Intent (in)  :: opt
    Real (kind=SP), Intent (out) :: val
    Logical :: GetOpt

    Integer :: I, Narg, Ist
    Character (len=MAXLN) :: arg
    Character (len=10000) :: foo

    Narg   = Command_argument_count()

    GetOpt = .False.
    If (Narg == 0) Return

    I = 0
    Do 
       I = I + 1
       Call Get_command_argument(I, arg)
       
       If ( (arg == opt).and.(len(Trim(arg))==len(Trim(opt)))) Then
          I = I+1
          Call Get_command_argument(I, foo, STATUS=Ist)
          If (Ist /= 0) CALL Aabort('GetOptCh', &
               & 'Fail retriving the argument '//Trim(opt))

          Read(foo,*,IOSTAT=Ist)val
          If (Ist /= 0) CALL Aabort('GetOptCh', &
               & 'Argument error '//Trim(opt))

          GetOpt = .True.
          Return
       End If
          
       If (I == Narg) Return
    End Do
    

    Return
  End Function GetOptSP

! ***************************************
! *
  Function GetOptDP(opt, val) Result (GetOpt)
! *
! ***************************************

    Character (len=*), Intent (in)  :: opt
    Real (kind=DP), Intent (out) :: val
    Logical :: GetOpt

    Integer :: I, Narg, Ist
    Character (len=MAXLN) :: arg
    Character (len=10000) :: foo

    Narg   = Command_argument_count()

    GetOpt = .False.
    If (Narg == 0) Return

    I = 0
    Do 
       I = I + 1
       Call Get_command_argument(I, arg)
       
       If ( (arg == opt).and.(len(Trim(arg))==len(Trim(opt)))) Then
          I = I+1
          Call Get_command_argument(I, foo, STATUS=Ist)
          If (Ist /= 0) CALL Aabort('GetOptCh', &
               & 'Fail retriving the argument '//Trim(opt))

          Read(foo,*,IOSTAT=Ist)val
          If (Ist /= 0) CALL Aabort('GetOptCh', &
               & 'Argument error '//Trim(opt))

          GetOpt = .True.
          Return
       End If
          
       If (I == Narg) Return
    End Do
    

    Return
  End Function GetOptDP

! ***************************************
! *
  Function GetOptTest(opt) Result (GetOpt)
! *
! ***************************************

    Character (len=*), Intent (in)  :: opt
    Logical :: GetOpt

    Integer :: I, Narg
    Character (len=MAXLN) :: arg

    Narg   = Command_argument_count()

    GetOpt = .False.
    If (Narg == 0) Return

    I = 0
    Do 
       I = I + 1
       Call Get_command_argument(I, arg)
       
       If (arg == opt) Then
          GetOpt = .True.
          Return
       End If
          
       If (I == Narg) Return
    End Do
    

    Return
  End Function GetOptTest

! ***************************************************
! *
  Subroutine aAbort(Routine, Msg)
! *
! ***************************************************
  
    Character (len=*), Intent (in), Optional :: Routine
    Character (len=*), Intent (in) :: Msg
    
    If (Present(Routine)) Then
       Write(error_unit, *)'  Abort: IN ', Trim(Routine),': ',Msg
    Else
       Write(error_unit, *)'  Abort: ', Msg
    End If

    Stop
    
    Return
  End Subroutine AAbort

! ***************************************************
! *
  Subroutine pPerror(Routine, Msg)
! *
! ***************************************************
  
    Character (len=*), Intent (in), Optional :: Routine
    Character (len=*), Intent (in) :: Msg
    
    If (Present(Routine)) Then
       Write(error_unit, *)'  Error: IN ', Trim(Routine),': ',Msg
    Else
       Write(error_unit, *)'  Error: ', Msg
    End If

    Return
  End Subroutine PPerror

! ***************************************************
! *
  Subroutine addtofile(orig, dest)
! *
! ***************************************************
    character (len=*), intent (in) :: orig, dest

    logical :: ok
    integer :: iszo, iszd, ifnd, ifno, ios
    character (len=1) :: c
    
    inquire(file=trim(orig), exist=ok)
    if (.not.ok) call aabort('addtofile', 'file '//trim(orig)//' does not exist')
    inquire(file=trim(dest), exist=ok)
    if (.not.ok) call aabort('addtofile', 'file '//trim(dest)//' does not exist')

    inquire(file=trim(orig), size=iszo)
    inquire(file=trim(dest), size=iszd)
    if (iszd>iszo) call aabort('addtofile', 'destination file larger than origin file')
    
    open(file=Trim(orig), newunit=ifno, ACTION="READ", &
              & Form='UNFORMATTED', Access='STREAM')
    open(file=Trim(dest), newunit=ifnd, ACTION="READWRITE", &
              & Form='UNFORMATTED', Access='STREAM')

    read(ifno, pos=iszd)
    read(ifnd, pos=iszd)
    do
       read(ifno, iostat=ios)c
       if (ios==iostat_end) exit
       write(ifnd)c
    end do
    
    close(ifno)
    close(ifnd)

    return
  end Subroutine addtofile
  
End MODULE NonNumeric

