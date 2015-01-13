!
! A MODULE for RANLUX Random Number Generator
!
! This module implements the RANLUX Random Number Generator in 
! standard fortran 2008. The relevant references are:
!
! (*) A portable high-quality random number generator for 
!     lattice field theory simulations.
!     M. Lüscher, Computer Physics Communications 79 (1994) 100 
!
! (*) RANLUX number generator.
!     M. Lüscher, http://luscher.web.cern.ch/luscher/ranlux/
!
! The code has been compared with the reference implementation,
! but as always the correctness of the code can not be guaranteed. 
! Use with caution.
!
! "THE BEER-WARE LICENSE":
! Alberto Ramos wrote this file. As long as you retain this 
! notice you can do whatever you want with this stuff. If we meet some 
! day, and you think this stuff is worth it, you can buy me a beer in 
! return. <alberto.ramos@desy.de>
!

! ***************************************************
! *
MODULE RanLux48
! *
! ***************************************************

  USE ISO_FORTRAN_ENV, Only : error_unit, output_unit
  IMPLICIT NONE

  private

  type :: state
     integer (kind=8) :: V(96)
     integer (kind=8) :: c(8)
     integer          :: cnt, pr, ir, is
     integer          :: ilv, SEED_TYPE = -1
     integer (kind=8) :: sd(4) = -1, otherseed = -1
     integer (kind=8) :: cyc = 0_8
  end type state

  integer (kind=8), parameter :: BASE = 281474976710656_8
  real (kind=8), parameter :: DINV = 3.5527136788005135511301048E-15_8

  type (state), save :: rnd
  integer, parameter :: next(1:12) = (/12,1,2,3,4,5,6,7,8,9,10,11/)

  Integer, Parameter, Private :: SEED_LCG=0, SEED_MINE=1, &
       & SEED_V3=2

  interface rnlx48
     module procedure :: rnlx48_f64, rnlx48_f64_v, rnlx48_f32, &
          rnlx48_f32_v, rnlx48_z32, rnlx48_z32_v, rnlx48_z64, &
          rnlx48_z64_v
  end interface rnlx48

  public :: rnlx48_init, rnlx48, rnlx48_seed_lcg, rnlx48_seed_v3, &
       rnlx48_seed_mine, rnlx48_print_info, rnlx48_size, rnlx48_get, rnlx48_reset, &
       rnlx48_f32, rnlx48_f32_v, rnlx48_f64, rnlx48_f64_v, &
       rnlx48_z32, rnlx48_z32_v, rnlx48_z64, rnlx48_z64_v

  !$OMP THREADPRIVATE(rnd)

CONTAINS

! ***************************************************
! *
  Function rnlx48_size()
! *
! ***************************************************
    integer :: rnlx48_size

    rnlx48_size = 116

    return
  end Function rnlx48_size

! ***************************************************
! *
  Subroutine rnlx48_get(state)
! *
! ***************************************************
    
    Integer (kind=8), intent (out) :: state(:)

    if (size(state)<rnlx48_size()) call rnlx48_error('rnlx48_get', &
         'Array too small to store the state of RANLUX')
    
    state(1:96)    = int(rnd%V(1:96)  ,kind=8)
    state(97:104)  = int(rnd%c(1:8)   ,kind=8)
    state(105)     = int(rnd%cnt      ,kind=8)
    state(106)     = int(rnd%pr       ,kind=8)
    state(107)     = int(rnd%ir       ,kind=8)
    state(108)     = int(rnd%is       ,kind=8)

    state(109)     = int(rnd%ilv      ,kind=8)
    state(110)     = int(rnd%otherseed,kind=8)
    state(111:114) = int(rnd%sd(1:4)  ,kind=8)
    state(115)     = int(rnd%SEED_TYPE,kind=8)
    state(116)     = int(rnd%cyc      ,kind=8)

    Return
  End Subroutine rnlx48_get

! ***************************************************
! *
  Subroutine rnlx48_reset(state)
! *
! ***************************************************
    
    Integer (kind=8), intent (in) :: state(:)

    if (size(state)<rnlx48_size()) call rnlx48_error('rnlx48_reset', &
         'Array too small to store the state of RANLUX')
    
    rnd%V(1:96)    = state(1:96) 
    rnd%c(1:8)     = state(97:104)
    rnd%cnt        = int(state(105)    ,kind=4)
    rnd%pr         = int(state(106)    ,kind=4)
    rnd%ir         = int(state(107)    ,kind=4)
    rnd%is         = int(state(108)    ,kind=4)

    rnd%ilv        = int(state(109)    ,kind=4)
    rnd%otherseed  = state(110)
    rnd%sd(1:4)    = state(111:114)
    rnd%SEED_TYPE  = int(state(115)    ,kind=4)
    rnd%cyc        = state(116)

    Return
  End Subroutine rnlx48_reset


! ***************************************************
! *
  Subroutine rnlx48_seed_lcg(nseed)
! *
! ***************************************************
    
    Integer (kind=8), intent (in) :: nseed
    Integer :: I

    rnd%V(1) = mod(nseed,BASE)
    Do I = 2, 96
       rnd%V(I) = mod(rnd%V(I-1)*13,BASE)
    End Do
    call fill_rnd()
    rnd%cnt = 97
    rnd%SEED_TYPE = SEED_LCG
    rnd%otherseed = nseed
    rnd%cyc = 0_8
    
    Return
  End Subroutine rnlx48_seed_lcg


! ***************************************************
! *
  Subroutine rnlx48_seed_v3(nseed)
! *
! ***************************************************
    
    Integer, intent (in) :: nseed
    
    integer (kind=8) :: ix
    integer :: ibt, jbt, bits(0:30), ns, i, j, k, ib

    if ((nseed>2147483647).or.(nseed<1)) call rnlx48_error(&
         'rnlx48_seed_mine',&
         'First component of seed too large (1<seed(1)<2147483647)')

    ns = nseed
    do i = 0, 30
       bits(i) = mod(ns,2)
       ns = ns/2
    end do
    
    ibt = 0
    jbt = 18

    do k = 0, 7
       do j = 0, 11
          ix = 0
          do i = 0, 47
             ib = bits(ibt)
             ix = 2_8*ix+int(ib,kind=8)
             
             bits(ibt) = mod(bits(ibt)+bits(jbt),2)
             ibt = mod(ibt+1,31)
             jbt = mod(jbt+1,31)
          end do
          if (mod(j,8)==k) ix = BASE-1_8-ix
          rnd%V(12-j+k*12) = ix
       end do
    end do
    rnd%cnt = 97
    rnd%SEED_TYPE = SEED_V3
    rnd%otherseed = nseed
    rnd%cyc = 0_8
    
    Return
  End Subroutine rnlx48_seed_v3


! ***************************************************
! *
  Subroutine rnlx48_seed_mine(nseed)
! *
! ***************************************************
    Integer, intent (in) :: nseed(:)
    
    integer (kind=8) :: ix
    integer :: ibt, jbt, bits(0:30), ns(4), ib, i, j, k, &
         isz, jj, ic

    isz   = size(nseed)
    ns(1) = nseed(1)

    call checkseed()
    write(*,*)isz, ns
    rnd%SEED_TYPE = SEED_MINE
    rnd%sd(:) = ns(:) 

    rnd%V = 0
    do ic=1,4
       do i = 0, 30
          bits(i) = mod(ns(ic),2)
          ns(ic) = ns(ic)/2
       end do
    
       ibt = 0
       jbt = 18
       
       do k = 0, 7
          do jj = 0, 2
             j = jj+3*(ic-1)
             ix = 0
             do i = 0, 47
                ib = bits(ibt)
                ix = 2*ix+ib
                
                bits(ibt) = mod(bits(ibt)+bits(jbt),2)
                ibt = mod(ibt+1,31)
                jbt = mod(jbt+1,31)
             end do
             if (mod(j,8)==k) ix = BASE-1_8-ix
             rnd%V(12-j+k*12) = ix
          end do
       end do
       
    end do

    rnd%cnt = 97

    Return
  contains
    
    subroutine checkseed()

      if (isz>4) call rnlx48_error('rnlx48_seed_mine','Too many numbers in seed')
      if (nseed(1)>2147483647) call rnlx48_error('rnlx48_seed_mine',&
           'First component of seed too large (1<seed(1)<2147483647)')
      if ((isz>1)) then 
         if ((nseed(2)>2147483646)) call rnlx48_error('rnlx48_seed_mine',&
              'Second component of seed out of range (1<seed(2)<2147483647)')
         ns(2) = nseed(2)+1
      else
         ns(2) = 1
      end if
      if ((isz>2)) then 
         if ((nseed(3)>2147483646)) call rnlx48_error('rnlx48_seed_mine',&
              'Third component of seed out of range (1<seed(3)<2147483647)')
         ns(3) = nseed(3)+1
      else
         ns(3) = 1
      end if
      if ((isz>3)) then 
         if ((nseed(4)>2147483646)) call rnlx48_error('rnlx48_seed_mine',&
              'Fourth component of seed out of range (1<seed(4)<2147483647)')
         ns(4) = nseed(4)+1
      else
         ns(4) = 1
      end if

      return
    end subroutine checkseed

  End Subroutine rnlx48_seed_mine

! ***************************************************
! *
  Subroutine rnlx48_init(n)
! *
! ***************************************************
    integer, intent (in), optional :: n
    integer :: np

    rnd%is  = 5
    rnd%ir  = 12
    rnd%c   = 0
    rnd%cnt = 1

    if (present(n)) then
       if (n>12) then
          rnd%pr = n
          rnd%ilv = 3
          return
       else
          np = n
       end if
    else
       np = 1
    end if

    select case (np)
    case (-1)
       rnd%pr = 12
       rnd%ilv = -1
    case (0)
       rnd%pr = 109
       rnd%ilv = 0
    case (1)
       rnd%pr = 202
       rnd%ilv = 1
    case (2)
       rnd%pr = 397
       rnd%ilv = 2
    end select

    return
  end Subroutine rnlx48_init

! ***************************************************
! *
  Subroutine fill_rnd()
! *
! ***************************************************
    integer :: j, ic
    integer (kind=8) :: delta(8)

    do ic = 1, rnd%pr
       do concurrent (j=0:7) 
          delta(j+1) = rnd%V(rnd%is+12*j) - rnd%V(rnd%ir+12*j) - rnd%c(j+1)
       end do

       where (delta < 0)
          rnd%c = 1
          delta = delta + BASE
       elsewhere
          rnd%c = 0
       end where

       do concurrent (j=0:7) 
          rnd%V(rnd%ir+12*j) = delta(j+1)
       end do

       rnd%ir = next(rnd%ir)
       rnd%is = next(rnd%is)
    end do
    rnd%cnt = 1
    rnd%cyc = rnd%cyc+1_8

    return
  end Subroutine fill_rnd

! ***************************************************
! *
  Subroutine rnlx48_int(n)
! *
! ***************************************************

    Integer (kind=4), Intent (out) :: n
    Integer (kind=4) :: ipos
    
    ipos = rnd%cnt
    If (ipos > 96) Then
       ipos = 1
       CALL fill_rnd()
    End If
    rnd%cnt = rnd%cnt + 1
    n = int(rnd%V(ipos),kind=4)

    Return
  End Subroutine rnlx48_int

! ***************************************************
! *
  Subroutine rnlx48_int_v(d)
! *
! ***************************************************

    Real (kind=4), Intent (out) :: d(:)
    Integer (kind=4) :: i, nleft, nst, ncyc
    
    nleft = Size(d)
    nst = 96-rnd%cnt+1
    If (nst > nleft) Then
       d(1:nleft) = int(rnd%V(rnd%cnt:rnd%cnt+nleft-1),kind=4)
       rnd%cnt = rnd%cnt+nleft
       Return
    Else if (nst>0) then
       d(1:nst) = int(rnd%V(rnd%cnt:rnd%cnt+nst-1),kind=4)
       nleft = nleft - nst
    else
       nst = 0
    End If
    CALL fill_rnd() 
    
    ncyc = Int(nleft/96)
    Do I = 1, ncyc
       d(nst+1+(I-1)*96:nst+I*96) = int(rnd%V(1:96),kind=4)
       CALL fill_rnd()
    End Do
    nleft = Modulo(nleft,96)
    d(nst+1+ncyc*96:Size(d)) = int(rnd%V(1:nleft),kind=4)
    rnd%cnt = nleft+1

    Return
  End Subroutine rnlx48_int_v

! ***************************************************
! *
  Subroutine rnlx48_f64(d)
! *
! ***************************************************

    Real (kind=8), Intent (out) :: d
    Integer (kind=4) :: ipos
    integer (kind=8) :: i1
    
    ipos = rnd%cnt
    If (ipos > 96) Then
       CALL fill_rnd()
       i1 = rnd%V(1)
    else
       i1 = rnd%V(ipos)
    End If
    rnd%cnt = rnd%cnt + 1
    d = Real(i1,kind=8)*Real(DINV,kind=8)

    Return
  End Subroutine rnlx48_f64

! ***************************************************
! *
  Subroutine rnlx48_f32(d)
! *
! ***************************************************

    Real (kind=4), Intent (out) :: d
    Integer (kind=4) :: ipos
    integer (kind=8) :: i1
    
    ipos = rnd%cnt
    If (ipos > 96) Then
       CALL fill_rnd()
       i1 = rnd%V(1)
    else
       i1 = rnd%V(ipos)
    End If
    rnd%cnt = rnd%cnt + 1
    d = Real(i1,kind=4)*Real(DINV,kind=4)

    Return
  End Subroutine rnlx48_f32

! ***************************************************
! *
  Subroutine rnlx48_z32(d)
! *
! ***************************************************

    Complex (kind=4), Intent (out) :: d
    Real (kind=4) :: x, y

    CALL rnlx48(x)
    CALL rnlx48(y)
    d = Cmplx(x,y)
    
    Return
  End Subroutine rnlx48_z32

! ***************************************************
! *
  Subroutine rnlx48_z64(d)
! *
! ***************************************************

    Complex (kind=8), Intent (out) :: d
    Real (kind=8) :: x, y

    CALL rnlx48(x)
    CALL rnlx48(y)
    d = Cmplx(x,y,kind=8)
    
    Return
  End Subroutine rnlx48_z64

! ***************************************************
! *
  Subroutine rnlx48_z32_v(d)
! *
! ***************************************************

    Complex (kind=4), Intent (out) :: d(:)
    Real (kind=4) :: x(Size(d))
    
    CALL rnlx48(x)
    d = Cmplx(x,0.0_4)
    CALL rnlx48(x)
    d = d + Cmplx(0.0_4,x)

    Return
  End Subroutine rnlx48_z32_v

! ***************************************************
! *
  Subroutine rnlx48_z64_v(d)
! *
! ***************************************************

    Complex (kind=8), Intent (out) :: d(:)
    Real (kind=8) :: x(Size(d))
    
    CALL rnlx48(x)
    d = Cmplx(x,0.0_8,kind=8)
    CALL rnlx48(x)
    d = d + Cmplx(0.0_8,x,kind=8)

    Return
  End Subroutine rnlx48_z64_v

! ***************************************************
! *
  Subroutine rnlx48_f64_v(d)
! *
! ***************************************************

    Real (kind=8), Intent (out) :: d(:)
    Integer (kind=4) :: i, nleft, nst, ncyc, iex

    nleft = Size(d)
    nst = 96-rnd%cnt+1
    If (nst > nleft) Then
       d(1:nleft) = Real(rnd%V(rnd%cnt:rnd%cnt+nleft-1),kind=8)*DINV
       rnd%cnt = rnd%cnt+nleft
       Return
    Else if (nst>0) then
       d(1:nst) = Real(rnd%V(rnd%cnt:rnd%cnt+nst-1),kind=8)*DINV
       nleft = nleft - nst
    else
       nst = 0
    End If
    CALL fill_rnd() 
    
    ncyc = Int(nleft/96)
    Do I = 1, ncyc
       d(nst+1+(I-1)*96:nst+I*96) = &
            & Real(rnd%V(1:96),kind=4)*DINV
       CALL fill_rnd()
    End Do
    nleft = Modulo(nleft,96)
    d(nst+1+ncyc*96:Size(d)) = &
         & Real(rnd%V(1:nleft),kind=4)*DINV
    rnd%cnt = nleft+1

    Return
  End Subroutine rnlx48_f64_v

! ***************************************************
! *
  Subroutine rnlx48_f32_v(d)
! *
! ***************************************************

    Real (kind=4), Intent (out) :: d(:)
    Integer (kind=4) :: i, nleft, nst, ncyc
    
    nleft = Size(d)
    nst = 96-rnd%cnt+1
    If (nst > nleft) Then
       d(1:nleft) = Real(rnd%V(rnd%cnt:rnd%cnt+nleft-1),kind=4)&
            & * Real(DINV,kind=4)
       rnd%cnt = rnd%cnt+nleft
       Return
    Else if (nst>0) then
       d(1:nst) = Real(rnd%V(rnd%cnt:rnd%cnt+nst-1),kind=4) &
            & * Real(DINV,kind=4)
       nleft = nleft - nst
    else
       nst = 0
    End If
    CALL fill_rnd() 
    
    ncyc = Int(nleft/96)
    Do I = 1, ncyc
       d(nst+1+(I-1)*96:nst+I*96) = &
            & Real(rnd%V(1:96),kind=4)*Real(DINV,kind=4)
       CALL fill_rnd()
    End Do
    nleft = Modulo(nleft,96)
    d(nst+1+ncyc*96:Size(d)) = &
         & Real(rnd%V(1:nleft),kind=4)*Real(DINV,kind=4)
    rnd%cnt = nleft+1

    Return
  End Subroutine rnlx48_f32_v

! ***************************************************
! *
  Subroutine rnlx48_print_info(ifn)
! *
! ***************************************************

    Integer, Intent (in), Optional :: ifn
    Integer :: iout
    integer (kind=8) :: nnumbers
    
    iout = output_unit
    if (present(ifn)) iout = ifn

    Write(iout,'(1A,1I4)')'RANLUX Working with luxury level: ', rnd%ilv
    Write(iout,'(1A,1I6)')'RANLUX luxury parameter:  ', 2*rnd%pr
    Write(iout,'(1A)', ADVANCE="NO")'RANLUX seeded with method: '

    Select Case (rnd%SEED_TYPE)
    Case (SEED_LCG)
       Write(iout,'(1A)', ADVANCE="NO")'seed LCG'
       Write(iout,'(1I11,1A)')rnd%otherseed, ' )'
    Case (SEED_V3)
       Write(iout,'(1A)', ADVANCE="NO")'seed ranlux v3.0 ('
       Write(iout,'(1I11,1A)')rnd%otherseed, ' )'
    Case (SEED_MINE)
       Write(iout,'(1A)', ADVANCE="NO")'seed 4-component ('
       Write(iout,'(4I11,1A)')rnd%sd(:), ' )'
    End Select

    nnumbers = 0_8
    if (rnd%cyc > 0_8) nnumbers = 96_8*(rnd%cyc-1_8)
    if (rnd%cnt <= 96) nnumbers = nnumbers + rnd%cnt - 1_8
    Write(iout,'(1A,1I22,1A)')'RANLUX since seeded generated: ', &
         & nnumbers*48_8,  ' bits' 
    
    Return
  End Subroutine rnlx48_print_info

! ********************************
! *
  Subroutine rnlx48_error(routine, msg)
! *
! ********************************
    
    Character (len=*), Intent (in) :: routine, msg
    
    Write(error_unit,*)'In '//Trim(routine)// ' :'
    Write(error_unit,'(5X,1A)')Trim(msg)
    Write(error_unit,*)
    
    Stop
  End Subroutine rnlx48_error
  
End MODULE RanLux48

