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
MODULE RanLux
! *
! ***************************************************

  USE ISO_FORTRAN_ENV, Only : error_unit, output_unit
  IMPLICIT NONE

  private

  type :: state
     integer :: V(96)
     integer :: c(4)
     integer :: cnt, pr, ir, is
     integer :: ilv, otherseed = -1, sd(4) = -1, SEED_TYPE = -1
     integer (kind=8) :: cyc = 0_8
  end type state

  integer, parameter :: BASE = 16777216
  real (kind=8), parameter :: DINV = 5.9604644775390625E-8_8

  type (state), save :: rnd
  integer :: next(1:24) = (/24,1,2,3,4,5,6,7,8,9,10,11,12,13,14, &
       15,16,17,18,19,20,21,22,23/)

  Integer, Parameter, Private :: SEED_LCG=0, SEED_MINE=1, &
       & SEED_V3=2

  interface rnlx
     module procedure :: rnlx_f64, rnlx_f64_v, rnlx_f32, rnlx_f32_v, &
          rnlx_z32, rnlx_z32_v, rnlx_z64, rnlx_z64_v
  end interface rnlx

  public :: rnlx_init, rnlx, rnlx_seed_lcg, rnlx_seed_v3, &
       rnlx_seed_mine, rnlx_print_info, rnlx_size, rnlx_get, rnlx_reset, &
       rnlx_f32, rnlx_f32_v, rnlx_f64, rnlx_f64_v, &
       rnlx_z32, rnlx_z32_v, rnlx_z64, rnlx_z64_v

  !$OMP THREADPRIVATE(rnd)

CONTAINS

! ***************************************************
! *
  Function rnlx_size()
! *
! ***************************************************
    integer :: rnlx_size

    rnlx_size = 112

    return
  end Function rnlx_size

! ***************************************************
! *
  Subroutine rnlx_get(state)
! *
! ***************************************************
    
    Integer (kind=8), intent (out) :: state(:)

    if (size(state)<rnlx_size()) call rnlx_error('rnlx_get', &
         'Array too small to store the state of RANLUX')
    
    state(1:96)    = int(rnd%V(1:96)  ,kind=8)
    state(97:100)  = int(rnd%c(1:4)   ,kind=8)
    state(101)     = int(rnd%cnt      ,kind=8)
    state(102)     = int(rnd%pr       ,kind=8)
    state(103)     = int(rnd%ir       ,kind=8)
    state(104)     = int(rnd%is       ,kind=8)

    state(105)     = int(rnd%ilv      ,kind=8)
    state(106)     = int(rnd%otherseed,kind=8)
    state(107:110) = int(rnd%sd(1:4)  ,kind=8)
    state(111)     = int(rnd%SEED_TYPE,kind=8)
    state(112)     = int(rnd%cyc      ,kind=8)

    Return
  End Subroutine rnlx_get

! ***************************************************
! *
  Subroutine rnlx_reset(state)
! *
! ***************************************************
    
    Integer (kind=8), intent (in) :: state(:)

    if (size(state)<rnlx_size()) call rnlx_error('rnlx_reset', &
         'Array too small to store the state of RANLUX')
    
    rnd%V(1:96)    = int(state(1:96)   ,kind=4)  
    rnd%c(1:4)     = int(state(97:100) ,kind=4)
    rnd%cnt        = int(state(101)    ,kind=4)
    rnd%pr         = int(state(102)    ,kind=4)
    rnd%ir         = int(state(103)    ,kind=4)
    rnd%is         = int(state(104)    ,kind=4)

    rnd%ilv        = int(state(105)    ,kind=4)
    rnd%otherseed  = int(state(106)    ,kind=4)
    rnd%sd(1:4)    = int(state(107:110),kind=4)
    rnd%SEED_TYPE  = int(state(111)    ,kind=4)
    rnd%cyc        = int(state(112)    ,kind=8)

    Return
  End Subroutine rnlx_reset


! ***************************************************
! *
  Subroutine rnlx_seed_lcg(nseed)
! *
! ***************************************************
    
    Integer, intent (in) :: nseed
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
  End Subroutine rnlx_seed_lcg


! ***************************************************
! *
  Subroutine rnlx_seed_v3(nseed)
! *
! ***************************************************
    
    Integer, intent (in) :: nseed
    
    integer :: ibt, jbt, bits(0:30), ns, ib, ix, i, j, k

    if ((nseed>2147483647).or.(nseed<1)) call rnlx_error(&
         'rnlx_seed_mine',&
         'First component of seed too large (1<seed(1)<2147483647)')

    ns = nseed
    do i = 0, 30
       bits(i) = mod(ns,2)
       ns = ns/2
    end do
    
    ibt = 0
    jbt = 18

    do k = 0, 3
       do j = 0, 23
          ix = 0
          do i = 0, 23
             ib = bits(ibt)
             ix = 2*ix+ib
             
             bits(ibt) = mod(bits(ibt)+bits(jbt),2)
             ibt = mod(ibt+1,31)
             jbt = mod(jbt+1,31)
          end do
          if (mod(j,4)==k) ix = BASE-1-ix
          rnd%V(24-j+k*24) = ix
       end do
    end do
    rnd%cnt = 97
    rnd%SEED_TYPE = SEED_V3
    rnd%otherseed = nseed
    rnd%cyc = 0_8
    
    Return
  End Subroutine rnlx_seed_v3


! ***************************************************
! *
  Subroutine rnlx_seed_mine(nseed)
! *
! ***************************************************
    Integer, intent (in) :: nseed(:)
    
    integer :: ibt, jbt, bits(0:30), ns(4), ib, ix, i, j, k, &
         isz, jj, ic

    isz = size(nseed)
    ns(1) = nseed(1)

    call checkseed()
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
       
       do k = 0, 3
          do jj = 0, 5
             j = jj+6*(ic-1)
             ix = 0
             do i = 0, 23
                ib = bits(ibt)
                ix = 2*ix+ib
                
                bits(ibt) = mod(bits(ibt)+bits(jbt),2)
                ibt = mod(ibt+1,31)
                jbt = mod(jbt+1,31)
             end do
             if (mod(j,4)==k) ix = BASE-1-ix
             rnd%V(24-j+k*24) = ix
          end do
       end do
       
    end do

    rnd%cnt = 97

    Return
  contains
    
    subroutine checkseed()

      if (isz>4) call rnlx_error('rnlx_seed_mine','Too many numbers in seed')
      if (nseed(1)>2147483647) call rnlx_error('rnlx_seed_mine',&
           'First component of seed too large (1<seed(1)<2147483647)')
      if ((isz>1)) then 
         if ((nseed(2)<2).or.(nseed(2)>2147483647)) call rnlx_error('rnlx_seed_mine',&
              'Second component of seed out of range (1<seed(2)<2147483647)')
         ns(2) = nseed(2)
      else
         ns(2) = 1
      end if
      if ((isz>2)) then 
         if ((nseed(3)<2).or.(nseed(3)>2147483647)) call rnlx_error('rnlx_seed_mine',&
              'Third component of seed out of range (1<seed(3)<2147483647)')
         ns(3) = nseed(3)
      else
         ns(3) = 1
      end if
      if ((isz>3)) then 
         if ((nseed(4)<2).or.(nseed(4)>2147483647)) call rnlx_error('rnlx_seed_mine',&
              'Fourth component of seed out of range (1<seed(4)<2147483647)')
         ns(4) = nseed(4)
      else
         ns(4) = 1
      end if

      return
    end subroutine checkseed

  End Subroutine rnlx_seed_mine

! ***************************************************
! *
  Subroutine rnlx_init(n)
! *
! ***************************************************
    integer, intent (in), optional :: n
    integer :: np

    rnd%is  = 10
    rnd%ir  = 24
    rnd%c   = 0
    rnd%cnt = 1

    if (present(n)) then
       if (n>24) then
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
       rnd%pr = 24
       rnd%ilv = -1
    case (0)
       rnd%pr = 218
       rnd%ilv = 0
    case (1)
       rnd%pr = 404
       rnd%ilv = 1
    case (2)
       rnd%pr = 794
       rnd%ilv = 2
    end select

    return
  end Subroutine rnlx_init


! ***************************************************
! *
  Subroutine fill_rnd()
! *
! ***************************************************
    integer :: delta(4), j, ic

    do ic = 1, rnd%pr
       do concurrent (j=0:3) 
          delta(j+1) = rnd%V(rnd%is+24*j) - rnd%V(rnd%ir+24*j) - rnd%c(j+1)
       end do

       where (delta < 0)
          rnd%c = 1
          delta = delta + BASE
       elsewhere
          rnd%c = 0
       end where

       do concurrent (j=0:3) 
          rnd%V(rnd%ir+24*j) = delta(j+1)
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
  Subroutine rnlx_int(n)
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
  End Subroutine rnlx_int

! ***************************************************
! *
  Subroutine rnlx_int_v(d)
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
  End Subroutine rnlx_int_v

! ***************************************************
! *
  Subroutine rnlx_f64(d)
! *
! ***************************************************

    Real (kind=8), Intent (out) :: d
    Integer (kind=4) :: ipos, i1, i2
    
    ipos = rnd%cnt
    If (ipos > 96) Then
       CALL fill_rnd()
       i1 = rnd%V(1)
       i2 = rnd%V(2)
       rnd%cnt = rnd%cnt + 2
    else if (ipos == 96) then
       i1 = rnd%V(96)
       CALL fill_rnd()
       i2 = rnd%V(1)
       rnd%cnt = rnd%cnt + 1
    else
       i1 = rnd%V(ipos)
       i2 = rnd%V(ipos+1)
       rnd%cnt = rnd%cnt + 2
    End If
    d = (Real(i1,kind=8)*DINV + Real(i2,kind=8))*DINV

    Return
  End Subroutine rnlx_f64

! ***************************************************
! *
  Subroutine rnlx_f32(d)
! *
! ***************************************************

    Real (kind=4), Intent (out) :: d
    Integer (kind=4) :: ipos, i1
    
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
  End Subroutine rnlx_f32

! ***************************************************
! *
  Subroutine rnlx_z32(d)
! *
! ***************************************************

    Complex (kind=4), Intent (out) :: d
    Real (kind=4) :: x, y

    CALL rnlx(x)
    CALL rnlx(y)
    d = Cmplx(x,y)
    
    Return
  End Subroutine rnlx_z32

! ***************************************************
! *
  Subroutine rnlx_z64(d)
! *
! ***************************************************

    Complex (kind=8), Intent (out) :: d
    Real (kind=8) :: x, y

    CALL rnlx(x)
    CALL rnlx(y)
    d = Cmplx(x,y,kind=8)
    
    Return
  End Subroutine rnlx_z64

! ***************************************************
! *
  Subroutine rnlx_z32_v(d)
! *
! ***************************************************

    Complex (kind=4), Intent (out) :: d(:)
    Real (kind=4) :: x(Size(d))
    
    CALL rnlx(x)
    d = Cmplx(x,0.0_4)
    CALL rnlx(x)
    d = d + Cmplx(0.0_4,x)

    Return
  End Subroutine rnlx_z32_v

! ***************************************************
! *
  Subroutine rnlx_z64_v(d)
! *
! ***************************************************

    Complex (kind=8), Intent (out) :: d(:)
    Real (kind=8) :: x(Size(d))
    
    CALL rnlx(x)
    d = Cmplx(x,0.0_8,kind=8)
    CALL rnlx(x)
    d = d + Cmplx(0.0_8,x,kind=8)

    Return
  End Subroutine rnlx_z64_v

! ***************************************************
! *
  Subroutine rnlx_f64_v(d)
! *
! ***************************************************

    Real (kind=8), Intent (out) :: d(:)
    Integer (kind=4) :: i, nleft, nst, ncyc, iex
    
    nleft = Size(d)
    nst = int( (96-rnd%cnt+1)/2 )
    iex = mod( (96-rnd%cnt+1),2 )
    If (nst > nleft) Then
       d(1:nleft) = (Real(rnd%V(rnd%cnt:rnd%cnt+nleft-1),kind=8)&
            & * DINV + &
            Real(rnd%V(rnd%cnt+nleft:rnd%cnt+2*nleft-1),kind=8)) * DINV
       rnd%cnt = rnd%cnt+2*nleft
       Return
    else if (nst>0) then
       d(1:nst) = (Real(rnd%V(rnd%cnt:rnd%cnt+nst-1),kind=8) &
            & * DINV + &
            Real(rnd%V(rnd%cnt+nst:rnd%cnt+2*nst-1),kind=8)) * DINV
       nleft = nleft - nst
    else
       nst = 0
    End If
    CALL fill_rnd() 
    
    ncyc = Int(nleft/48)
    Do I = 1, ncyc
       d(nst+1+(I-1)*48:nst+I*48) = &
            & ( Real(rnd%V(1:48),kind=8) * DINV + &
            Real(rnd%V(49:96),kind=8)  ) * DINV
       CALL fill_rnd()
    End Do
    nleft = Modulo(nleft,48)
    d(nst+1+ncyc*48:Size(d)) = &
         & (Real(rnd%V(1:nleft),kind=8) * DINV + &
         Real(rnd%V(nleft+1:2*nleft),kind=8) ) * DINV
    rnd%cnt = 2*nleft+1

    Return
  End Subroutine rnlx_f64_v

! ***************************************************
! *
  Subroutine rnlx_f32_v(d)
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
  End Subroutine rnlx_f32_v

! ***************************************************
! *
  Subroutine rnlx_print_info(ifn)
! *
! ***************************************************

    Integer, Intent (in), Optional :: ifn
    Integer :: iout
    integer (kind=8) :: nnumbers
    
    iout = output_unit
    if (present(ifn)) iout = ifn

    Write(iout,'(1A,1I4)')'RANLUX Working with luxury level: ', rnd%ilv
    Write(iout,'(1A,1I6)')'RANLUX luxury parameter:  ', rnd%pr
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
         & nnumbers*24_8,  ' bits' 
    
    Return
  End Subroutine rnlx_print_info

! ********************************
! *
    Subroutine rnlx_error(routine, msg)
! *
! ********************************

      Character (len=*), Intent (in) :: routine, msg

      Write(error_unit,*)'In '//Trim(routine)// ' :'
      Write(error_unit,'(5X,1A)')Trim(msg)
      Write(error_unit,*)

      Stop
    End Subroutine rnlx_error

End MODULE RanLux

