!
! A MODULE for RANDOM numbers
!
! "THE BEER-WARE LICENSE":
! Alberto Ramos wrote this file. As long as you retain this 
! notice you can do whatever you want with this stuff. If we meet some 
! day, and you think this stuff is worth it, you can buy me a beer in 
! return. <alberto.ramos@desy.de>
!

! ***************************************************
! *
MODULE random
! *
! ***************************************************

  USE ISO_FORTRAN_ENV, Only : error_unit, output_unit
  USE NumTypes
  USE Constants
  USE RANLUX48
  USE MIXMAX

  IMPLICIT NONE

  private
  
  interface
     subroutine f32(r)
       USE Numtypes
       real (kind=SP), intent (out) :: r
     end subroutine f32
     subroutine f32_v(r)
       USE Numtypes
       real (kind=SP), intent (out) :: r(:)
     end subroutine f32_v
     subroutine f64(r)
       USE Numtypes
       real (kind=DP), intent (out) :: r
     end subroutine f64
     subroutine f64_v(r)
       USE Numtypes
       real (kind=DP), intent (out) :: r(:)
     end subroutine f64_v

     subroutine z32(r)
       USE Numtypes
       complex (kind=SP), intent (out) :: r
     end subroutine z32
     subroutine z32_v(r)
       USE Numtypes
       complex (kind=SP), intent (out) :: r(:)
     end subroutine z32_v
     subroutine z64(r)
       USE Numtypes
       complex (kind=DP), intent (out) :: r
     end subroutine z64
     subroutine z64_v(r)
       USE Numtypes
       complex (kind=DP), intent (out) :: r(:)
     end subroutine z64_v

     subroutine seed(ns)
       integer, intent (in) :: ns(:)
     end subroutine seed
     subroutine get(ns)
       integer (kind=8), intent (out) :: ns(:)
     end subroutine get
     subroutine reset(ns)
       integer (kind=8), intent (in) :: ns(:)
     end subroutine reset

     subroutine info(ifn)
       integer, intent (in), optional :: ifn
     end subroutine info
     function rsize()
       integer :: rsize
     end function rsize
  end interface

  procedure (f32),   pointer :: rndm_f32   => null()
  procedure (f32_v), pointer :: rndm_f32_v => null()
  procedure (f64),   pointer :: rndm_f64   => null()
  procedure (f64_v), pointer :: rndm_f64_v => null()
  procedure (z32),   pointer :: rndm_z32   => null()
  procedure (z32_v), pointer :: rndm_z32_v => null()
  procedure (z64),   pointer :: rndm_z64   => null()
  procedure (z64_v), pointer :: rndm_z64_v => null()

  procedure (seed),  pointer :: rndm_seed  => null()
  procedure (rsize), pointer :: rndm_size  => null()

  procedure (info),  pointer :: rndm_print_info => null()

  interface rndm
     procedure rndm_f32, rndm_f32_v, rndm_f64, rndm_f64_v, &
          rndm_z32, rndm_z32_v, rndm_z64, rndm_z64_v
  end interface

  public :: rndm, rnlx_start, mxmx_start, rndm_seed, rndm_get, &
       rndm_reset, rndm_size, rndm_print_info

  interface normal 
     Module Procedure NormalS, NormalV, NormalS2, NormalV2, &
          & NormalS_SP, NormalV_SP, NormalS2_SP, NormalV2_SP
  end interface

  interface laplace
     Module Procedure Laplace_DP, Laplace_SP,Laplace2_DP, Laplace2_SP
  end interface

  interface levy
     Module Procedure Levy_DP, Levy_SP, LevyV_DP, LevyV_SP
  end interface levy

  interface cauchy
     Module Procedure Cauchy_DP, Cauchy_SP, CauchyV_DP, CauchyV_SP
  end interface

  interface lorentz
     Module Procedure  Lorentz_DP, Lorentz_SP, LorentzV_DP, LorentzV_SP
  end interface

  interface fishtipp
     Module Procedure FishTipp_DP, FishTipp_SP,FishTipp2_DP, FishTipp2_SP
  end interface


  integer (kind=8), parameter :: IS_RNLX=0, IS_MXMX=1
  integer (kind=8) :: IRNDM

  public :: normal, laplace, levy, cauchy, lorentz, fishtipp, metropolis

CONTAINS


! ***************************************************
! *
  Subroutine mxmx_start(nin)
! *
! ***************************************************

    Integer, Intent (in), Optional :: nin
    Integer :: nmat, ilv

    If (Present(nin)) Then
       nmat = nin
    Else
       nmat = 256
    End If

    ilv = 1
    call mxmx_init(nmat)
    call mxmx_setlux(ilv)
    if (ilv > 0) then
       rndm_f32   => mxmx2_f32
       rndm_f32_v => mxmx2_f32_v
       rndm_f64   => mxmx2_f64
       rndm_f64_v => mxmx2_f64_v
       rndm_z32   => mxmx2_z32
       rndm_z32_v => mxmx2_z32_v
       rndm_z64   => mxmx2_z64
       rndm_z64_v => mxmx2_z64_v
    else
       rndm_f32   => mxmx_f32
       rndm_f32_v => mxmx_f32_v
       rndm_f64   => mxmx_f64
       rndm_f64_v => mxmx_f64_v
       rndm_z32   => mxmx_z32
       rndm_z32_v => mxmx_z32_v
       rndm_z64   => mxmx_z64
       rndm_z64_v => mxmx_z64_v
    end if

    rndm_seed  => mxmx_seed_skip
    rndm_size  => mxmx_size

    rndm_print_info => mxmx_print_info
    
    IRNDM = IS_MXMX

    return
  end Subroutine mxmx_start

! ***************************************************
! *
  Subroutine rnlx_start(nin)
! *
! ***************************************************

    Integer, Intent (in), Optional :: nin
    Integer :: nmat

    If (Present(nin)) Then
       nmat = nin
    Else
       nmat = 1
    End If

    call rnlx48_init(nmat)
    rndm_f32   => rnlx48_f32
    rndm_f32_v => rnlx48_f32_v
    rndm_f64   => rnlx48_f64
    rndm_f64_v => rnlx48_f64_v
    rndm_z32   => rnlx48_z32
    rndm_z32_v => rnlx48_z32_v
    rndm_z64   => rnlx48_z64
    rndm_z64_v => rnlx48_z64_v

    rndm_seed  => rnlx48_seed_mine
    rndm_size  => rnlx48_size

    rndm_print_info => rnlx48_print_info

    IRNDM = IS_RNLX

    return
  end Subroutine rnlx_start

! ***************************************************
! *
  Subroutine rndm_get(iseed)
! *
! ***************************************************
    integer (kind=8), intent (out) :: iseed(0:)

    iseed(0) = IRNDM
    if (iseed(0)==IS_RNLX) then
       call rnlx48_get(iseed(1:))
    else if (iseed(0)==IS_MXMX) then
       call mxmx_get(iseed(1:))
    end if

    return
  end Subroutine rndm_get

! ***************************************************
! *
  Subroutine rndm_reset(iseed)
! *
! ***************************************************
    integer (kind=8), intent (in) :: iseed(0:)

    if (iseed(0)==IS_RNLX) then
       call rnlx48_reset(iseed(1:))
    else if (iseed(0)==IS_MXMX) then
       call mxmx_reset(iseed(1:))
    end if

    return
  end Subroutine rndm_reset

!  *********************************************
!  *                                           *
  subroutine normals(X)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X
    Real (kind=DP) :: U1, U2
    
    call rndm(u1)
    call rndm(u2)

    x = sqrt(-2.0_dp*log(1.0_DP-u1)) * cos(twopi_dp*u2)

    return
  end subroutine  normals

!  *********************************************
!  *                                           *
  subroutine normalv(X)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP) :: U1, U2
    
    integer :: i, isz
    complex (kind=DPC) :: z
    real (kind=DP)     :: r

    isz = int(size(X)/2)
    do i = 1, isz
       CALL Rndm(U1)
       CALL Rndm(U2)

       r = Sqrt(-2.0_DP*Log(1.0_DP-U1))
       z = exp(cmplx(0.0_DP,TWOPI_DP,kind=DP)*U2)
       X(2*i-1) = r*Real(z)
       X(2*i)   = r*Aimag(z)
    end do
    
    if (mod(size(X),2) == 1) call normals(X(size(X)))

    Return
  End Subroutine  NormalV

!  *********************************************
!  *                                           *
  subroutine normals2(X, Rmed, Rsig)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent (in) :: Rmed, Rsig
    Real (kind=DP) :: U1, U2
    
    CALL Rndm(U1)
    CALL Rndm(U2)

    X = Rsig * Sqrt(-2.0_DP*Log(1.0_DP-U1)) * Cos(TWOPI_DP*U2) + Rmed 

    Return
  End Subroutine  NormalS2

!  *********************************************
!  *                                           *
  subroutine normalv2(X, rmed, rsig)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent (in) :: Rmed, Rsig
    Real (kind=DP) :: U1, U2
    
    integer :: i, isz
    complex (kind=DPC) :: z
    real (kind=DP)     :: r

    isz = int(size(X)/2)
    do i = 1, isz
       CALL Rndm(U1)
       CALL Rndm(U2)

       r = Sqrt(-2.0_DP*Log(1.0_DP-U1))
       z = exp(cmplx(0.0_DP,TWOPI_DP,kind=DP)*U2)
       X(2*i-1) = r*Real(z)
       X(2*i)   = r*Aimag(z)
    end do
    
    if (mod(size(X),2) == 1) call normals(X(size(X)))

    X = rsig*X + Rmed

    return
  end subroutine  normalv2

!  *********************************************
!  *                                           *
  subroutine normals_SP(X)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X
    Real (kind=SP) :: U1, U2
    
    call rndm(u1)
    call rndm(u2)

    x = sqrt(-2.0_SP*log(1.0_SP-u1)) * cos(twopi_SP*u2)

    return
  end subroutine  normals_SP

!  *********************************************
!  *                                           *
  subroutine normalv_SP(X)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP) :: U1, U2
    
    integer :: i, isz
    complex (kind=SPC) :: z
    real (kind=SP)     :: r

    isz = int(size(X)/2)
    do i = 1, isz
       CALL Rndm(U1)
       CALL Rndm(U2)

       r = Sqrt(-2.0_SP*Log(1.0_SP-U1))
       z = exp(cmplx(0.0_SP,TWOPI_SP)*U2)
       X(2*i-1) = r*Real(z)
       X(2*i)   = r*Aimag(z)
    end do
    
    if (mod(size(X),2) == 1) call normals_SP(X(size(X)))

    Return
  End Subroutine  Normalv_SP

!  *********************************************
!  *                                           *
  subroutine normals2_SP(X, Rmed, Rsig)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent (in) :: Rmed, Rsig
    Real (kind=SP) :: U1, U2
    
    CALL Rndm(U1)
    CALL Rndm(U2)

    X = Rsig * Sqrt(-2.0_SP*Log(1.0_SP-U1)) * Cos(TWOPI_SP*U2) + Rmed 

    Return
  End Subroutine  Normals2_SP

!  *********************************************
!  *                                           *
  subroutine normalv2_SP(X, rmed, rsig)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent (in) :: Rmed, Rsig
    Real (kind=SP) :: U1, U2
    
    integer :: i, isz
    complex (kind=SPC) :: z
    real (kind=SP)     :: r

    isz = int(size(X)/2)
    do i = 1, isz
       CALL Rndm(U1)
       CALL Rndm(U2)

       r = Sqrt(-2.0_SP*Log(1.0_SP-U1))
       z = exp(cmplx(0.0_SP,TWOPI_SP)*U2)
       X(2*i-1) = r*Real(z)
       X(2*i)   = r*Aimag(z)
    end do
    
    if (mod(size(X),2) == 1) call normals_SP(X(size(X)))
    X = rsig*X + Rmed

    return
  end subroutine  normalv2_SP

!  *********************************************
!  *                                           *
  Subroutine Laplace_SP(X, Rmu, Rb)
!  *                                           *
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
    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent (in) :: Rmu, Rb
    Real (kind=SP) :: U

    integer :: i
    
    do i = 1, size(X)
       CALL Random_Number(U)
       U = U - 0.5_SP
       X(i) = Rmu - Rb * Sign(Log(1-2.0_SP*Abs(U)), U)
    end do

    Return
  End Subroutine  Laplace2_SP

!  *********************************************
!  *                                           *
  Subroutine Laplace2_DP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent (in) :: Rmu, Rb
    Real (kind=DP) :: U
    
    integer :: i
    
    do i = 1, size(X)
       CALL Random_Number(U)
       U = U - 0.5_SP
       X(i) = Rmu - Rb * Sign(Log(1-2.0_SP*Abs(U)), U)
    end do

    Return
  End Subroutine  Laplace2_DP

!  *********************************************
!  *                                           *
  Subroutine Levy_DP(X, a, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent(in)  :: a
    Real (kind=DP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=DP) :: b, sig, mu, w, phi, ainv, z
    Real (kind=DP) :: U1, U2
    
    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_DP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_DP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_DP
    End If
   
    CALL rndm(U1)
    CALL rndm(U2)
    
    w = -Log(1.0_DP-U1)
    phi = PI_DP * (U2 - 0.5_DP)
    z = b*Tan(a*HALFPI_DP)
    ainv = 1.0_DP/a

    X = (Sin(a*phi) + z*Cos(a*phi))/Cos(phi) * &
         & ( (Cos((1.0_DP-a)*phi) + z*Sin((1.0_DP-a)*phi))/(w*Cos(phi)) )**(ainv-1.0_DP)
    X = sig*X + mu

    Return
  End Subroutine  Levy_DP

!  *********************************************
!  *                                           *
  Subroutine LevyV_DP(X, a, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent(in)  :: a
    Real (kind=DP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=DP) :: b, sig, mu, w, phi, ainv, z
    Real (kind=DP) :: U1, U2
    Integer :: I
  
    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_DP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_DP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_DP
    End If
   
    z = b*Tan(a*HALFPI_DP)
    ainv = 1.0_DP/a

    Do I = 1, Size(X)
       CALL rndm(U1)
       CALL rndm(U2)
       
       w = -Log(1.0_DP-U1)
       phi = PI_DP * (U2 - 0.5_DP)
       
       X(I) = (Sin(a*phi) + z*Cos(a*phi))/Cos(phi) * &
            & ( (Cos((1.0_DP-a)*phi) + &
            &    z*Sin((1.0_DP-a)*phi))/(w*Cos(phi)) )**(ainv-1.0_DP)

    End Do

    X = sig*X + mu

    Return
  End Subroutine  LevyV_DP
  
!  *********************************************
!  *                                           *
  Subroutine Levy_SP(X, a, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent(in)  :: a
    Real (kind=SP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=SP) :: b, sig, mu, w, phi, ainv, z
    Real (kind=SP) :: U1, U2
    
    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_SP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_SP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_SP
    End If
   
    CALL rndm(U1)
    CALL rndm(U2)
    
    w = -Log(1.0_SP-U1)
    phi = PI_SP * (U2 - 0.5_SP)
    z = b*Tan(a*HALFPI_SP)
    ainv = 1.0_SP/a

    X = (Sin(a*phi) + z*Cos(a*phi))/Cos(phi) * &
         & ( (Cos((1.0_SP-a)*phi) + z*Sin((1.0_SP-a)*phi))/(w*Cos(phi)) )**(ainv-1.0_SP)
    X = sig*X + mu

    Return
  End Subroutine  Levy_SP

!  *********************************************
!  *                                           *
  Subroutine LevyV_SP(X, a, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent(in)  :: a
    Real (kind=SP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=SP) :: b, sig, mu, w, phi, z, ainv
    Real (kind=SP) :: U1, U2
    Integer :: I
  
    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_SP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_SP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_SP
    End If
   
    z = b*Tan(a*HALFPI_SP)
    ainv = 1.0_SP/a

    Do I = 1, Size(X)
       CALL rndm(U1)
       CALL rndm(U2)
       
       w = -Log(1.0_SP-U1)
       phi = PI_SP * (U2 - 0.5_SP)
       
       X(I) = (Sin(a*phi) + z*Cos(a*phi))/Cos(phi) * &
            & ( (Cos((1.0_SP-a)*phi) + &
            &    z*Sin((1.0_SP-a)*phi))/(w*Cos(phi)) )**(ainv-1.0_SP)

    End Do

    X = sig*X + mu

    Return
  End Subroutine  LevyV_SP

!  *********************************************
!  *                                           *
  Subroutine Cauchy_DP(X, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=DP) :: b, sig, mu, w, phi
    Real (kind=DP) :: U1, U2
    
    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_DP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_DP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_DP
    End If
   
    CALL rndm(U1)
    CALL rndm(U2)
    
    w = -Log(1.0_DP-U1)
    phi = PI_DP * (U2 - 0.5_DP)

    X = 1.0_DP/HALFPI_DP * ( (HALFPI_DP+b*phi)*Tan(phi) - &
         & b*log( (HALFPI_DP*w*cos(phi))/(HALFPI_DP + b*phi) ) )

    X = sig*X + mu

    Return
  End Subroutine Cauchy_DP

!  *********************************************
!  *                                           *
  Subroutine CauchyV_DP(X, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=DP) :: b, sig, mu, w, phi
    Real (kind=DP) :: U1, U2
    Integer :: I

    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_DP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_DP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_DP
    End If
   
    Do I = 1, Size(X)
       CALL rndm(U1)
       CALL rndm(U2)
       
       w = -Log(1.0_DP-U1)
       phi = PI_DP * (U2 - 0.5_DP)
       
       X(I) = 1.0_DP/HALFPI_DP * ( (HALFPI_DP+b*phi)*Tan(phi) - &
            & b*log( (HALFPI_DP*w*cos(phi))/(HALFPI_DP + b*phi) ) )
       
    End Do
    X = sig*X + mu

    Return
  End Subroutine CauchyV_DP

!  *********************************************
!  *                                           *
  Subroutine Cauchy_SP(X, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=SP) :: b, sig, mu, w, phi
    Real (kind=SP) :: U1, U2
    
    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_SP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_SP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_SP
    End If
   
    CALL rndm(U1)
    CALL rndm(U2)
    
    w = -Log(1.0_SP-U1)
    phi = PI_SP * (U2 - 0.5_SP)

    X = 1.0_SP/HALFPI_SP * ( (HALFPI_SP+b*phi)*Tan(phi) - &
         & b*log( (HALFPI_SP*w*cos(phi))/(HALFPI_SP + b*phi) ) )

    X = sig*X + mu

    Return
  End Subroutine Cauchy_SP

!  *********************************************
!  *                                           *
  Subroutine CauchyV_SP(X, bi, sigi, mui)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent(in), Optional  :: bi, sigi, mui
    
    Real (kind=SP) :: b, sig, mu, w, phi
    Real (kind=SP) :: U1, U2
    Integer :: I

    If (Present(bi)) Then
       b = bi
    Else
       b = 0.0_SP
    End If

    If (Present(sigi)) Then
       sig = sigi
    Else
       sig = 1.0_SP
    End If

    If (Present(bi)) Then
       mu = mui
    Else
       mu = 0.0_SP
    End If
   
    Do I = 1, Size(X)
       CALL rndm(U1)
       CALL rndm(U2)
       
       w = -Log(1.0_SP-U1)
       phi = PI_SP * (U2 - 0.5_SP)
       
       X(I) = 1.0_SP/HALFPI_SP * ( (HALFPI_SP+b*phi)*Tan(phi) - &
            & b*log( (HALFPI_SP*w*cos(phi))/(HALFPI_SP + b*phi) ) )
       
    End Do
    X = sig*X + mu

    Return
  End Subroutine CauchyV_SP

!  *********************************************
!  *                                           *
  Subroutine Lorentz_DP(X, mu, gamma)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent(in), Optional  :: mu, gamma
    
    Real (kind=DP) :: U
    
    CALL rndm(U)
    
    X = Tan(TWOPI_DP*U)

    If (Present(gamma)) X = gamma * X 
    If (Present(mu))    X = X + mu

    Return
  End Subroutine Lorentz_DP

!  *********************************************
!  *                                           *
  Subroutine LorentzV_DP(X, mu, gamma)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent(in), Optional  :: mu, gamma
    
    Real (kind=DP) :: U
    Integer :: I

    Do I = 1, Size(X)
       CALL rndm(U)
       X(I) = Tan(TWOPI_DP*U)
    End Do
    
    If (Present(gamma)) X(:) = gamma * X(:) 
    If (Present(mu))    X(:) = X(:) + mu

    Return
  End Subroutine LorentzV_DP

!  *********************************************
!  *                                           *
  Subroutine Lorentz_SP(X, mu, gamma)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent(in), Optional  :: mu, gamma
    
    Real (kind=SP) :: U
    
    CALL rndm(U)
    
    X = Tan(TWOPI_SP*U)

    If (Present(gamma)) X = gamma * X 
    If (Present(mu))    X = X + mu


    Return
  End Subroutine Lorentz_SP

!  *********************************************
!  *                                           *
  Subroutine LorentzV_SP(X, mu, gamma)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent(in), Optional  :: mu, gamma
    
    Real (kind=SP) :: U
    Integer :: I

    Do I = 1, Size(X)
       CALL rndm(U)
       X(I) = Tan(TWOPI_SP*U)
    End Do
    
    If (Present(gamma)) X(:) = gamma * X(:) 
    If (Present(mu))    X(:) = X(:) + mu

    Return
  End Subroutine LorentzV_SP

!  *********************************************
!  *                                           *
  subroutine fishtipp_SP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X
    Real (kind=SP), Intent (in) :: Rmu, Rb
    Real (kind=SP) :: U
    
    CALL rndm(U)

    X = Rmu - Rb * Log(-Log(1.0_SP-U))

    Return
  End Subroutine  FishTipp_SP

!  *********************************************
!  *                                           *
  Subroutine FishTipp_DP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X
    Real (kind=DP), Intent (in) :: Rmu, Rb
    Real (kind=DP) :: U
    
    CALL rndm(U)

    X = Rmu - Rb * Log(-Log(1.0_DP-U))

    Return
  End Subroutine  FishTipp_DP

!  *********************************************
!  *                                           *
  Subroutine FishTipp2_SP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
    Real (kind=SP), Intent(out) :: X(:)
    Real (kind=SP), Intent (in) :: Rmu, Rb
    Real (kind=SP) :: U
    
    integer :: i

    do i = 1, size(X)
       CALL rndm(U)
       X(i) = Rmu - Rb * Log(-Log(1.0_SP-U))
    end do

    Return
  End Subroutine  FishTipp2_SP

!  *********************************************
!  *                                           *
  Subroutine FishTipp2_DP(X, Rmu, Rb)
!  *                                           *
!  *********************************************
    Real (kind=DP), Intent(out) :: X(:)
    Real (kind=DP), Intent (in) :: Rmu, Rb
    Real (kind=DP) :: U
    
    integer :: i

    do i = 1, size(X)
       CALL rndm(U)
       X(i) = Rmu - Rb * Log(-Log(1.0_DP-U))
    end do

    Return
  End Subroutine  FishTipp2_DP

!  *********************************************
!  *                                           *
  function metropolis(Pfin,Pini)
!  *                                           *
!  *********************************************
    real (kind=DP), intent (in) :: pfin, pini
    logical :: metropolis

    real (kind=DP) :: Pacc, r

    Pacc = Pfin / Pini
    if (pacc > 1.0_DP) then
       metropolis = .true.
       return
    else
       call rndm(r)
       if (pacc > r) then
          metropolis = .true.
          return
       end if
    end if
    metropolis = .false.

    return
  end function metropolis

end MODULE Random
