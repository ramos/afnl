
MODULE LapackAPI

!  USE MODMinuit
  USE NumTypes
  USE ISO_C_BINDING

  IMPLICIT NONE

  Interface GEVP
     Module Procedure GEVP_D, GEVP_C
  End Interface

  Interface SVD
     Module Procedure SVD, SVD_C
  End Interface

  Interface PseudoInverse
     Module Procedure PseudoInverse, PseudoInverse_C
  End Interface

  Interface Diag
     Module Procedure Diag_DPC
  End Interface Diag

  Interface expM
     Module Procedure expM_C
  End Interface ExpM

  Private GEVP_D, GEVP_C, expM_C

CONTAINS

! ***********************************************
! *
  Subroutine GEVP_D(A, B, Dv, U, Info) 
! *
! ***********************************************
! * Solves the Generalized Eigenvalue Problem 
! *   Av = lBv
! * "eigenvalues" returned in Dv(:)
! *  
! ***********************************************

    Real (kind=DP), Intent (in) :: A(:,:), B(:,:)
    Real (kind=DP), Intent (out), Optional :: U(:,:)
    Real (kind=DP), Intent (out) :: Dv(:)
    Integer, Intent (out), Optional :: Info

!    Real (kind=DP) :: SA(Size(A,1),Size(A,2))
    Real (kind=DP), Allocatable :: Work(:), SA(:,:), SB(:,:), &
         & SSv(:)
    Integer :: N, NWork, Istat

    N = Size(A,1)
    NWork = 100*Max(N**2,3*N-1)
    Allocate(Work(NWork), SA(N,N), Sb(N,N), SSV(N))

    SA = A
    Sb = B
    CALL DSYGV(1, "V", "L", N, SA, N, SB, N, SSv, Work, NWork, Istat)

    If (Present(Info)) Info = Istat
    Dv = SSv
    
    If (Present(U)) Then
       U = SA
    End If

    Deallocate(Work, SA, SSv, SB)

    Return
  End Subroutine GEVP_D


! ***********************************************
! *
  Subroutine GEVP_C(A, B, Dv, U, Info) 
! *
! ***********************************************
! * Solves the Generalized Eigenvalue Problem 
! *   Av = lBv
! * "eigenvalues" returned in Dv(:)
! *  
! ***********************************************

    Complex (kind=DPC), Intent (in) :: A(:,:), B(:,:)
    Complex (kind=DPC), Intent (out), Optional :: U(:,:)
    Real (kind=DP), Intent (out) :: Dv(:)
    Integer, Intent (out), Optional :: Info

!    Real (kind=DP) :: SA(Size(A,1),Size(A,2))
    Complex (kind=DPC), Allocatable :: Work(:), SA(:,:), SB(:,:)
    Real (kind=DP), Allocatable :: SSv(:), Rwork(:)
    Integer :: N, NWork, Istat

    N = Size(A,1)
    NWork = 100*Max(N**2,3*N-1)
    Allocate(Work(NWork), SA(N,N), Sb(N,N), SSV(N), Rwork(3*N-2))

    SA = A
    Sb = B
    CALL ZHEGV(1, "V", "L", N, SA, N, SB, N, SSv, Work, NWork, Rwork, Istat)
!        ZHEGV(I, JOB, UPL, N, A, DA, B, DB, W,   WORK, LWORK, RWORK, INFO )

    If (Present(Info)) Info = Istat
    Dv = SSv
    
    If (Present(U)) Then
       U = SA
    End If

    Deallocate(Work, SA, SSv, SB)

    Return
  End Subroutine GEVP_C

! ***********************************************
! *
  Subroutine SVD(A, Sv, U, VT, Info) 
! *
! ***********************************************

    Real (kind=DP), Intent (in) :: A(:,:)
    Real (kind=DP), Intent (out), Optional :: U(:,:), VT(:,:)
    Real (kind=DP), Intent (out) :: Sv(:)
    Integer, Intent (out), Optional :: Info

!    Real (kind=DP) :: SA(Size(A,1),Size(A,2))
    Real (kind=DP), Allocatable :: Work(:), SA(:,:), SU(:,:), &
         & SVT(:,:), SSv(:)
    Integer :: N, M, NWork, Istat

    N = Size(A,1)
    M = Size(A,2)
    NWork = 100*Max(N**2,3*N-1)
    Allocate(Work(NWork), SA(N,M), SU(N,N), SVT(M,M), SSV(Min(N,M)))

    SA = A
    CALL DGESVD("A", "A", N, M, SA, &
         & N, SSv, SU, N, SVT, M, Work, NWork, Istat)

    If (Present(U)) U = SU
    If (Present(U)) VT = SVT
    If (Present(Info)) Info = Istat
    Sv = SSv

    Deallocate(Work, SA, SSv, SU, SVT)

    Return
  End Subroutine SVD

! ***********************************************
! *
  Function PseudoInverse(A, Ikeep) Result(Ainv)
! *
! ***********************************************
! *
! * Computes the Pseudo Inverse of matrix A(:,:)
! * by keeping Ikeep singular values
! * 
! ***********************************************

    Real (kind=DP), Intent (in) :: A(:,:)
    Integer, Intent (in), Optional :: Ikeep 
    Real (kind=DP) :: Ainv(Size(A,2),Size(A,1))

    Real (kind=DP), Allocatable :: U(:,:), VT(:,:), S(:)
    Integer :: N, M, I, Ns, Ikp

    N = Size(A,1)
    M = Size(A,2)
    Ns = Min(N,M)
    Allocate(U(N,N), VT(M,M), S(Ns))

    If (Present(Ikeep)) Then
       Ikp = Ikeep
    Else
       Ikp = Size(A,1)
    End If

    CALL SVD(A, S, U, VT)
    S(1:Ikp) = 1.0_DP/S(1:Ikp)
    If (Ikp < Ns) S(Ikp+1:Ns) = 0.0_DP
    
    U = Transpose(U)
    Do I = 1, N
       U(I,:) = S(I)*U(I,:)
    End Do

    VT = Transpose(VT)
    Ainv = MatMul(VT,U)

    Deallocate(U, VT, S)

    Return
  End Function PseudoInverse

! ***********************************************
! *
  Subroutine SVD_C(A, Sv, U, VT, Info) 
! *
! ***********************************************

    Complex (kind=DPC), Intent (in) :: A(:,:)
    Complex (kind=DPC), Intent (out), Optional :: U(:,:), VT(:,:)
    Real (kind=DP), Intent (out) :: Sv(:)
    Integer, Intent (out), Optional :: Info

!    Real (kind=DP) :: SA(Size(A,1),Size(A,2))
    Complex (kind=DPC), Allocatable :: SU(:,:), SVT(:,:), &
         & SA(:,:), Work(:)
    Real (kind=DP), Allocatable :: SSv(:), Rwork(:)
    Integer :: N, M, NWork, Istat

    N = Size(A,1)
    M = Size(A,2)
    NWork = 100*Max(N**2,3*N-1)
    Allocate(Work(NWork), RWork(5*Min(M,N)), SA(N,M), &
         & SU(N,N), SVT(M,M), SSV(Min(N,M)))

    SA = A
    CALL ZGESVD("A", "A", N, M, SA, &
         & N, SSv, SU, N, SVT, M, Work, NWork, RWork, Istat)

    If (Present(U)) U = SU
    If (Present(U)) VT = SVT
    If (Present(Info)) Info = Istat
    Sv = SSv

    Deallocate(Work, SA, SSv, SU, SVT, RWork)

    Return
  End Subroutine SVD_C

! ***********************************************
! *
  Function PseudoInverse_C(A, Ikeep) Result(Ainv)
! *
! ***********************************************
! *
! * Computes the Pseudo Inverse of matrix A(:,:)
! * by keeping Ikeep singular values
! * 
! ***********************************************

    Complex (kind=DPC), Intent (in) :: A(:,:)
    Integer, Intent (in), Optional :: Ikeep 
    Complex (kind=DPC) :: Ainv(Size(A,2),Size(A,1))

    Complex (kind=DPC), Allocatable :: U(:,:), VT(:,:)
    Real (kind=DP), Allocatable :: S(:)
    Integer :: N, M, I, Ns, Ikp

    N = Size(A,1)
    M = Size(A,2)
    Ns = Min(N,M)
    Allocate(U(N,N), VT(M,M), S(Ns))

    If (Present(Ikeep)) Then
       Ikp = Ikeep
    Else
       Ikp = Size(A,1)
    End If

    CALL SVD(A, S, U, VT)
    S(1:Ikp) = 1.0_DP/S(1:Ikp)
    If (Ikp < Ns) S(Ikp+1:Ns) = 0.0_DP
    
    U = Transpose(Conjg(U))
    Do I = 1, N
       U(I,:) = S(I)*U(I,:)
    End Do

    VT = Transpose(Conjg(VT))
    Ainv = MatMul(VT,U)

    Deallocate(U, VT, S)

    Return
  End Function PseudoInverse_C

! ***********************************************
! *
  Function expM_C(A, t) Result(Aexp)
! *
! ***********************************************
! *
! * Computes the exponential of matrix A(:,:)
! * 
! ***********************************************

    Complex (kind=DPC), Intent (in) :: A(:,:)
    Real (kind=DP), Intent (in), Optional :: t
    Complex (kind=DPC) :: Aexp(Size(A,1),Size(A,2))

    Complex (kind=DPC), Allocatable :: U(:,:)
    Real (kind=DP), Allocatable :: S(:)
    Real (kind=DP) :: tf
    Integer :: N, I

    N = Size(A,1)
    Allocate(U(N,N), S(N))
    tf = 1.0_DP
    If (Present(t)) tf = t


    Aexp = A
    CALL Diag(Aexp, S, U)
    S(1:N) = exp(tf*S(1:N))
    
    Aexp = Cmplx(0.0_DP)
    Do I = 1, N
       Aexp(I,:) = S(I)*Conjg(U(:,I))
    End Do
    Aexp = MatMul(U,Aexp)

    Deallocate(U, S)

    Return
  End Function expM_C

! ***********************************************
! *
  Subroutine Diag_DPC(A, Lm, U, Info) 
! *
! ***********************************************

    Complex (kind=DPC), Intent (inout) :: A(:,:)
    Real (kind=DP), Intent (out) :: Lm(Size(A,1))
    Complex (kind=DPC), Intent (out) :: U(Size(A,1),Size(A,2))
    Integer, Intent (out), Optional :: Info

    Real (kind=DP) :: Tol
    Complex (kind=DPC) :: Work(2*Size(A,1)), RWork(7*Size(A,1))
    Integer :: Iwork(5*Size(A,1)), Ifail(Size(A,1))
    Integer :: Nfoo, Jinfo

    Tol = 2.0_DP*Tiny(0.0_DP)
    CALL ZHEEVX( 'V', 'A', 'U', Size(A,1), A, Size(A,1), 0.0_DP, &
         & 0.0_DP, 0, 0, Tol, Nfoo, Lm, U, Size(A,1), Work, &
         & 2*Size(A,1), Rwork, Iwork, IFAIL, Jinfo )

    If (Present(Info)) Info = Jinfo

    Return
  End Subroutine Diag_DPC
  
END MODULE LapackAPI
  
