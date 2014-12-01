
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
     Module Procedure expM_pade
  End Interface ExpM

  Private GEVP_D, GEVP_C, expM_C, expM_pade, pade3exp, pade5exp, &
       pade7exp, pade9exp, pade13exp, ludec
  

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
    
    Aexp = Cmplx(0.0_DP,kind=DPC)
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


! ***********************************************
! *
  Function expM_pade(A, t) Result(Aexp)
! *
! ***********************************************
! *
! * Computes the exponential of matrix A(:,:)
! * 
! ***********************************************

    Complex (kind=DPC), Intent (in) :: A(:,:)
    Real (kind=DP), Intent (in), Optional :: t
    Complex (kind=DPC) :: Aexp(Size(A,1),Size(A,2))

    Complex (kind=DP) :: X(Size(A,1),Size(A,2))
    Real (kind=DP) :: fac
    integer :: ns, k, ip
    Real (kind=DP), parameter :: th(5) = &
         (/1.495585217958292e-2_DP, 2.539398330063230e-1_DP, &
         9.504178996162932e-1_DP, 2.097847961257068e0_DP, &
         5.371920351148152e0_DP /), &
         log2 = 0.693147180559945309417232_DP

    if (present(t)) then
       X = t*A
    else
       X = A
    end if

    ns = size(A,1)
    fac = Sum(Abs(X(:,1)))
    do k = 2, ns
       fac = max(fac,Sum(Abs(X(:,k))))
    end do
    
    if (fac < th(1)) then
       call pade3exp(X,Aexp)
    else if (fac < th(2)) then
       call pade5exp(X,Aexp)
    else if (fac < th(3)) then
       call pade7exp(X,Aexp)
    else if (fac < th(4)) then
       call pade9exp(X,Aexp)
    else if (fac < th(5)) then
       call pade13exp(X,Aexp)
    else
       ip = int(log(fac/th(5))/log2)+1
       X = X/2.0_DP**ip

       call pade13exp(X,Aexp)
       do k = 1, ip
          Aexp = Matmul(Aexp,Aexp)
       end do
    end if

    Return
  End Function expM_pade

! *******************************************
! *
 Pure Subroutine pade3exp(X, M)
! *
! *******************************************
    
    Complex (kind=DP), Intent (in)  :: X(:,:)
    Complex (kind=DP), Intent (out) :: M(Size(X,1),Size(X,2))
    
    Real (kind=DP), parameter :: b(0:3) = &
         & (/120.0_DP, 60.0_DP, 12.0_DP, 1.0_DP/)
    Complex (kind=DP) :: Xpow(1,Size(X,1),Size(X,2)), &
         U(Size(X,1),Size(X,2)), V(Size(X,1),Size(X,2)), &
         Q(Size(X,1),Size(X,2))
    integer :: i, j, idim

    idim = Size(X,1)

    U = (0.0_DP,0.0_DP)
    V = (0.0_DP,0.0_DP)
    forall (i=1:idim) 
       U(i,i) = Cmplx(b(1),kind=DP)
       V(i,i) = Cmplx(b(0),kind=DP)
    end forall

    Xpow(1,:,:) = Matmul(X,X)
    U = Matmul(X,U + b(3)*Xpow(1,:,:))
    V = V + b(2)*Xpow(1,:,:)

    M = U + V
    Q = V - U
    CALL ludec(Q)

    Do j = 1, idim
       M(1,j) = M(1,j)
       Do I = 2, idim
          M(I,j) = M(I,j) - Sum(Q(I,1:I-1)*M(1:I-1,j))
       End Do
       
       M(idim,j) = M(idim,j) / Q(idim, idim)
       Do I = idim - 1, 1, -1
          M(I,j) = M(I,j) - Sum(Q(I, I+1:Idim)*M(I+1:Idim,j))
          M(I,j) = M(I,j) / Q(I,I)
       End Do
    End Do
    
    Return
  End Subroutine pade3exp
  

! *******************************************
! *
  Pure Subroutine pade5exp(X, M)
! *
! *******************************************
    
    Complex (kind=DP), Intent (in)  :: X(:,:)
    Complex (kind=DP), Intent (out) :: M(Size(X,1),Size(X,2))
    
    Real (kind=DP), parameter :: b(0:5) = &
         (/ 30240.0_DP, 15120.0_DP, 3360.0_DP, 420.0_DP, &
         30.0_DP, 1.0_DP /)

    Complex (kind=DP) :: Xpow(2,Size(X,1),Size(X,2)), &
         U(Size(X,1),Size(X,2)), V(Size(X,1),Size(X,2)), &
         Q(Size(X,1),Size(X,2))
    integer :: i, j, idim

    idim = Size(X,1)

    U = (0.0_DP,0.0_DP)
    V = (0.0_DP,0.0_DP)
    forall (i=1:idim) 
       U(i,i) = Cmplx(b(1),kind=DP)
       V(i,i) = Cmplx(b(0),kind=DP)
    end forall

    Xpow(1,:,:) = Matmul(X,X)
    Xpow(2,:,:) = Matmul(Xpow(1,:,:),Xpow(1,:,:))

    do i = 1,2
       V = V + b(2*i)*Xpow(i,:,:)
       U = U + b(2*i+1)*Xpow(i,:,:)
    end do
    U = Matmul(X,U)

    M = U + V
    Q = V - U
    CALL ludec(Q)

    Do j = 1, idim
       M(1,j) = M(1,j)
       Do I = 2, idim
          M(I,j) = M(I,j) - Sum(Q(I,1:I-1)*M(1:I-1,j))
       End Do
       
       M(idim,j) = M(idim,j) / Q(idim, idim)
       Do I = idim - 1, 1, -1
          M(I,j) = M(I,j) - Sum(Q(I, I+1:Idim)*M(I+1:Idim,j))
          M(I,j) = M(I,j) / Q(I,I)
       End Do
    End Do
    
    Return
  End Subroutine pade5exp
  
! *******************************************
! *
  Pure Subroutine pade7exp(X, M)
! *
! *******************************************
    
    Complex (kind=DP), Intent (in)  :: X(:,:)
    Complex (kind=DP), Intent (out) :: M(Size(X,1),Size(X,2))
    
    Real (kind=DP), parameter :: b(0:7) = &
         (/ 17297280.0_DP, 8648640.0_DP, 1995840.0_DP, 277200.0_DP, &
         25200.0_DP, 1512.0_DP, 56.0_DP, 1.0_DP/)

    Complex (kind=DP) :: Xpow(3,Size(X,1),Size(X,2)), &
         U(Size(X,1),Size(X,2)), V(Size(X,1),Size(X,2)), &
         Q(Size(X,1),Size(X,2))
    integer :: i, j, idim

    idim = Size(X,1)

    U = (0.0_DP,0.0_DP)
    V = (0.0_DP,0.0_DP)
    forall (i=1:idim) 
       U(i,i) = Cmplx(b(1),kind=DP)
       V(i,i) = Cmplx(b(0),kind=DP)
    end forall

    Xpow(1,:,:) = Matmul(X,X)
    Xpow(2,:,:) = Matmul(Xpow(1,:,:),Xpow(1,:,:))
    Xpow(3,:,:) = Matmul(Xpow(2,:,:),Xpow(1,:,:))

    do i = 1,3
       V = V + b(2*i)*Xpow(i,:,:)
       U = U + b(2*i+1)*Xpow(i,:,:)
    end do
    U = Matmul(X,U)

    M = U + V
    Q = V - U
    CALL ludec(Q)

    Do j = 1, idim
       M(1,j) = M(1,j)
       Do I = 2, idim
          M(I,j) = M(I,j) - Sum(Q(I,1:I-1)*M(1:I-1,j))
       End Do
       
       M(idim,j) = M(idim,j) / Q(idim, idim)
       Do I = idim - 1, 1, -1
          M(I,j) = M(I,j) - Sum(Q(I, I+1:Idim)*M(I+1:Idim,j))
          M(I,j) = M(I,j) / Q(I,I)
       End Do
    End Do
    
    Return
  End Subroutine pade7exp
  
! *******************************************
! *
  Pure Subroutine pade9exp(X, M)
! *
! *******************************************
    
    Complex (kind=DP), Intent (in)  :: X(:,:)
    Complex (kind=DP), Intent (out) :: M(Size(X,1),Size(X,2))
    
    Real (kind=DP), parameter :: b(0:9) = &
         (/17643225600.0_DP, 8821612800.0_DP, 2075673600.0_DP, &
         302702400.0_DP, 30270240.0_DP, 2162160.0_DP, 110880.0_DP, &
         3960.0_DP, 90.0_DP, 1.0_DP/)

    Complex (kind=DP) :: Xpow(4,Size(X,1),Size(X,2)), &
         U(Size(X,1),Size(X,2)), V(Size(X,1),Size(X,2)), &
         Q(Size(X,1),Size(X,2))
    integer :: i, j, idim

    idim = Size(X,1)

    U = (0.0_DP,0.0_DP)
    V = (0.0_DP,0.0_DP)
    forall (i=1:idim) 
       U(i,i) = Cmplx(b(1),kind=DP)
       V(i,i) = Cmplx(b(0),kind=DP)
    end forall

    Xpow(1,:,:) = Matmul(X,X)
    Xpow(2,:,:) = Matmul(Xpow(1,:,:),Xpow(1,:,:))
    Xpow(3,:,:) = Matmul(Xpow(2,:,:),Xpow(1,:,:))
    Xpow(4,:,:) = Matmul(Xpow(3,:,:),Xpow(1,:,:))

    do i = 1,4
       V = V + b(2*i)*Xpow(i,:,:)
       U = U + b(2*i+1)*Xpow(i,:,:)
    end do
    U = Matmul(X,U)

    M = U + V
    Q = V - U
    CALL ludec(Q)

    Do j = 1, idim
       M(1,j) = M(1,j)
       Do I = 2, idim
          M(I,j) = M(I,j) - Sum(Q(I,1:I-1)*M(1:I-1,j))
       End Do
       
       M(idim,j) = M(idim,j) / Q(idim, idim)
       Do I = idim - 1, 1, -1
          M(I,j) = M(I,j) - Sum(Q(I, I+1:Idim)*M(I+1:Idim,j))
          M(I,j) = M(I,j) / Q(I,I)
       End Do
    End Do
    
    Return
  End Subroutine pade9exp
  
! *******************************************
! *
  Pure Subroutine pade13exp(X, M)
! *
! *******************************************
    
    Complex (kind=DP), Intent (in)  :: X(:,:)
    Complex (kind=DP), Intent (out) :: M(Size(X,1),Size(X,2))
    
    Real (kind=DP), parameter :: b(0:13) = &
         (/ 64764752532480000.0_DP, 32382376266240000.0_DP, &
         7771770303897600.0_DP, 1187353796428800.0_DP, &
         129060195264000.0_DP, 10559470521600.0_DP, & 
         670442572800.0_DP, 33522128640.0_DP, 1323241920.0_DP, &
         40840800.0_DP, 960960.0_DP, 16380.0_DP, 182.0_DP, 1.0_DP /)

    Complex (kind=DP) :: Xpow(3,Size(X,1),Size(X,2)), &
         U(Size(X,1),Size(X,2)), V(Size(X,1),Size(X,2)), &
         Q(Size(X,1),Size(X,2))
    integer :: i, j, idim

    idim = Size(X,1)

    U = (0.0_DP,0.0_DP)
    V = (0.0_DP,0.0_DP)
    Xpow(1,:,:) = Matmul(X,X)
    do i=2, 3
       Xpow(i,:,:) = Matmul(Xpow(1,:,:),Xpow(i-1,:,:))
    end do

    do i = 1, 3
       V = V + b(2*i+6)*Xpow(i,:,:)
       U = U + b(2*i+7)*Xpow(i,:,:)
    end do
    V = Matmul(Xpow(3,:,:),V)
    U = Matmul(Xpow(3,:,:),U)
    do i = 1, 3
       V = V + b(2*i)*Xpow(i,:,:)
       U = U + b(2*i+1)*Xpow(i,:,:)
    end do
    do i = 1, idim
       V(i,i) = V(i,i) + b(0)
       U(i,i) = U(i,i) + b(1)
    end do
    U = Matmul(X,U)
    
    M = U + V
    Q = V - U
    CALL ludec(Q)

    Do j = 1, idim
       M(1,j) = M(1,j)
       Do I = 2, idim
          M(I,j) = M(I,j) - Sum(Q(I,1:I-1)*M(1:I-1,j))
       End Do
       
       M(idim,j) = M(idim,j) / Q(idim, idim)
       Do I = idim - 1, 1, -1
          M(I,j) = M(I,j) - Sum(Q(I,I+1:Idim)*M(I+1:Idim,j))
          M(I,j) = M(I,j) / Q(I,I)
       End Do
    End Do
    
    Return
  End Subroutine pade13exp
  
!  *********************************************
!  *                                           *
  Pure Subroutine ludec(M)
!  *                                           *
!  *********************************************
!  * Makes LU decomposition of Matrix M
!  *********************************************

    Complex (kind=DP), Intent (inout) :: M(:,:)

    Integer :: I, J, idim

    idim = Size(M,1)
    Do J = 2, Idim
       M(J, 1) = M(J, 1) / M(1, 1)
    End Do

    Do I = 2, Idim
       M(I, I) = M(I, I) - Sum(M(I,1:I-1)*M(1:I-1, I)) 
       Do J = i+1, Idim
          M(I, J) = M(I, J) - Sum(M(I, 1:I-1)*M(1:I-1, J))
          M(J, I) = M(J, I) - Sum(M(J, 1:I-1)*M(1:I-1, I))
          M(J, I) = M(J, I) / M(I, I)
       End Do
    End Do

    Return
  End Subroutine ludec


  
END MODULE LapackAPI
  
