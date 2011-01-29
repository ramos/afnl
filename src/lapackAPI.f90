
MODULE LapackAPI

!  USE MODMinuit
  USE NumTypes
  USE ISO_C_BINDING


  IMPLICIT NONE


CONTAINS

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
  
END MODULE LapackAPI
  
