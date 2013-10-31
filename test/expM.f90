
Program TestExpM

  USE NumTypes
  USE LapackAPI

  Complex (kind=DPC), Allocatable :: M(:,:), Res(:,:), P(:,:), R2(:,:)
  Real (kind=DP), Allocatable :: V(:)
  Real (kind=DP) :: fac = 1.0_DP
  Integer :: N = 8
  
  Allocate(M(N,N), V(N), P(N,N), Res(N,N), R2(N,N))

  Do I = 1, N
     CALL Random_Number(V)
     M(I,:) = Cmplx(V(:))/10.0_DP

     CALL Random_Number(V)
     M(I,:) = M(I,:) + Cmplx(0.0_DP,V(:))/10.0_DP
  End Do
!  M = Cmplx(0.0_DP)
!  Forall (I=1:N) M(I,I) = Cmplx(V(I))/10.0_DP



  M = M + Transpose(Conjg(M))
  CALL Show(M)
  R2 = expM(M)

  Res = Cmplx(0.0_DP)
  P = Cmplx(0.0_DP)
  Forall (I=1:N) P(I,I) = Cmplx(1.0_DP)
  Do I = 0, 1000
     Res = Res + P/fac
     P = MatMul(P,M)
     fac = fac*Real(I+1,kind=DP)
  End Do
  
  Res = Abs(R2-Res)
  Write(*,*)Sum(Res)
  CALL Show(Res)


  Stop
CONTAINS
  
  Subroutine Show(AA)

    Complex (kind=DPC), Intent (in) :: AA(:,:)

    Write(*,*)'[ '
    Do I = 1, Size(AA,1)
       Write(*,'(1000ES13.5)')(AA(I,J), J=1,Size(AA,2))
    End Do

    Return
  End Subroutine Show


End Program TestExpM
