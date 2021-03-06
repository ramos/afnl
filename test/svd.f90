
Program SVDtest

  USE NumTypes
  USE LapackAPI

  Real (kind=DP), Allocatable :: A(:,:), U(:,:), VT(:,:), Sv(:),&
       & SA(:,:)

  Real (kind=DP), Allocatable :: Work(:)
  Integer :: NWork
  Character (len=100) :: foo

  Write(*,*)'Enter matrix size: '
  Read(*,*)N
  M = N

  ALLOCATE(A(N,N), SA(N,N), U(N,N), VT(N,N), Sv(N))
  CALL Random_Number(A)
!  A(:,N) = A(:,1)
  Write(*,*)'A = ['
  Do K1 = 1, N
     Write(foo,*)'(1X, ', N,'ES13.5,1A1)'
     Write(*,foo)(A(K1,K2), K2=1, N), ';'
  End Do
  Write(*,*)']'

  SA = A
!  CALL DGESVD("A", "A", N, M, SA, &
!       & N, Sv, U, N, VT, M, WORK, NWORK, INFO)

  Write(*,*)Sv

  CALL SVD(A, Sv, U, VT)
  Write(*,*)Sv

  Write(*,*)'U = ['
  Do K1 = 1, N
     Write(foo,*)'(1X, ', N,'ES13.5,1A1)'
     Write(*,foo)(U(K1,K2), K2=1, N), ';'
  End Do
  Write(*,*)']'

  Write(*,*)'VT = ['
  Do K1 = 1, N
     Write(foo,*)'(1X, ', N,'ES13.5,1A1)'
     Write(*,foo)(VT(K1,K2), K2=1, N), ';'
  End Do
  Write(*,*)']'


  Do J = 1, N
     Write(*,*)'Keeping', J, ' singular values.'
     A = PseudoInverse(SA)
     U = MatMul(A,SA)
     Write(*,*)'Id = ['
     Do K1 = 1, N
        Write(foo,*)'(1X, ', N,'ES13.5,1A1)'
        Write(*,foo)(U(K1,K2), K2=1, N), ';'
     End Do
     Write(*,*)']'
  End Do


  Stop
End Program SVDtest
